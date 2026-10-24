"""Behavior checks for mixed-batch attribution and formula archiving."""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path


ENTRY = os.environ.get("YEL004_ENTRY", "/app/mixed_batch.py")
R = 8.31446261815324


def _components() -> list[dict]:
    return [
        {
            "component_id": "C1",
            "formula": "C2H4O2",
            "spectrum": [1.0, 0.2, 0.1, 0.0, 0.0],
            "rate_factor": 1.0,
        },
        {
            "component_id": "C2",
            "formula": "C3H6O3",
            "spectrum": [0.1, 1.0, 0.2, 0.1, 0.0],
            "rate_factor": 1.2,
        },
        {
            "component_id": "C3",
            "formula": "CH2O",
            "spectrum": [0.2, 0.3, 1.0, 0.4, 0.1],
            "rate_factor": 0.8,
        },
    ]


def _case(sample_count: int = 6, batch_count: int = 2, run_id: str = "MIX-RUN-1") -> dict:
    components = _components()
    weights_list = [
        [0.50, 0.30, 0.20],
        [0.40, 0.40, 0.20],
        [0.60, 0.25, 0.15],
        [0.30, 0.50, 0.20],
        [0.45, 0.35, 0.20],
        [0.25, 0.45, 0.30],
        [0.55, 0.20, 0.25],
        [0.35, 0.30, 0.35],
    ]
    masses = [10.0, 20.0, 15.0, 12.0, 18.0, 14.0, 16.0, 11.0]
    temperature_sets = [
        [380.0, 400.0, 420.0, 440.0],
        [390.0, 410.0, 430.0, 450.0],
        [370.0, 395.0, 415.0, 435.0],
        [385.0, 405.0, 425.0, 445.0],
    ]
    parameters = [
        (3.0, 28000.0),
        (3.4, 32000.0),
        (3.2, 30000.0),
        (3.6, 34000.0),
    ]
    groups = [[] for _ in range(batch_count)]
    for index in range(sample_count):
        groups[index % batch_count].append(index)
    samples = []
    for index in range(sample_count):
        weights = weights_list[index]
        spectrum = [
            sum(weights[component] * components[component]["spectrum"][point] for component in range(3))
            for point in range(5)
        ]
        batch_index = index % batch_count
        ln_a, energy = parameters[batch_index]
        factor = sum(weights[component] * components[component]["rate_factor"] for component in range(3))
        conversion = []
        for temperature in temperature_sets[batch_index]:
            rate = math.exp(ln_a - energy / (R * temperature)) * factor
            conversion.append(1.0 - math.exp(-rate * 600.0))
        samples.append(
            {
                "sample_id": f"S-{index + 1:02d}",
                "mass_g": masses[index],
                "spectrum": spectrum,
                "conversion": conversion,
            }
        )
    manifest = []
    for batch_index, group in enumerate(groups):
        total = sum(samples[index]["mass_g"] for index in group)
        aggregate = [
            sum(samples[index]["mass_g"] * samples[index]["spectrum"][point] for index in group) / total
            for point in range(5)
        ]
        manifest.append(
            {
                "batch_id": f"B-{batch_index + 1:02d}",
                "total_mass_g": total,
                "aggregate_spectrum": aggregate,
                "temperatures_c": temperature_sets[batch_index],
            }
        )
    return {
        "schema_version": 1,
        "run_id": run_id,
        "campaign_id": "CAMP-001",
        "components": components,
        "samples": samples,
        "batch_manifest": manifest,
    }


def _run(case: dict, tmp_path: Path):
    input_path = tmp_path / "campaign.json"
    output_dir = tmp_path / "result"
    input_path.write_text(json.dumps(case), encoding="utf-8")
    output_dir.mkdir(exist_ok=True)
    proc = subprocess.run(
        [sys.executable, ENTRY, str(input_path), str(output_dir)],
        capture_output=True,
        text=True,
        timeout=45,
    )
    report_path = output_dir / "report.json"
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
    return proc, report, output_dir


def _run_raw(raw: str, tmp_path: Path):
    input_path = tmp_path / "campaign.json"
    output_dir = tmp_path / "result"
    input_path.write_text(raw, encoding="utf-8")
    output_dir.mkdir(exist_ok=True)
    return subprocess.run(
        [sys.executable, ENTRY, str(input_path), str(output_dir)],
        capture_output=True,
        text=True,
        timeout=45,
    )


def test_output_contract_creates_three_artifacts(tmp_path):
    """题面：「归档完成后目录内包含配方报告、追加账本与恢复记录。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert proc.stdout == ""
    assert proc.stderr == ""
    assert (output_dir / "ledger.jsonl").is_file()
    assert (output_dir / "checkpoint.json").is_file()
    assert report["summary"]["batch_count"] == 2


def test_success_output_channels_are_silent(tmp_path):
    """题面：「归档完成时终端没有文本输出。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert proc.stdout == ""
    assert proc.stderr == ""


def test_top_level_keys_reject_missing_field(tmp_path):
    """题面：「混样归档 JSON 最外层字段为 schema_version,run_id,campaign_id,components,samples,batch_manifest。」"""
    case = _case()
    del case["campaign_id"]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_component_keys_reject_extra_field(tmp_path):
    """题面：「成分定义携带 component_id,formula,spectrum,rate_factor。」"""
    case = _case()
    case["components"][0]["note"] = "extra"
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_sample_keys_reject_missing_field(tmp_path):
    """题面：「样品记录携带 sample_id,mass_g,spectrum,conversion。」"""
    case = _case()
    del case["samples"][0]["conversion"]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_batch_keys_reject_extra_field(tmp_path):
    """题面：「批次汇总携带 batch_id,total_mass_g,aggregate_spectrum,temperatures_c。」"""
    case = _case()
    case["batch_manifest"][0]["note"] = "extra"
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_spectrum_length_accepts_five_values(tmp_path):
    """题面：「spectrum 与 aggregate_spectrum 各有 5 个数值。」"""
    case = _case()
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 0
    assert len(case["samples"][0]["spectrum"]) == 5
    assert report["summary"]["sample_count"] == 6


def test_spectrum_length_rejects_four_values(tmp_path):
    """题面：「spectrum 与 aggregate_spectrum 各有 5 个数值。」"""
    case = _case()
    case["samples"][0]["spectrum"] = case["samples"][0]["spectrum"][:-1]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_conversion_length_accepts_four_values(tmp_path):
    """题面：「conversion 与 temperatures_c 各有 4 个数值。」"""
    case = _case()
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 0
    assert len(case["samples"][0]["conversion"]) == 4
    assert report["summary"]["sample_count"] == 6


def test_duplicate_identifiers_are_refused(tmp_path):
    """题面：「component_id,sample_id 与 batch_id 不得重复。」"""
    case = _case()
    case["samples"][1]["sample_id"] = case["samples"][0]["sample_id"]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_boolean_schema_version_is_refused(tmp_path):
    """题面：「质量、谱图和温度字段都需有限数值，布尔值不是有效数值。」"""
    case = _case()
    case["schema_version"] = True
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert not (output_dir / "report.json").exists()


def test_maximum_sample_count_is_supported(tmp_path):
    """题面：「托盘只接收八件样品，批次上限四组。」"""
    proc, report, output_dir = _run(_case(sample_count=8, batch_count=2), tmp_path)
    assert proc.returncode == 0
    assert report["summary"]["sample_count"] == 8


def test_maximum_batch_count_is_supported(tmp_path):
    """题面：「托盘只接收八件样品，批次上限四组。」"""
    proc, report, output_dir = _run(_case(sample_count=8, batch_count=4), tmp_path)
    assert proc.returncode == 0
    assert report["summary"]["batch_count"] == 4


def test_aggregate_mass_must_match_partition(tmp_path):
    """题面：「batch_manifest 的质量和聚合谱图必须与实际样品分组一致。」"""
    case = _case()
    case["batch_manifest"][0]["total_mass_g"] += 1.0
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert "input refused:" in proc.stderr
    assert not (output_dir / "report.json").exists()


def test_report_bytes_end_with_newline(tmp_path):
    """题面：「文件末尾只有一个换行，正文按 UTF-8 紧凑写出。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    raw = (output_dir / "report.json").read_bytes()
    assert proc.returncode == 0
    assert raw.endswith(b"\n")
    assert b'": ' not in raw


def test_report_key_set_is_exact(tmp_path):
    """题面：「归档报告在 schema_version 与 run_id,campaign_id 占据归档报告前三项，余下四项是 assignments,batches,formulas,summary。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert set(report) == {
        "schema_version",
        "run_id",
        "campaign_id",
        "assignments",
        "batches",
        "formulas",
        "summary",
    }


def test_assignment_item_key_set_is_exact(tmp_path):
    """题面：「归属行保留 sample_id,batch_id,weights。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert all(set(item) == {"sample_id", "batch_id", "weights"} for item in report["assignments"])


def test_batch_item_key_set_is_exact(tmp_path):
    """题面：「批次参数行保留 batch_id,ln_a,activation_energy_j_mol,rmse。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert all(
        set(item) == {"batch_id", "ln_a", "activation_energy_j_mol", "rmse"}
        for item in report["batches"]
    )


def test_formula_item_key_set_is_exact(tmp_path):
    """题面：「配方行保留 batch_id,formula。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert all(set(item) == {"batch_id", "formula"} for item in report["formulas"])


def test_summary_key_set_is_exact(tmp_path):
    """题面：「归档报告尾部三个字段为 sample_count,batch_count,converged。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert set(report["summary"]) == {"sample_count", "batch_count", "converged"}


def test_assignment_recovers_hidden_batch_membership(tmp_path):
    """题面：「归属步骤比较质量与混合谱图，从而还原样品所属批次。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assigned = {item["sample_id"]: item["batch_id"] for item in report["assignments"]}
    assert proc.returncode == 0
    assert assigned == {
        "S-01": "B-01",
        "S-02": "B-02",
        "S-03": "B-01",
        "S-04": "B-02",
        "S-05": "B-01",
        "S-06": "B-02",
    }


def test_component_weights_are_nonnegative_and_normalized(tmp_path):
    """题面：「weights 是成分顺序固定的非负数组，所有分量的和为 1。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    assert proc.returncode == 0
    assert all(all(value >= 0 for value in item["weights"]) for item in report["assignments"])
    assert all(abs(sum(item["weights"]) - 1.0) < 1e-6 for item in report["assignments"])


def test_arrhenius_parameters_are_recovered_for_each_batch(tmp_path):
    """题面：「同一批次内的转化率与温度程序联合拟合 Arrhenius 截距和活化能。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    fitted = {item["batch_id"]: item for item in report["batches"]}
    assert proc.returncode == 0
    assert abs(fitted["B-01"]["ln_a"] - 3.0) < 0.01
    assert abs(fitted["B-01"]["activation_energy_j_mol"] - 28000.0) < 20.0
    assert abs(fitted["B-02"]["ln_a"] - 3.4) < 0.01
    assert abs(fitted["B-02"]["activation_energy_j_mol"] - 32000.0) < 20.0


def _published_conversion_deltas(case: dict, report: dict) -> list[float]:
    assignments = {item["sample_id"]: item for item in report["assignments"]}
    fitted = {item["batch_id"]: item for item in report["batches"]}
    components = case["components"]
    deltas = []
    for sample in case["samples"]:
        assignment = assignments[sample["sample_id"]]
        batch = fitted[assignment["batch_id"]]
        factor = sum(
            assignment["weights"][index] * components[index]["rate_factor"]
            for index in range(len(components))
        )
        batch_case = next(
            item
            for item in case["batch_manifest"]
            if item["batch_id"] == assignment["batch_id"]
        )
        for temperature, observed in zip(
            batch_case["temperatures_c"], sample["conversion"]
        ):
            rate = math.exp(
                batch["ln_a"]
                - batch["activation_energy_j_mol"] / (R * temperature)
            ) * factor
            predicted = 1.0 - math.exp(-rate * 600.0)
            deltas.append(abs(predicted - observed))
    return deltas


def test_published_rate_law_reproduces_fixture_conversions(tmp_path):
    """题面：「`rate_factor` 不是待拟合参数，也不能被忽略或再次缩放。」"""
    case = _case()
    proc, report, output_dir = _run(case, tmp_path)
    deltas = _published_conversion_deltas(case, report)
    assert proc.returncode == 0
    assert deltas
    assert max(deltas) < 2e-5


def test_formula_uses_hill_order_and_compact_counts(tmp_path):
    """题面：「formula 使用 Hill 顺序，先写 C,再写 H，其余元素按字母顺序。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    formulas = {item["batch_id"]: item["formula"] for item in report["formulas"]}
    assert proc.returncode == 0
    assert formulas["B-01"].startswith("C")
    assert "H" in formulas["B-01"]
    assert re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formulas["B-01"])


def test_repeated_run_keeps_report_bytes(tmp_path):
    """题面：「同一批归档重放时，报告文件字节不变。」"""
    proc_first, report_first, output_dir = _run(_case(), tmp_path)
    first = (output_dir / "report.json").read_bytes()
    proc_second, report_second, _ = _run(_case(), tmp_path)
    assert proc_first.returncode == 0
    assert proc_second.returncode == 0
    assert (output_dir / "report.json").read_bytes() == first


def test_repeated_run_does_not_append_ledger(tmp_path):
    """题面：「重复导入后，账本仍保持原来的行数。」"""
    proc_first, report_first, output_dir = _run(_case(), tmp_path)
    first = (output_dir / "ledger.jsonl").read_bytes()
    proc_second, report_second, _ = _run(_case(), tmp_path)
    assert proc_first.returncode == 0
    assert proc_second.returncode == 0
    assert (output_dir / "ledger.jsonl").read_bytes() == first


def test_changed_digest_appends_ledger(tmp_path):
    """题面：「归档内容一变，账本就多一条批次交接。」"""
    proc_first, report_first, output_dir = _run(_case(), tmp_path)
    proc_second, report_second, _ = _run(_case(run_id="MIX-RUN-2"), tmp_path)
    lines = (output_dir / "ledger.jsonl").read_text(encoding="utf-8").splitlines()
    assert proc_first.returncode == 0
    assert proc_second.returncode == 0
    assert len(lines) == 2


def test_ledger_line_has_exact_keys(tmp_path):
    """题面：「sample_count 与 batch_count 在混样账本末尾，前面依次是 run_id,input_digest,status。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    row = json.loads((output_dir / "ledger.jsonl").read_text(encoding="utf-8").splitlines()[0])
    assert proc.returncode == 0
    assert set(row) == {"run_id", "input_digest", "status", "sample_count", "batch_count"}
    assert len(row["input_digest"]) == 64


def test_checkpoint_has_exact_keys(tmp_path):
    """题面：「恢复对象开头是 batch_count，随后是 run_id、input_digest、report_digest。」"""
    proc, report, output_dir = _run(_case(), tmp_path)
    checkpoint = json.loads((output_dir / "checkpoint.json").read_text(encoding="utf-8"))
    assert proc.returncode == 0
    assert set(checkpoint) == {"run_id", "input_digest", "report_digest", "batch_count"}
    assert len(checkpoint["report_digest"]) == 64


def test_missing_key_is_refused(tmp_path):
    """题面：「缺键、多余键、重复键都会让混样分群停止。」"""
    case = _case()
    del case["samples"]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert "input refused:" in proc.stderr
    assert not (output_dir / "report.json").exists()


def test_duplicate_json_key_is_refused(tmp_path):
    """题面：「缺键、多余键、重复键都会让混样分群停止。」"""
    raw = json.dumps(_case())[:-1] + ',"run_id":"MIX-RUN-2"}'
    proc = _run_raw(raw, tmp_path)
    assert proc.returncode == 2
    assert "input refused:" in proc.stderr
    assert not (tmp_path / "result" / "report.json").exists()


def test_refusal_uses_fixed_stderr_prefix(tmp_path):
    """题面：「拒收原因只在标准错误出现一次，首字符序列为 input refused:。」"""
    case = _case()
    case["batch_manifest"][0]["temperatures_c"] = [100.0, 200.0, 300.0, 400.0]
    proc, report, output_dir = _run(case, tmp_path)
    assert proc.returncode == 2
    assert proc.stderr.splitlines()[0].startswith("input refused:")
    assert not (output_dir / "report.json").exists()


def test_refusal_preserves_existing_artifacts(tmp_path):
    """题面：「拒收事件发生后，既有的报告文件、账本文件和检查点文件均不变化。」"""
    proc_first, report_first, output_dir = _run(_case(), tmp_path)
    old_report = (output_dir / "report.json").read_bytes()
    old_ledger = (output_dir / "ledger.jsonl").read_bytes()
    old_checkpoint = (output_dir / "checkpoint.json").read_bytes()
    broken = _case()
    del broken["components"]
    proc_second, report_second, _ = _run(broken, tmp_path)
    assert proc_first.returncode == 0
    assert proc_second.returncode == 2
    assert (output_dir / "report.json").read_bytes() == old_report
    assert (output_dir / "ledger.jsonl").read_bytes() == old_ledger
    assert (output_dir / "checkpoint.json").read_bytes() == old_checkpoint
