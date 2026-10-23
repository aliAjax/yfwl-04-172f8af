#!/usr/bin/env python3
"""Reference solver for mixed-batch attribution and formula archiving."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import sys
import tempfile
from fractions import Fraction
from pathlib import Path
from typing import Any


TOP_KEYS = {"schema_version", "run_id", "campaign_id", "components", "samples", "batch_manifest"}
COMPONENT_KEYS = {"component_id", "formula", "spectrum", "rate_factor"}
SAMPLE_KEYS = {"sample_id", "mass_g", "spectrum", "conversion"}
BATCH_KEYS = {"batch_id", "total_mass_g", "aggregate_spectrum", "temperatures_c"}
REPORT_KEYS = {"schema_version", "run_id", "campaign_id", "assignments", "batches", "formulas", "summary"}
LEDGER_KEYS = {"run_id", "input_digest", "status", "sample_count", "batch_count"}
CHECKPOINT_KEYS = {"run_id", "input_digest", "report_digest", "batch_count"}
GAS_CONSTANT = 8.31446261815324
FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def _number(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _no_duplicate(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate object key")
        result[key] = value
    return result


def _read(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_no_duplicate)
    if not isinstance(value, dict):
        raise ValueError("top-level object required")
    return value


def _keys(obj: Any, expected: set[str], label: str) -> dict[str, Any]:
    if not isinstance(obj, dict):
        raise ValueError(f"{label} must be an object")
    if set(obj) != expected:
        raise ValueError(f"{label} key set does not match")
    return obj


def _formula_counts(formula: str) -> dict[str, int]:
    if not isinstance(formula, str) or not formula:
        raise ValueError("formula must be non-empty")
    counts: dict[str, int] = {}
    position = 0
    for match in FORMULA_TOKEN.finditer(formula):
        if match.start() != position:
            raise ValueError("formula contains an invalid token")
        element = match.group(1)
        count = int(match.group(2) or "1")
        if count <= 0:
            raise ValueError("formula count must be positive")
        counts[element] = counts.get(element, 0) + count
        position = match.end()
    if position != len(formula):
        raise ValueError("formula contains an invalid token")
    return counts


def _validate(payload: dict[str, Any]):
    _keys(payload, TOP_KEYS, "input")
    if payload["schema_version"] != 1 or isinstance(payload["schema_version"], bool):
        raise ValueError("schema_version must be integer 1")
    if not isinstance(payload["run_id"], str) or not payload["run_id"]:
        raise ValueError("run_id must be non-empty")
    if not isinstance(payload["campaign_id"], str) or not payload["campaign_id"]:
        raise ValueError("campaign_id must be non-empty")
    components = payload["components"]
    samples = payload["samples"]
    batches = payload["batch_manifest"]
    if not isinstance(components, list) or not components:
        raise ValueError("components must be a non-empty array")
    if not isinstance(samples, list) or not samples or len(samples) > 8:
        raise ValueError("samples must contain one to eight entries")
    if not isinstance(batches, list) or not batches or len(batches) > 4:
        raise ValueError("batch_manifest must contain one to four batches")

    component_ids: set[str] = set()
    component_counts: list[dict[str, int]] = []
    component_spectra: list[list[float]] = []
    for component in components:
        _keys(component, COMPONENT_KEYS, "component")
        component_id = component["component_id"]
        if not isinstance(component_id, str) or not component_id or component_id in component_ids:
            raise ValueError("component ids must be unique")
        component_ids.add(component_id)
        spectrum = component["spectrum"]
        if not isinstance(spectrum, list) or len(spectrum) != 5 or not all(_number(value) for value in spectrum):
            raise ValueError("component spectrum must contain five finite values")
        if not _number(component["rate_factor"]) or component["rate_factor"] <= 0:
            raise ValueError("component rate_factor must be positive")
        component_counts.append(_formula_counts(component["formula"]))
        component_spectra.append([float(value) for value in spectrum])

    sample_ids: set[str] = set()
    for sample in samples:
        _keys(sample, SAMPLE_KEYS, "sample")
        sample_id = sample["sample_id"]
        if not isinstance(sample_id, str) or not sample_id or sample_id in sample_ids:
            raise ValueError("sample ids must be unique")
        sample_ids.add(sample_id)
        if not _number(sample["mass_g"]) or sample["mass_g"] <= 0:
            raise ValueError("sample mass must be positive")
        spectrum = sample["spectrum"]
        if not isinstance(spectrum, list) or len(spectrum) != 5 or not all(_number(value) for value in spectrum):
            raise ValueError("sample spectrum must contain five finite values")
        conversion = sample["conversion"]
        if not isinstance(conversion, list) or len(conversion) != 4 or not all(_number(value) for value in conversion):
            raise ValueError("sample conversion must contain four finite values")
        if any(value < 0 or value > 1 for value in conversion):
            raise ValueError("sample conversion must be within [0, 1]")

    batch_ids: set[str] = set()
    for batch in batches:
        _keys(batch, BATCH_KEYS, "batch_manifest")
        batch_id = batch["batch_id"]
        if not isinstance(batch_id, str) or not batch_id or batch_id in batch_ids:
            raise ValueError("batch ids must be unique")
        batch_ids.add(batch_id)
        if not _number(batch["total_mass_g"]) or batch["total_mass_g"] <= 0:
            raise ValueError("batch total mass must be positive")
        spectrum = batch["aggregate_spectrum"]
        if not isinstance(spectrum, list) or len(spectrum) != 5 or not all(_number(value) for value in spectrum):
            raise ValueError("aggregate spectrum must contain five finite values")
        temperatures = batch["temperatures_c"]
        if not isinstance(temperatures, list) or len(temperatures) != 4 or not all(_number(value) for value in temperatures):
            raise ValueError("temperatures_c must contain four finite values")
        if any(value < 250 or value > 800 for value in temperatures):
            raise ValueError("temperatures_c is outside the supported range")
    return components, component_counts, component_spectra, samples, batches


def _weights(spectrum: list[float], component_spectra: list[list[float]]) -> list[float]:
    count = len(component_spectra)
    weights = [1.0 / count for _ in range(count)]
    observed = [float(value) for value in spectrum]
    for _ in range(360):
        for index in range(count):
            remainder = [
                sum(weights[other] * component_spectra[other][point] for other in range(count) if other != index)
                for point in range(5)
            ]
            basis = component_spectra[index]
            numerator = sum(basis[point] * (observed[point] - remainder[point]) for point in range(5))
            denominator = sum(value * value for value in basis)
            if denominator > 1e-12:
                weights[index] = max(0.0, numerator / denominator)
        total = sum(weights)
        weights = [value / total for value in weights]
    return [round(value, 9) for value in weights]


def _partition(samples: list[dict[str, Any]], batches: list[dict[str, Any]]) -> dict[str, str]:
    order = sorted(range(len(samples)), key=lambda index: -float(samples[index]["mass_g"]))
    assigned: dict[str, str] = {}
    best: tuple[float, dict[str, str]] | None = None
    sums = {batch["batch_id"]: 0.0 for batch in batches}
    vectors = {batch["batch_id"]: [0.0] * 5 for batch in batches}
    batch_ids = [batch["batch_id"] for batch in batches]

    def search(position: int) -> None:
        nonlocal best
        if position == len(order):
            error = 0.0
            for batch in batches:
                total = sums[batch["batch_id"]]
                if abs(total - float(batch["total_mass_g"])) > 1e-6:
                    return
                if total <= 0:
                    return
                vector = [value / total for value in vectors[batch["batch_id"]]]
                error += sum(
                    (vector[index] - float(batch["aggregate_spectrum"][index])) ** 2
                    for index in range(5)
                )
            signature = dict(sorted(assigned.items()))
            candidate = (error, signature)
            if best is None or candidate < best:
                best = candidate
            return
        sample = samples[order[position]]
        mass = float(sample["mass_g"])
        for batch in batches:
            batch_id = batch["batch_id"]
            if sums[batch_id] + mass > float(batch["total_mass_g"]) + 1e-6:
                continue
            if assigned and best is not None and sums[batch_id] + mass > float(batch["total_mass_g"]) + 1e-6:
                continue
            assigned[sample["sample_id"]] = batch_id
            sums[batch_id] += mass
            for index in range(5):
                vectors[batch_id][index] += mass * float(sample["spectrum"][index])
            search(position + 1)
            for index in range(5):
                vectors[batch_id][index] -= mass * float(sample["spectrum"][index])
            sums[batch_id] -= mass
            del assigned[sample["sample_id"]]

    search(0)
    if best is None or best[0] > 1e-8:
        raise ValueError("batch mass and aggregate spectrum do not match a partition")
    return best[1]


def _fit_batch(
    samples: list[dict[str, Any]],
    assigned: dict[str, str],
    batch: dict[str, Any],
    component_spectra: list[list[float]],
    components: list[dict[str, Any]],
) -> dict[str, Any]:
    xs: list[float] = []
    ys: list[float] = []
    for sample in samples:
        if assigned[sample["sample_id"]] != batch["batch_id"]:
            continue
        weights = _weights(sample["spectrum"], component_spectra)
        factor = sum(weights[index] * float(components[index]["rate_factor"]) for index in range(len(weights)))
        for temperature, conversion in zip(batch["temperatures_c"], sample["conversion"]):
            if 1e-8 < conversion < 1.0 - 1e-8:
                rate = -math.log(1.0 - conversion) / 600.0 / factor
                if rate > 0:
                    ys.append(math.log(rate))
                    xs.append(1.0 / float(temperature))
    if len(xs) < 2:
        raise ValueError("not enough conversion points to fit Arrhenius parameters")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denominator = sum((value - mean_x) ** 2 for value in xs)
    if denominator <= 1e-18:
        raise ValueError("temperature program has no variation")
    slope = sum((xs[index] - mean_x) * (ys[index] - mean_y) for index in range(len(xs))) / denominator
    intercept = mean_y - slope * mean_x
    energy = -slope * GAS_CONSTANT
    predictions = [intercept + slope * value for value in xs]
    rmse = math.sqrt(sum((ys[index] - predictions[index]) ** 2 for index in range(len(xs))) / len(xs))
    return {
        "batch_id": batch["batch_id"],
        "ln_a": round(intercept, 6),
        "activation_energy_j_mol": round(energy, 6),
        "rmse": round(rmse, 9),
    }


def _canonical_multiplier(counts: list[int]) -> list[int]:
    # 配方规范化：在元素计数集合上按固定取整规则归一化。
    positive = [value for value in counts if value > 0]
    if not positive:
        raise ValueError("formula has no atoms")
    divisor = 0
    for value in positive:
        divisor = math.gcd(divisor, value)
    return [value // max(1, divisor) for value in counts]


def _formula(
    batch_samples: list[dict[str, Any]],
    component_counts: list[dict[str, int]],
    component_spectra: list[list[float]],
) -> str:
    if not batch_samples:
        raise ValueError("cannot normalize an empty batch formula")
    totals: dict[str, Fraction] = {}
    for sample in batch_samples:
        weights = _weights(sample["spectrum"], component_spectra)
        for index, counts in enumerate(component_counts):
            for element, value in counts.items():
                totals[element] = totals.get(element, Fraction(0, 1)) + Fraction(
                    str(weights[index])
                ) * value
    ordered_elements = sorted(totals, key=lambda element: (element != "C", element != "H", element))
    integer_counts = [
        int(round(float(totals[element]) * 1000)) for element in ordered_elements
    ]
    normalized = _canonical_multiplier(integer_counts)
    pieces = []
    for element, count in zip(ordered_elements, normalized):
        pieces.append(element if count == 1 else f"{element}{count}")
    return "".join(pieces)


def _canonical(report: dict[str, Any]) -> bytes:
    return (json.dumps(report, ensure_ascii=True, separators=(",", ":")) + "\n").encode("utf-8")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _append(path: Path, row: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, ensure_ascii=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _build(payload: dict[str, Any]) -> tuple[dict[str, Any], str]:
    components, component_counts, component_spectra, samples, batches = _validate(payload)
    assigned = _partition(samples, batches)
    assignments = []
    for sample in samples:
        assignments.append(
            {
                "sample_id": sample["sample_id"],
                "batch_id": assigned[sample["sample_id"]],
                "weights": _weights(sample["spectrum"], component_spectra),
            }
        )
    assignments.sort(key=lambda item: item["sample_id"])
    fitted = [
        _fit_batch(samples, assigned, batch, component_spectra, components)
        for batch in batches
    ]
    fitted.sort(key=lambda item: item["batch_id"])
    formulas = []
    for batch in batches:
        batch_samples = [
            sample for sample in samples if assigned[sample["sample_id"]] == batch["batch_id"]
        ]
        formulas.append(
            {
                "batch_id": batch["batch_id"],
                "formula": _formula(batch_samples, component_counts, component_spectra),
            }
        )
    formulas.sort(key=lambda item: item["batch_id"])
    report = {
        "schema_version": 1,
        "run_id": payload["run_id"],
        "campaign_id": payload["campaign_id"],
        "assignments": assignments,
        "batches": fitted,
        "formulas": formulas,
        "summary": {
            "sample_count": len(samples),
            "batch_count": len(batches),
            "converged": True,
        },
    }
    if set(report) != REPORT_KEYS:
        raise ValueError("internal report shape mismatch")
    return report, hashlib.sha256(_canonical(report)).hexdigest()


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("input refused: expected input and output paths", file=sys.stderr)
        return 2
    input_path = Path(argv[1])
    output_dir = Path(argv[2])
    try:
        payload = _read(input_path)
        report, report_digest = _build(payload)
        output_dir.mkdir(parents=True, exist_ok=True)
        input_digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
        ledger_row = {
            "run_id": payload["run_id"],
            "input_digest": input_digest,
            "status": "complete",
            "sample_count": len(payload["samples"]),
            "batch_count": len(payload["batch_manifest"]),
        }
        _keys(ledger_row, LEDGER_KEYS, "ledger")
        checkpoint = {
            "run_id": payload["run_id"],
            "input_digest": input_digest,
            "report_digest": report_digest,
            "batch_count": len(payload["batch_manifest"]),
        }
        _keys(checkpoint, CHECKPOINT_KEYS, "checkpoint")
        ledger_path = output_dir / "ledger.jsonl"
        previous = ledger_path.read_text(encoding="utf-8").splitlines() if ledger_path.exists() else []
        row_text = json.dumps(ledger_row, ensure_ascii=True, separators=(",", ":"))
        if not previous or previous[-1] != row_text:
            _append(ledger_path, ledger_row)
        _atomic_write(output_dir / "report.json", _canonical(report))
        _atomic_write(
            output_dir / "checkpoint.json",
            (json.dumps(checkpoint, ensure_ascii=True, separators=(",", ":")) + "\n").encode("utf-8"),
        )
        return 0
    except (OSError, ValueError, TypeError, KeyError, ZeroDivisionError) as exc:
        print(f"input refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
