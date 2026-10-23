# 多批混样的配方与动力学归档

材料实验组把几批混样放在同一个接收托盘里，样品标签与批次汇总在扫描时被打散。每批留下了总质量,混合谱图和一组温度程序，样品本身保留质量,谱图和转化率曲线。实验组要在投料评审前恢复样品归属,联合得到每批 Arrhenius 参数，并把配方整理成可以直接比较的规范式。

mixed_batch.py 读取归档 JSON，并把结果写入给定目录。

归档完成时终端没有文本输出。

混样归档 JSON 最外层字段为 schema_version,run_id,campaign_id,components,samples,batch_manifest。

成分定义携带 component_id,formula,spectrum,rate_factor。

样品记录携带 sample_id,mass_g,spectrum,conversion。

批次汇总携带 batch_id,total_mass_g,aggregate_spectrum,temperatures_c。

spectrum 与 aggregate_spectrum 各有 5 个数值。

conversion 与 temperatures_c 各有 4 个数值。

component_id,sample_id 与 batch_id 不得重复。

质量、谱图和温度字段都需有限数值，布尔值不是有效数值。

托盘只接收八件样品，批次上限四组。

归档完成后目录内包含配方报告、追加账本与恢复记录。

文件末尾只有一个换行，正文按 UTF-8 紧凑写出。

归档报告在 schema_version 与 run_id,campaign_id 占据归档报告前三项，余下四项是 assignments,batches,formulas,summary。

归属行保留 sample_id,batch_id,weights。

批次参数行保留 batch_id,ln_a,activation_energy_j_mol,rmse。

配方行保留 batch_id,formula。

归档报告尾部三个字段为 sample_count,batch_count,converged。

归属步骤比较质量与混合谱图，从而还原样品所属批次。

weights 是成分顺序固定的非负数组，所有分量的和为 1。

同一批次内的转化率与温度程序联合拟合 Arrhenius 截距和活化能。

formula 使用 Hill 顺序，先写 C,再写 H，其余元素按字母顺序。

formula 是元素计数的规范化结果，计数为 1 时省略数字。

同一批归档重放时，报告文件字节不变。

重复导入后，账本仍保持原来的行数。

归档内容一变，账本就多一条批次交接。

sample_count 与 batch_count 在混样账本末尾，前面依次是 run_id,input_digest,status。

恢复对象开头是 batch_count，随后是 run_id、input_digest、report_digest。

缺键、多余键、重复键都会让混样分群停止。

拒收原因只在标准错误出现一次，首字符序列为 input refused:。

拒收事件发生后，既有的报告文件、账本文件和检查点文件均不变化。

batch_manifest 的质量和聚合谱图必须与实际样品分组一致。

交付结果应让下一班同时核对配方与动力学参数。
