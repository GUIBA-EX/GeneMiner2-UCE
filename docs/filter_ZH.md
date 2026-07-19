# 1. Filter

[English version](filter_EN.md)

Filter 是其余三条流程共用的 reads 招募层：它根据参考 k-mer 扫描 reads，并写出分配至各 locus 的 reads；它本身不决定 contig 是否接受、marker 丰度或群体结构。

## 输入与参考格式

样本表为 tab 分隔：样本名、R1 和可选 R2。`original` 与 `uce` 的 `-r` 应为目录，其中每个 locus 一个 FASTA，文件主名即 locus 名。profiling 是例外：它直接接收一个 `.fa` 或 `.fasta` marker 参考库。

## 首轮招募与 refilter

`filter` 以 `-kf` 和扫描步长 `-s` 招募参考相关 reads，临时结果写入 `filtered_pe/`。`refilter` 随后细化逐 locus 分配并施加深度和文件大小限制，保留 reads 写入 `filtered/`。

`uce` 模式下，只要一对 reads 中任一 mate 通过，整对都会保留，使 core 命中的 mate 可保留另一端的 flank 证据。profiling 只执行首轮招募，不运行 `refilter`。

## 缓存与可复现性

重复使用未改变的参考时，加 `--reuse-reference-cache`。缓存指纹包含参考文件、过滤 k-mer 和扫描步长。`--reference-cache-dir` 可指定共享目录；默认位于 `output/.gm2_reference_cache`。

## 输出与检查

- `ref_reads_count_dict.txt`：首轮逐 locus read 数。
- `filtered_pe/`：首轮分配结果；refilter 成功后通常删除。
- `filtered/`：传递给组装器的 reads。
- `large_files/`：只有超过 refilter 深度或大小限制时才出现。

某 locus 未招募到 reads 本身不能直接证明不存在；应检查参考距离、`-kf`、测序深度和 read-count 表。文件字段见[输出说明](../manual/ZH_CN/output.md)，参数定义见[命令行指南](../manual/ZH_CN/command_line.md)。
