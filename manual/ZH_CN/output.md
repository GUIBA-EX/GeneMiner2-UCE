# 输出文件

本文档说明当前命令行版本会生成的输出文件。旧 GUI、细胞器组装和 Windows 辅助脚本相关输出已经不属于本 CLI fork 的支持范围。

## 单样本目录

输入样本表中的每个样本都会在输出目录下生成一个同名文件夹。

**filtered_pe**：第一轮参考 k-mer 过滤后分配到各 locus 的临时 paired-end reads。如果同一次运行继续完成 re-filter，该目录通常会在 re-filter 成功后被删除。

**filtered**：进一步过滤后保留的 reads。UCE 模式下，只要 paired-end 的任一端通过 locus 过滤，整对 reads 都会被保留。

**large_files**：进一步过滤时超过深度限制或文件大小限制的 reads。只有产生这类文件时才会出现。

**results**：主组装结果。每个 locus 的最佳 contig 写为 `<locus>.fasta`。

**contigs_all**：组装器评估过的候选 contigs。

**contigs_all_low**：UCE 模式下保留的低支持延伸候选 contigs，仅用于检查，不会提升到主 `results` 目录。

**consensus**：可选 `consensus` 命令生成的一致性序列。

**blast**：可选 `trim` 命令生成的参考切齐结果。

**log.txt**：单样本日志文件。

**result_dict.txt**：每个 locus 的组装状态和 reads 支持摘要。

**ref_reads_count_dict.txt**：第一轮过滤阶段每个 locus 分配到的 reads 数量。

**uce_assembly_summary.csv**：UCE 模式下的单样本逐 locus 组装摘要，包含状态、最佳 contig 长度、reads 支持跨度、read count、read density、支持比例、侧翼平衡度、k-mer 深度指标、候选数和低质量标记。

**uce_rescue_summary.csv**：使用 `--assembly-mode uce --uce-rescue-reads` 时生成的单样本 raw-read rescue 摘要。

## 合并输出

**combined_results**：按 locus 合并不同样本恢复序列后的文件。

**combined_results/aligned**：`combine` 阶段生成的多序列比对结果。

**combined_trimed**：比对列过滤后的逐 locus alignment。默认由 trimAl 生成；使用 `--alignment-filter alifilter` 时由 AliFilter 生成；使用 `--alignment-filter none` 或 `--no-trimal` 时不会生成。

**combined_results.fasta**：比对列过滤前的串联 alignment。

**combined_trimed.fasta**：比对列过滤后的串联 alignment。

**failed_gene_trees.tsv**：溯祖树流程中单基因树构建失败的 locus 列表。

**failed_samples.tsv**：当任一样本在 filter、refilter、assemble 或 UCE rescue 阶段失败时生成。命令行流程会在写出该文件后停止，避免下游步骤在样本结果不完整时继续运行。

**Coalescent.tree**：溯祖流程生成的物种树。

**Concatenation.tree**：串联流程生成的系统树。

## UCE 专用输出

**uce_contigs**：UCE 组装模式生成的 phyluce 兼容 contig 输出。每个样本一个 `*.contigs.fasta` 文件。`sample_name_map.tsv` 记录 GeneMiner2 样本名与 phyluce 安全样本名的对应关系。

**uce_rescue_summary.csv**：跨样本合并后的 rescue 摘要。使用 `--assembly-mode uce --uce-rescue-reads` 时生成。

## 统计输出

`stats` 子命令会从已有输出目录生成类似 HybPiper 的 UCE 恢复统计。

**uce_stats.tsv**：样本级统计表，汇总过滤 reads、成功 loci、低质量 loci、主要失败类型、按参考长度比例统计的恢复数量、rescue 回退数量、总恢复碱基数、contig 长度、reads 支持跨度和 read density。

**uce_locus_stats.tsv**：locus 级统计表，汇总参考长度、样本占有率、恢复长度、reads 支持、侧翼平衡度和候选 contig 数量。

**uce_seq_lengths.tsv**：sample-by-locus contig 长度矩阵，包含根据参考目录计算的 `MeanLength` 行。

**uce_read_counts.tsv**：sample-by-locus 矩阵，记录最佳 contig 的 read count。

**uce_filtered_read_counts.tsv**：sample-by-locus 矩阵，根据 `ref_reads_count_dict.txt` 记录第一轮过滤阶段分配到各 locus 的 reads 数量。

**uce_rescue_stats.tsv**：由已有 `uce_rescue_summary.csv` 转换得到的 rescue 明细表。

**uce_recovery_heatmap.png** 和 **uce_read_counts_heatmap.png**：`stats` 子命令生成的热图；如果环境中没有 `pandas`、`seaborn` 或 `matplotlib`，或使用了 `--stats-no-heatmap`，则不会生成。
