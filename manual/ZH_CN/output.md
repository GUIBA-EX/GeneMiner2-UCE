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

**assembly_graphs**：仅在使用 `--assembler-graph-format gfa`、`dot` 或 `both` 时生成的逐 locus 压缩组装图；默认不生成。

## ITS2 专用输出

**its2_assembly_summary.csv**：ITS2 模式的逐 locus 汇总表，记录组装状态、接受状态、代表序列长度、reads 支持、深度指标和候选数。

**results/<locus>.fasta**：该 locus 所有已接受 ITS2 变体，而不只是单条最佳序列。标题包含 fragment、paired-fragment、diagnostic-fragment、EM support 和 abundance。

**results/<locus>.its2_support.tsv**：已接受和低支持候选的明细表，包含 equivalence members、总 fragment support、paired support、diagnostic support、EM support、EM abundance 和 `PASS`、`EQUIVALENCE_GROUP` 或 `LOW_SUPPORT` 状态。

**contigs_all/<locus>.fasta**：所有已接受 ITS2 变体；`contigs_all_low/<locus>.fasta` 保存未达到 fragment support 标准的候选。

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

## Population 输出

`population` 子命令在已有输出目录下创建 `population/`：

- `sample_manifest.tsv`：原始样本名、GeneMiner2 内部目录名、VCF 样本名、reads 路径和 SE/PE 布局的对应关系。
- `reference/population_reference.fasta`：所有样本统一 mapping 使用的公共 UCE 参考。使用 `--population-reference-fasta` 时为复制后的外部参考。
- `reference/population_reference_provenance.tsv`：内部构建参考时，每个 locus 的来源样本、选择策略、候选数、长度及 reads 支持指标。
- `reference/reference_contribution.tsv`：内部构建参考时，每个样本贡献的参考 loci 数及比例，用于检查公共参考是否由少数样本主导。
- `reference/locus_name_map.tsv`：内部构建参考时原始 locus 名和 VCF 安全名称的对应关系。
- `reference/reference_source.tsv`：使用固定外部参考时，记录源文件与复制后的参考路径。
- `mapping/<sample>.bam` 和索引：minibwa 统一 mapping 后经 samtools 处理的 BAM。
- `mapping/mapping_qc.tsv`：每个样本的 mapped/properly-paired reads、mapping rate、覆盖广度和平均深度。
- `variants/cohort.raw.bcf`：联合检测的原始 cohort BCF。
- `variants/cohort.biallelic.snps.vcf.gz`：拆分并保留的双等位 SNP。
- `variants/cohort.genotype_filtered.vcf.gz`：低 DP/GQ 基因型设为缺失后的 VCF。
- `variants/cohort.tagged.vcf.gz`：补充群体统计标签后的 VCF。
- `variants/cohort.filtered.vcf.gz`：经过 QUAL、call rate 和 MAC 过滤的分析起点。
- `variants/variant_qc.tsv`：raw calling、双等位、基因型过滤、标签补充和位点过滤各阶段的位点数与文件路径。
- `structure/all_snps.vcf.gz` 和同名前缀 PLINK/PCA 文件：保留同一 UCE 内多个 SNP 的高敏感度面板。
- `structure/one_snp_per_uce.vcf.gz`、`population.{bed,bim,fam}` 和 `population_pca.*`：每个 UCE 按 call rate、QUAL 和 MAC 选择一个代表 SNP 的主面板。
- `structure/ld_pruned.vcf.gz` 和同名前缀 PLINK/PCA 文件：对全部 SNP 做 LD pruning 的敏感性面板。
- `structure/selected_snps.tsv`：每个 UCE 被选中代表 SNP 的位置和筛选统计。
- `structure/panel_summary.tsv`：三种面板的 SNP 数量、路径和建议用途。
- `structure/admixture/K<K>.log`、`population.<K>.Q` 和 `population.<K>.P`：各 K 的 ADMIXTURE 日志、个体遗传成分和祖先群体等位基因频率。
- `structure/admixture/cv_errors.tsv`：各 K 的交叉验证误差以及最低误差 K 标记。
- `structure/admixture/status.tsv`：ADMIXTURE 的 `complete`、`skipped`、`unavailable` 或 `failed` 状态。

PCA 和 ADMIXTURE 的主解释应优先使用每个 UCE 一个 SNP 的面板，并与 all-SNP 和 LD-pruned PCA 比较。若样本 mapping rate、coverage breadth 明显偏低，应先排查缺失数据；仅在内部构建参考时，才根据 `reference_contribution.tsv` 的参考来源不均衡排查参考偏倚。

## 统计输出

`stats` 子命令会从已有输出目录生成类似 HybPiper 的 UCE 恢复统计。

**uce_stats.tsv**：样本级统计表，汇总过滤 reads、成功 loci、低质量 loci、主要失败类型、按参考长度比例统计的恢复数量、rescue 回退数量、总恢复碱基数、contig 长度、reads 支持跨度和 read density。

**uce_locus_stats.tsv**：locus 级统计表，汇总参考长度、样本占有率、恢复长度、reads 支持、侧翼平衡度和候选 contig 数量。

**uce_seq_lengths.tsv**：sample-by-locus contig 长度矩阵，包含根据参考目录计算的 `MeanLength` 行。

**uce_read_counts.tsv**：sample-by-locus 矩阵，记录最佳 contig 的 read count。

**uce_filtered_read_counts.tsv**：sample-by-locus 矩阵，根据 `ref_reads_count_dict.txt` 记录第一轮过滤阶段分配到各 locus 的 reads 数量。

**uce_rescue_stats.tsv**：由已有 `uce_rescue_summary.csv` 转换得到的 rescue 明细表。

**uce_recovery_heatmap.png** 和 **uce_read_counts_heatmap.png**：`stats` 子命令生成的热图；如果环境中没有 `pandas`、`seaborn` 或 `matplotlib`，或使用了 `--stats-no-heatmap`，则不会生成。
