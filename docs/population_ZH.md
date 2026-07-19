# 4. Population

前置章节：[Filter](filter_ZH.md)、[Assembler](assembler_ZH.md)。参数与完整文件清单见[命令行指南](../manual/ZH_CN/command_line.md)和[输出说明](../manual/ZH_CN/output.md)。

`population` 将多个已完成的二倍体 UCE 样本统一到一套参考坐标，输出联合 SNP、PCA 和可选 ADMIXTURE。它产生未定相基因型，不替代单倍型或单 locus 基因树分析。

## 输入与运行

每个样本需要保留：`uce_assembly_summary.csv`、`results/`中的接受 contig，以及样本表中的原始 R1/R2。运行依赖 minibwa、samtools、bcftools 和 PLINK 1.9；ADMIXTURE 缺失时，其余结果仍会完成。

推荐以 PanRefV2 构建公共 UCE 参考：

```bash
cli/geneminer2 population \
  -f samples.tsv -r baits_by_locus -o output -p 32 \
  --assembly-mode uce --engine panrefv2
```

`baits_by_locus`中每个 FASTA 文件代表一个 locus。PanRefV2 不接受 `--population-reference-fasta`；若已有固定参考，使用 `--engine pseudoref --population-reference-fasta FILE`。

## 参考引擎

| 引擎 | 用途 |
| --- | --- |
| `pseudoref` | 默认；从已接受 contig 选择每个 locus 的代表序列。可用 `--population-reference-strategy sqcl-longest|supported`。 |
| `panref` | 旧实验实现，仅为兼容保留。 |
| `panrefv2` | 流式、局部图参考构建；推荐用于本项目的 UCE Population 流程。 |

PanRefV2使用唯一 bait minimizer 将read pair归入locus：先以双端共同命中建立冻结核心图，再只接收连接核心图的单端救援read。它不写按locus拆分的FASTQ，不跨UCE locus连接序列。

每条backbone按以下顺序选择：

```text
bait一致性 > 支持样本数 > 当前图边的read/PE支持 > 每样本深度稳定性 > 长度
```

## 阶段与结果

```text
reference → mapping → calling → selection
```

- `reference`：构建或复制公共参考；
- `mapping`：minibwa 与 samtools 统一比对、去重复和QC；
- `calling`：bcftools 联合检测与过滤；
- `selection`：生成 all-SNP、每UCE一个SNP和LD-pruned面板，运行PCA与可选ADMIXTURE。

用 `--population-start-at` 和 `--population-stop-after` 重启或截停阶段；后续阶段只复用已检查的已有输出。

PanRefV2额外输出到 `population/reference/panrefv2/`：

- `index_metadata.tsv`：bait minimizer索引摘要；
- `recruitment_summary.tsv`：强证据、候选救援、实际救援和歧义pair；
- `population_graph.gfa`：按locus分区的局部DBG；
- `locus_summary.tsv`：状态、长度、unitig数、样本支持和backbone路径证据。

默认只有`pass` locus写入`population_reference.fasta`。`short`、`low_sample_support`、`low_coverage`、`complex`和`no_core`保留在统计中但不参与后续比对。仅作诊断时可加入：

```bash
--population-panrefv2-include-low-confidence
```

## 默认过滤与QC

默认阈值：MAPQ 20、base quality 20、genotype DP 5、GQ 20、site QUAL 20、call rate 0.8、MAC 2。低DP/GQ基因型先设为缺失，再执行位点过滤。

解释群体结果前至少检查：

- `reference/panrefv2/locus_summary.tsv`：恢复数、短locus和复杂图；
- `mapping/mapping_qc.tsv`：mapping rate、breadth和depth；
- `variants/variant_qc.tsv`：各过滤阶段SNP数；
- `structure/panel_summary.tsv`：三套面板的位点数；
- 三套面板PCA是否给出一致的主要结构；以及ADMIXTURE不同K的CV误差。

主解释通常使用每UCE一个SNP的面板，并以all-SNP和LD-pruned结果作敏感性检查。
