# Population 流程说明

[English version](population_EN.md)

`population` 面向多个二倍体 UCE 样本，从已接受 contig 和原始 reads 构建公共伪参考、联合 SNP、PCA 与 ADMIXTURE 面板。它输出未定相基因型，不代替需要完整单倍型或单 locus 基因树的 phasing 分析。

参数全集见[中文命令行指南](../manual/ZH_CN/command_line.md)，文件列表见[输出文件说明](../manual/ZH_CN/output.md)。

## 输入要求

每个样本必须已经完成 UCE 组装，并保留：

- 样本目录中的 `uce_assembly_summary.csv`；
- `results/` 中已接受的 UCE contig；
- 样本表所列的原始 R1 和可选 R2。

`population` 只使用通过 UCE 接受标准且具有唯一 reads 支持的 contig。所有样本随后统一 mapping 到同一公共参考，避免分别对各自组装结果检测变异。

## 四个阶段

```text
reference → mapping → calling → selection
```

1. `reference`：构建或复制公共伪参考；
2. `mapping`：minibwa mapping，samtools 排序、去重复并统计 QC；
3. `calling`：bcftools 联合检测和过滤 SNP；
4. `selection`：构建三种 SNP 面板，运行 PLINK PCA 和可选 ADMIXTURE。

使用 `--population-stop-after` 可在指定阶段后停止。`--population-start-at mapping|calling|selection` 只复用经过检查的既有输出；所需参考、BAM 或过滤 VCF 缺失时程序会报错，不会静默混用文件。

## 公共伪参考

默认 `--population-reference-strategy sqcl-longest` 对每个 locus 选择最长的合格 contig，类似 SqCL `make_PRG.py` 的代表序列策略。`supported` 策略优先 reads 支持，再考虑长度和其他指标。

内部构建时会输出：

- `population_reference.fasta`：统一 mapping 参考；
- `population_reference_provenance.tsv`：每个 locus 的来源样本、候选数和支持指标；
- `reference_contribution.tsv`：每个样本贡献的参考 loci 数和比例；
- `locus_name_map.tsv`：原始 locus 名与 VCF 安全名称的对应关系。

使用固定参考验证伪参考或重复分析时，可指定：

```bash
--population-reference-fasta cohort_reference.fasta
```

外部参考会复制到 Population 输出目录，并在 `reference_source.tsv` 中记录来源。内部伪参考是从一个真实已接受 contig 中选择代表序列，不是把多个个体强行拼成共识，因此应结合贡献比例和 mapping QC 检查参考偏倚。

## Mapping 和变异检测

默认外部工具为 minibwa、samtools 和 bcftools。调用阶段默认使用：

```text
minimum MAPQ = 20
minimum base quality = 20
minimum genotype DP = 5
minimum genotype GQ = 20
minimum site QUAL = 20
minimum call rate = 0.8
minimum minor allele count = 2
```

低 DP/GQ 的基因型先设为缺失，再在位点层面应用 QUAL、call rate 和 MAC。阈值应按测序深度和缺失模式调整，不宜只根据最终 SNP 数量放宽。

## SNP 面板

流程生成三种互补面板：

- `all_snps`：保留同一 UCE 内多个 SNP，信息量高，但连锁和长 locus 贡献会更强；
- `one_snp_per_uce`：每个 UCE 选择一个代表 SNP，是 PCA 和 ADMIXTURE 的默认主面板；
- `ld_pruned`：从全部 SNP 中按 LD pruning 保留相对独立的位点，用于敏感性比较。

PLINK 对三种面板分别运行 PCA。ADMIXTURE 默认使用 `one_snp_per_uce`，在 `--population-admixture-k-min` 到 `--population-admixture-k-max` 范围内计算交叉验证误差。缺少 ADMIXTURE 时，伪参考、VCF、PLINK 和 PCA 仍会完成，并在状态文件中记录 `unavailable`。

## 必查 QC

正式解释 PCA 或遗传成分前至少检查：

- `mapping/mapping_qc.tsv`：各样本 mapping rate、coverage breadth 和 depth；
- `variants/variant_qc.tsv`：每个过滤阶段保留的变异数；
- `reference/reference_contribution.tsv`：内部参考是否由少数样本贡献过多；
- `structure/panel_summary.tsv`：三个面板的 SNP 数量；
- 样本和位点缺失率；
- 三种 SNP 面板 PCA 的主要结构是否一致；
- ADMIXTURE 不同 K 的 CV error 是否稳定。

低 mapping 或低覆盖样本应先排查参考距离、文库质量、测序深度和污染。只有内部构建公共参考时，才根据 `reference_contribution.tsv` 判断参考来源偏倚。

## 推荐验证伪参考

验证 Population 伪参考时，建议同时运行：

1. 默认 `sqcl-longest`；
2. `supported`；
3. 固定外部参考（如果存在可靠参考基因组或已验证 UCE reference）。

比较三套结果的 mapping rate、coverage breadth、SNP 缺失率、PCA 主轴和每个 locus 的变异数。若主要群体结构只在某一种伪参考下出现，应优先视为参考敏感性信号，而不是直接作生物学解释。
