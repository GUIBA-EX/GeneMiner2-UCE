# GeneMiner2-UCE

**[English README](README_EN.md)**

GeneMiner2-UCE 是 GeneMiner2 的 UCE 拓展版本，面向 target-enrichment、ultraconserved element（UCE）以及相关二代测序数据。它保留参考引导的 reads 捕获框架，并扩展了 UCE 侧翼恢复、ITS2 多候选组装和 UCE 群体遗传分析。

使用本软件时请引用 [GeneMiner2-UCE GitHub 仓库](https://github.com/GUIBA-EX/GeneMiner2-UCE)；正式论文发布后将在此补充。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 主要功能

- 从 genome skimming 或 target-capture reads 中恢复目标分子标记。
- 在 UCE 模式中保留有 reads 支持的 core 与侧翼序列。
- 在 ITS2 模式中保留多个候选，并报告 paired、diagnostic 和 EM 丰度证据。
- 从多个 UCE 样本构建公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。
- 导出 PHYLUCE 兼容 contig，并汇总样本和 locus 层面的恢复质量。

## 选择运行模式

| 模式 | 适用数据 | 主要结果 |
| --- | --- | --- |
| `--assembly-mode reference` | genome skimming、常规基因恢复 | 参考引导 contig；默认流程包含参考裁切 |
| `--assembly-mode uce` | UCE target capture | UCE core 与有 reads 支持的 flanking sequence |
| `--assembly-mode its2` | 同一 locus 可能存在多个真实 ITS2 变体 | 多候选序列及 fragment-level 支持 |
| `population` 子命令 | 二倍体 UCE 群体重测序或 target capture | 公共伪参考、联合 SNP、PCA 和 ADMIXTURE 面板 |

## 安装

GeneMiner2-UCE 目前从源码构建。在仓库根目录运行：

```bash
make
```

构建后的入口为：

```bash
cli/geneminer2
```

每次拉取包含源码变更的新版本后应重新运行 `make`。完整构建依赖和外部程序要求见[中文命令行指南](manual/ZH_CN/command_line.md)。

## 快速开始

准备 tab 分隔的样本表；每行依次为样本名、R1 和可选的 R2：

```text
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz
Sample_B\t/path/to/Sample_B_R1.fq.gz
```

参考目录中每个 locus 使用一个独立 FASTA 文件，文件主名作为 locus 名：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

运行 UCE 恢复并执行一轮受控的 raw-read rescue：

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --uce-rescue-reads
```

首先检查以下结果：

- `output/uce_assembly_summary.csv`：每个样本和 locus 的组装状态与支持指标；
- `output/uce_rescue_summary.csv`：rescue 前后比较和回退状态；
- `output/uce_contigs/`：PHYLUCE 兼容的样本 contig。

完整目录结构见[输出文件说明](manual/ZH_CN/output.md)。

## UCE 模式

UCE 模式降低短 probe 边界对组装的限制，默认跳过参考引导的 `trim`，并优先保留延伸更长且仍有 reads 支持的候选。refilter 阶段只要任一 mate 通过 locus 过滤，整对 paired-end reads 都会保留。

默认 Rust assembler 使用不反复 backtrack 的 backbone 策略；`--uce-rescue-reads` 以第一轮 contig 和原始参考再次招募 reads，并在 rescue 质量下降时恢复第一轮结果。参数、质量护栏、reference cache 和回退规则见 [UCE 流程说明](docs/uce-workflow_ZH.md)。

## ITS2 模式

ITS2 模式固定使用 Rust assembler 和 21-mer，允许每个 locus 输出多个候选，并报告 fragment、paired-fragment、diagnostic-fragment 和 EM abundance。无法由 reads 区分的候选保留为 equivalence group。

常规 `combine` 面向每 locus 单序列，因此 ITS2 应显式运行：

```bash
cli/geneminer2 filter refilter assemble \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode its2
```

ITS2 不支持原始 Python assembler，也不会在 Rust assembler 失败时回退。结果位于各样本的 `results/`、`its2_assembly_summary.csv` 和 `<locus>.its2_support.tsv`。

## Population 模式

`population` 从已接受的 UCE contig 和原始 reads 构建一致的未定相二倍体 SNP 矩阵。它先生成或复用公共伪参考，再统一 mapping、联合检测变异，并输出 all-SNP、每个 UCE 一个 SNP 和 LD-pruned 三种面板。

已有 UCE 组装结果后运行：

```bash
cli/geneminer2 population \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

运行时需要 minibwa、samtools、bcftools 和 PLINK 1.9；ADMIXTURE 为可选依赖。伪参考策略、断点续跑、SNP 面板和必查 QC 见 [Population 流程说明](docs/population_ZH.md)。

## 实现与文档

默认构建 Rust MainFilter、Refilter、Assembler、Population 和 Rust 辅助工具。`--assembler-implementation auto` 会优先使用 Rust；失败时回退到未修改的原始 Python assembler，但 ITS2 模式除外。主 CLI 编排器和 consensus 程序继续使用 Python。

- [中文命令行指南](manual/ZH_CN/command_line.md)
- [中文输出文件说明](manual/ZH_CN/output.md)
- [UCE 流程说明](docs/uce-workflow_ZH.md)
- [Population 流程说明](docs/population_ZH.md)
- [MainFilter 性能与兼容性说明](docs/mainfilter-performance.md)
- [版本更新记录](CHANGELOG.md)

## 引用与联系

GeneMiner2 主要引用：

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

相关工具包括 [Easy353](https://doi.org/10.1093/molbev/msac261) 和 [GeneMiner](https://doi.org/10.1111/1755-0998.13924)。使用 `--alignment-filter alifilter` 时还应引用 [AliFilter](https://doi.org/10.1093/molbev/msag097)。

如需联系项目维护者，请发送邮件至 [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp)。
