# GeneMiner2-UCE

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[English](README_EN.md)** · [更新记录](CHANGELOG.md) · [报告问题](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

GeneMiner2-UCE 是面向 genome skimming、target capture 和相关短读长数据的参考引导恢复工具。它以 Rust 为默认计算后端，支持 UCE、常规 marker、线粒体、marker profiling、UCE 群体遗传、核基因家族和 repeatome 工作流。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 选择工作流

| 目标 | 命令或模式 | 主要结果 |
| --- | --- | --- |
| 恢复 exon、SCO 或一般核/线粒体 marker | `--assembly-mode original`（默认） | 参考引导 contig，默认按参考裁切 |
| 恢复 UCE core 与有证据的侧翼 | `--assembly-mode uce` | UCE contig、组装与 rescue 汇总 |
| 只评估 marker 的 reads 支持 | `profiling` | 每条参考序列的支持证据 |
| 组装常规单环动物线粒体 | `mito` | 受 reads 支持的 closed 或 partial 结果 |
| 对多样本 UCE 做变异分析 | `population` | 公共伪参考、联合 VCF、PCA、ADMIXTURE 输入 |
| 恢复核基因家族候选 | `gene` | family 状态、候选 FASTA 与队列汇总 |
| 从短 reads 建立 repeatome | `te` | 保守 repeat library、注释与样本 RPM |
| 用 WGS 补充 ipyrad RAD loci | `rad` / `rad-validate` | 双 arm 恢复、状态表与验证后的严格矩阵 |

`original` 默认使用确定性的 `original-rust`；UCE 固定使用 Rust `ucefilter` 与 `uce-rust`，不会回退到 Python。

## 安装

需要 Rust/Cargo、Python 3.11、Cython、PyInstaller 与生物信息学依赖。请按[命令行指南](manual/ZH_CN/command_line.md)创建环境；随后在仓库根目录运行：

```bash
make build
```

入口程序为 `cli/geneminer2`。更新源码后请再次运行 `make build`。
该入口是原生 Rust 调度器；Python 实现仅保留为开发期的兼容性对照，不作为运行时回退。

## 快速开始：UCE

样本表为 tab 分隔的 `sample_id  R1  [R2]`；UCE 需要配对 reads。参考目录中每个 locus 一个 FASTA，文件名即 locus 名。

```text
samples.tsv
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz

references/
  uce-0001.fasta
  uce-0002.fasta
```

```bash
cli/geneminer2 filter assemble \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

### UCE 的两条 reads 路径

默认 UCE 不走通用的 `MainFilter + refilter`（常简称 **main + re**）两步路径。`ucefilter` 把宽松招募、paired-fragment 保留、方向/精确匹配证据和逐 locus 的深度/位置选择融合为一次原始 FASTQ 扫描：

```text
原始 paired FASTQ
  → ucefilter（招募 + 选择；融合 main + re 的职责）
  → <sample>/filtered/
  → uce-rust 主组装
  → contigs + uce_assembly_summary.csv
```

因此默认 UCE 使用 `filter assemble`，没有独立的 `refilter` 步骤，也不会写出 `filtered_pe/`。`--legacy-uce-filter` 才恢复兼容性路径：`MainFilter → filtered_pe/ → refilter → filtered/ → uce-rust`；它仅用于对照或诊断，不是推荐的正式路线。

`--uce-rescue-reads` 在主组装接受 contig 后才执行：第一轮以“原参考 + 已接受 contig”重新招募并重组装；第二轮仅对仍在增长的 locus 使用 contig 两端窗口。每轮按 locus 审计增量，任一侧或整条 contig 证据不足都会回退，因此 rescue 只延长已成立的主结果，不以参考填补缺口。

首先检查：

- `output/uce_assembly_summary.csv`：每个 locus 的组装选择与 QC；
- `output/uce_rescue_summary.csv`：可选 rescue 的增量与回退；
- `output/uce_contigs/`：最终 UCE contig。

没有 `--uce-rescue-reads` 时只执行主 UCE 恢复；该开关最多进行两轮受约束延伸。

## 其他常用入口

```bash
# 常规 marker；不写 --assembly-mode 也是 original
cli/geneminer2 filter assemble -f samples.tsv -r references -o original_out -p 8 \
  --assembly-mode original

# 核基因家族
cli/geneminer2 gene -f samples.tsv -r family_references -o gene_out -p 8

# 常规单环动物线粒体；需要带注释 GenBank 参考
cli/geneminer2 mito -f samples.tsv -o mito_out -p 8 \
  --mito-genbank mitochondrial_reference.gb

# marker reads 证据，不组装
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o profile_out -p 8

# 已完成 UCE 组装后的群体流程
cli/geneminer2 population -f samples.tsv -r references -o population_out -p 8 \
  --assembly-mode uce --engine panrefv2

# 参考无关 repeatome；使用独立的 TE 样本表
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

`mito` 仅适用于常规单条环形动物线粒体；不应用于多分子、显著重排或复杂植物/真菌线粒体。`profiling` 提供参考相容性证据，不等同于唯一物种鉴定或丰度估计。

## RAD：用 WGS 补充已有 RAD 矩阵

RAD 路径仅补充 **ipyrad `.loci` 中缺失的新 WGS 样本**。R1 与 R2 始终是两条独立的限制性位点 arm；流程不会推断它们之间未测序的插入区，也不会把 WGS 恢复结果当作 restriction-site allele dropout 的直接证据。

最稳妥的路线是：先由已完成的 `.loci` 建 probe，恢复 WGS arm，再独立验证。

```bash
cli/geneminer2 rad-probe --ipyrad-loci assembly.loci -o rad_probe
cli/geneminer2 rad --rad-probe rad_probe/rad_reference \
  -f wgs_samples.tsv -o rad_out -p 8
cli/geneminer2 rad-validate --rad-probe rad_probe/rad_reference \
  --rad-recovery rad_out/rad_recovery -o rad_validate_out
```

`rad-probe --ipyrad-params params.txt` 可调度 ipyrad 1--7 步；`--rad-denovo` 仅适用于已经 demultiplex 的 paired RAD reads，生成保守的 candidate probe，不替代 ipyrad 的完整 clustering。`rad-validate` 只接纳 R1、R2 都通过覆盖度、identity 与跨 locus 竞争检查的样本。完整输入、参数、输出和解释见 [RAD 文档](docs/rad_CN.md)。

## 文档

| 需要做什么 | 文档 |
| --- | --- |
| 安装、全部参数与外部依赖 | [命令行指南](manual/ZH_CN/command_line.md) |
| 理解输出目录和表格 | [输出文件说明](manual/ZH_CN/output.md) |
| reads 过滤与缓存 | [Filter](docs/filter_ZH.md) |
| original/UCE 组装与 rescue | [Assembler](docs/assembler_ZH.md) |
| marker reads profiling | [Profiling](docs/profiling_ZH.md) |
| UCE 群体遗传 | [Population](docs/population_ZH.md) |
| 核基因家族 | [Gene](docs/gene_ZH.md) |
| repeatome | [TE / repeatome](docs/te_ZH.md) |
| 线粒体 | [Mito](docs/mitochondria_CN.md) |
| WGS 补充 RAD loci | [RAD](docs/rad_CN.md) |
| 性能设计与边界 | [MainFilter 性能说明](docs/development/mainfilter-performance.md) |

加 `--workflow-profile` 可在输出根目录写入原子更新的 `workflow_profile.tsv`：其中包括参考索引、逐样本 filter/refilter/assemble/rescue，以及 cohort 与后续汇总阶段；失败阶段也会保留 `failed` 状态。支持的 Rust assembler 还会写入样本级 `assembly_profile.tsv`。两者只记录耗时与 I/O，不改变结果。

`combine` 中，`-p` 是并发 locus 数，`--msa-threads` 是每个 MSA 进程的线程数（不得大于 `-p`），`--filter-processes` 独立限制同时运行的 trimAl/AliFilter 数量。三个值均须至少为 1。

`--cleanup-intermediates` 是显式的完成后清理选项：仅在同一次调用完成 filter 与 assemble，且所有后续步骤成功后，删除可再生的 filtered reads 与 rescue 临时参考；最终 contig、汇总、原始 reads 和参考均保留。删除记录写入 `cleanup_manifest.tsv`。

## 引用与支持

请引用：Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

问题和功能建议请提交 [GitHub Issue](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)。项目以 [GPL-3.0-or-later](LICENSE) 发布；第三方与移植代码的来源边界见 [NOTICE](NOTICE)。
