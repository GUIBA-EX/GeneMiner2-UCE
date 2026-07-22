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

`original` 默认使用确定性的 `original-rust`；UCE 固定使用 Rust `ucefilter` 与 `uce-rust`，不会回退到 Python。

## 安装

需要 Rust/Cargo、Python 3.11、Cython、PyInstaller 与生物信息学依赖。请按[命令行指南](manual/ZH_CN/command_line.md)创建环境；随后在仓库根目录运行：

```bash
make build
```

入口程序为 `cli/geneminer2`。更新源码后请再次运行 `make build`。

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
| 性能设计与边界 | [MainFilter 性能说明](docs/development/mainfilter-performance.md) |

加 `--workflow-profile` 可在输出根目录写入 `workflow_profile.tsv`；支持的 Rust assembler 还会写入样本级 `assembly_profile.tsv`。两者只记录耗时与 I/O，不改变结果。

## 引用与支持

请引用：Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

问题和功能建议请提交 [GitHub Issue](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)。项目以 [GPL-3.0-or-later](LICENSE) 发布；第三方与移植代码的来源边界见 [NOTICE](NOTICE)。
