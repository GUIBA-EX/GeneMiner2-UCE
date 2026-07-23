# TStools (old GeneMiner2-UCE)

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[English](README_EN.md)** · [更新记录](CHANGELOG.md) · [报告问题](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

TStools 是一个 Rust 原生的短读长恢复工具：以参考序列招募 reads，再按任务组装或量化。它面向 genome skimming、target capture、UCE、RAD 补充和 repeatome 分析；不需要 Python 运行时。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 安装与最小示例

安装 Rust/Cargo 与所需生物信息学软件后，在仓库根目录运行：

```bash
cargo run -p xtask -- build
```

统一入口是 `cli/geneminer2`。样本表使用 tab 分隔的 `sample_id  R1  [R2]`；参考目录中每个 FASTA 文件代表一个 locus 或 bait。

```bash
# UCE：原始 paired reads → 选择性招募 → UCE 组装
cli/geneminer2 filter assemble \
  -f samples.tsv -r uce_references -o uce_out -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

首先查看 `uce_out/uce_assembly_summary.csv` 和 `uce_out/uce_contigs/`。`--uce-rescue-reads` 是可选的、最多两轮的受证据约束延伸；它不会用参考序列虚构缺口。

## 选择命令

| 目标 | 命令 | 主要结果 |
| --- | --- | --- |
| 常规 exon、SCO 或核 marker | `filter assemble` | 参考引导 contig |
| UCE core 与 reads 支持的侧翼 | `filter assemble --assembly-mode uce` | UCE contig 与 QC 汇总 |
| 线粒体 | `mito` | closed、linear 或明确的结构歧义判定 |
| marker 的 reads 支持 | `profiling` | 每条参考的支持证据 |
| UCE 群体遗传 | `population` | 公共伪参考、VCF、PCA 等 |
| 核基因家族 | `gene` | family 候选、copy 状态与后续解析输入 |
| WGS 补充既有 RAD 矩阵 | `rad-probe` → `rad` → `rad-validate` | 双 arm 恢复与严格矩阵 |
| 无参考 repeatome | `te` | 保守 repeat library、注释和 RPM |

默认的常规组装使用确定性的 `original-rust`；UCE 使用专用的 `ucefilter` 与 `uce-rust`。两条路线都是 Rust 实现。

## 常用命令

```bash
# 常规 marker；original 为默认模式
cli/geneminer2 filter assemble -f samples.tsv -r references -o marker_out -p 8

# 单环动物线粒体；需要带注释的 GenBank 参考
cli/geneminer2 mito -f samples.tsv -o mito_out -p 8 \
  --mito-genbank mitochondrial_reference.gb

# 核基因家族
cli/geneminer2 gene -f samples.tsv -r family_baits -o gene_out -p 8

# 不组装，只评估 reads 对 marker 的支持
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o profile_out -p 8

# UCE 组装完成后的群体流程
cli/geneminer2 population -f samples.tsv -r uce_references -o population_out -p 8 \
  --assembly-mode uce --engine panrefv2

# 无参考 repeatome；使用独立样本表
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

## 结果的边界

- UCE 默认将广泛招募与每 locus 的 reads 选择合并为一次 FASTQ 扫描；`--legacy-uce-filter` 仅用于对照与诊断。
- `mito` 只适用于常规单环动物线粒体。串联重复或超过 insert size 的完全重复不能由短 reads 可靠定拷贝数，结果会保留为 linear 或 ambiguous。
- `profiling` 是参考相容性证据，不是物种鉴定或丰度估计。
- RAD 中 R1/R2 是独立限制性位点 arm；WGS 恢复不直接证明 allele dropout。请以 `rad-validate` 的双 arm 检查为准。
- `--cleanup-intermediates` 只在同次完整流程成功后删除可再生的过滤 reads；最终 contig、汇总、原始 reads 和参考始终保留。

## 文档

| 主题 | 文档 |
| --- | --- |
| 安装、参数、外部依赖 | [命令行指南](manual/ZH_CN/command_line.md) |
| 输出目录与结果表 | [输出说明](manual/ZH_CN/output.md) |
| 过滤与缓存 | [Filter](docs/filter_ZH.md) |
| 常规/UCE 组装 | [Assembler](docs/assembler_ZH.md) |
| 线粒体、Gene、RAD、TE | [Mito](docs/mitochondria_CN.md) · [Gene](docs/gene_ZH.md) · [RAD](docs/rad_CN.md) · [TE](docs/te_ZH.md) |
| 群体与 profiling | [Population](docs/population_ZH.md) · [Profiling](docs/profiling_ZH.md) |

`--workflow-profile` 会写入仅记录时间与 I/O 的 `workflow_profile.tsv`，不会改变分析结果。

## 引用

请引用：Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

```bibtex
@software{GeneMiner2_UCE,
  author  = {XIA, Fei and TANG, Zizhen and XU, Yan},
  title   = {GeneMiner2-UCE: Reference-Guided Short-Read Recovery for UCE, Mitochondrial, Gene-Family, and RAD Workflows},
  year    = {2026},
  version = {1.5.6},
  url     = {https://github.com/GUIBA-EX/GeneMiner2-UCE},
  publisher = {GitHub},
  note    = {GPL-3.0-or-later licensed software}
}
```

项目以 [GPL-3.0-or-later](LICENSE) 发布；第三方与移植代码的来源边界见 [NOTICE](NOTICE)。
