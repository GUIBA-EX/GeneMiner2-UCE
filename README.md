# GeneMiner2-UCE

**[English README](README_EN.md)**

GeneMiner2-UCE 是 GeneMiner2 的 UCE 扩展版，面向 target-enrichment、UCE 和相关短读长数据。它提供参考引导 marker 恢复、UCE 侧翼恢复、marker profiling、常规动物线粒体组装和 UCE 群体遗传分析。

使用本软件请引用 [GeneMiner2-UCE GitHub 仓库](https://github.com/GUIBA-EX/GeneMiner2-UCE)。

![GeneMiner2-UCE 流程](docs/images/summary_ZH.png)

## 功能

- 从 genome skimming 或 target-capture reads 恢复 marker。
- 恢复 UCE core 及有 reads 支持的侧翼。
- 从 WGS 或宏基因组 reads 评估参考序列支持。
- 构建公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。
- 从带注释 GenBank 参考恢复常规环形动物线粒体。

## 模式选择

| 模式 | 适合的数据 | 主要结果 |
| --- | --- | --- |
| `--assembly-mode original` | exon、SCO、核或线粒体 marker | 参考引导 contig；默认按参考裁切 |
| `--assembly-mode uce` | genome skimming 或 target capture 的 UCE | UCE core 及有 reads 支持的侧翼 |
| `profiling` | WGS 或宏基因组中的任意 marker | 每条参考序列的支持证据 |
| `mito` | 常规环形动物线粒体 | 有 reads 支持的闭环或 partial 组装 |
| `population` | 已完成 UCE 组装的多个样本 | 公共伪参考、VCF、PCA 与 ADMIXTURE 输入 |

## 安装

在含 Cython、PyInstaller 和 Rust 的环境中，从仓库根目录构建：

```bash
conda activate geneminer2uce
make build
```

入口程序是 `cli/geneminer2`。拉取包含源码改动的版本后，再跑一次 `make build`。完整依赖见[中文命令行指南](manual/ZH_CN/command_line.md)。

## 快速开始

样本表为 tab 分隔：样本名、R1 和可选 R2。参考目录中每个 locus 单独一个 FASTA，文件主名就是 locus 名。

```text
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz
Sample_B\t/path/to/Sample_B_R1.fq.gz
```

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

先跑 UCE 恢复：

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

先检查 `output/uce_assembly_summary.csv`、`output/uce_rescue_summary.csv` 和 `output/uce_contigs/`。完整输出见[中文输出文件说明](manual/ZH_CN/output.md)。

## Original 模式

**用途。** 用于 exon、SCO 和核或线粒体 marker 的参考引导恢复；不写 `--assembly-mode` 时默认使用。默认后端为单线程确定性的 `original-rust`，结果按参考坐标裁切。

**运行。**

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

**说明。** `original` 是工作流模式。仅在与固定上游 Python 实现严格对照时，才加 `--assembler-implementation original`。见 [Assembler 章节](docs/assembler_ZH.md)。

## UCE 模式

**用途。** 用于 genome skimming 或 target capture 的 UCE 恢复。它弱化短 probe 边界、默认跳过参考引导 `trim`，并保留有 reads 支持的侧翼；任一 mate 通过过滤时，整对 reads 都保留。

**运行。**

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

**说明。** UCE 只使用 `uce-rust`。rescue 用第一轮 contig 加原始参考再招一次 reads；质量变差时退回第一轮。见 [Assembler 章节](docs/assembler_ZH.md)。

## Mito 模式

**用途。** 仅适用于**常规单条环形动物线粒体基因组**，并要求带注释 GenBank 参考；不适合多分子、严重重排、强异质性、复杂植物或真菌线粒体。

**运行。**

```bash
cli/geneminer2 mito -f samples.tsv -o mito_output -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

**说明。** 只有 overlap、局部 read graph 和 junction-spanning reads 都支持时才报告闭环；不会用参考坐标补洞。见 [线粒体章节](docs/5.mito.md)。

## Profiling 模式

**用途。** 不组装，只给 reads 证据：先招募 marker 相关 reads，再用 Themisto 对参考库做 k-mer 伪比对。适用于 ITS2、16S 或其他有参考库的 marker。

**运行。**

```bash
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o output -p 8
```

**说明。** 主结果是 `<output>/<sample>/marker_profile/marker_reference_support.tsv`，按参考序列报告命中、分数化共享支持和单例支持。一个 read 有 N 个候选时，每个候选只计 `1/N`；这是参考相容证据，不是唯一鉴定或生物量丰度。可选 `--profile-group-map` 只补充 `group` 注释。见 [Profiling 章节](docs/profiling_ZH.md)。

## Population 模式

**用途。** 用已完成 UCE 组装的多个样本及原始 reads，生成公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。

**运行。**

```bash
cli/geneminer2 population -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --engine panrefv2
```

**说明。** `panrefv2` 冻结 core、受限救援侧翼，默认只写入通过 QC 的 locus。需要 minibwa、samtools、bcftools 与 PLINK 1.9；ADMIXTURE 可选。见 [Population 章节](docs/population_ZH.md)。

## 文档

| 主题 | 文档 |
| --- | --- |
| 命令与依赖 | [中文命令行指南](manual/ZH_CN/command_line.md) |
| 输出文件 | [中文输出文件说明](manual/ZH_CN/output.md) |
| Reads 过滤 | [Filter 章节](docs/filter_ZH.md) |
| 组装器 | [Assembler 章节](docs/assembler_ZH.md) |
| 线粒体 | [线粒体章节](docs/5.mito.md) |
| Marker profiling | [Profiling 章节](docs/profiling_ZH.md) |
| UCE 群体遗传 | [Population 章节](docs/population_ZH.md) |
| 更新记录 | [CHANGELOG](CHANGELOG.md) |

## 引用与联系

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

相关工具包括 [Easy353](https://doi.org/10.1093/molbev/msac261) 和 [GeneMiner](https://doi.org/10.1111/1755-0998.13924)。使用 `--alignment-filter alifilter` 时，请同时引用 [AliFilter](https://doi.org/10.1093/molbev/msag097)。联系维护者：[xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp)。
