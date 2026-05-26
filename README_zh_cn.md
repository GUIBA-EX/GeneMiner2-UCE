# GeneMiner2 UCE CLI

**[English README](README.md)**

本仓库是面向 target-enrichment 和 UCE 数据流程的 GeneMiner2 命令行分支。仓库已经移除 GUI 工程、图形界面说明、截图和内置演示数据，只保留 CLI 源码、构建文件和命令行文档。

## 主要功能

- 基于参考序列从二代测序 reads 中恢复目标分子标记。
- 通过 `--assembly-mode uce` 启用 UCE 组装模式，优先保留更长且有 reads 支持的侧翼 contig。
- 通过 `--uce-rescue-reads` 启用一轮 UCE raw-read rescue。
- 在 `uce_contigs/` 下导出 phyluce 兼容的 UCE contig 文件。
- 通过 `--alignment-filter alifilter` 支持可选 AliFilter 比对列过滤。
- 通过 `--msa-threads` 和 `--filter-processes` 控制 combine 阶段并行。

## 构建

安装完整构建依赖后运行：

```bash
make
```

CLI 入口会生成在：

```bash
cli/geneminer2
```

完整构建方法和运行时依赖见 [manual/ZH_CN/command_line.md](manual/ZH_CN/command_line.md)。

## 最小用法

准备 tab 分隔的样本列表：

```text
Sample_A	/path/to/Sample_A_R1.fq.gz	/path/to/Sample_A_R2.fq.gz
Sample_B	/path/to/Sample_B_R1.fq.gz	/path/to/Sample_B_R2.fq.gz
```

准备参考序列目录，每个目标 locus 一个 FASTA 文件，例如：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

运行默认 UCE 流程：

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  --assembly-mode uce \
  --uce-rescue-reads
```

## 文档

- [命令行说明](manual/ZH_CN/command_line.md)
- [输出文件说明](manual/ZH_CN/output.md)
- [English command-line usage](manual/EN_US/command_line.md)
- [English output files](manual/EN_US/output.md)

## 引用

GeneMiner2 主要引用：

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. https://doi.org/10.1111/1755-0998.70111

相关前期工具：

Zhang Z, Xie PL, Guo YL, Zhou WB, Liu EY, Yu Y. 2022. **Easy353**: A tool to get Angiosperms353 genes for phylogenomic research. *Molecular Biology and Evolution* 39(12): msac261. https://doi.org/10.1093/molbev/msac261

Xie PL, Guo YL, Teng Y, Zhou WB, Yu Y. 2024. **GeneMiner**: A tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 24(3): e13924. https://doi.org/10.1111/1755-0998.13924

如果使用 `--alignment-filter alifilter`，也请引用：

Bianchini G, Zhu R, Cicconardi F, Moody ERR. 2026. **AliFilter: a machine learning approach to alignment filtering.** *Molecular Biology and Evolution* 43(4): msag097. https://doi.org/10.1093/molbev/msag097
