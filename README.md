# GeneMiner2 UCE CLI

**[查看中文版说明](README_zh_cn.md)**

This repository is a command-line focused fork of GeneMiner2 for target-enrichment and UCE workflows. GUI project files, graphical interface documentation, screenshots, and bundled demo datasets have been removed so the repository only contains the CLI source, build files, and command-line documentation.

## Main Features

- Reference-guided marker recovery from next-generation sequencing reads.
- UCE-oriented assembly mode through `--assembly-mode uce`, designed to retain longer read-supported flanking contigs.
- Optional one-round UCE raw-read rescue through `--uce-rescue-reads`.
- phyluce-compatible UCE contig export under `uce_contigs/`.
- Optional AliFilter alignment filtering through `--alignment-filter alifilter`.
- Controlled combine-stage parallelism through `--msa-threads` and `--filter-processes`.

## Build

Install the full build dependencies, then run:

```bash
make
```

The CLI entry point is generated as:

```bash
cli/geneminer2
```

For full build instructions and runtime dependencies, see [manual/EN_US/command_line.md](manual/EN_US/command_line.md).

## Minimal Usage

Prepare a tab-delimited sample list:

```text
Sample_A	/path/to/Sample_A_R1.fq.gz	/path/to/Sample_A_R2.fq.gz
Sample_B	/path/to/Sample_B_R1.fq.gz	/path/to/Sample_B_R2.fq.gz
```

Prepare a reference directory where each target locus is one FASTA file, for example:

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

Run the default UCE workflow:

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  --assembly-mode uce \
  --uce-rescue-reads
```

## Documentation

- [Command-line usage](manual/EN_US/command_line.md)
- [Output files](manual/EN_US/output.md)
- [中文命令行说明](manual/ZH_CN/command_line.md)
- [中文输出文件说明](manual/ZH_CN/output.md)

## Citations

Primary reference for GeneMiner2:

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. https://doi.org/10.1111/1755-0998.70111

Related earlier work:

Zhang Z, Xie PL, Guo YL, Zhou WB, Liu EY, Yu Y. 2022. **Easy353**: A tool to get Angiosperms353 genes for phylogenomic research. *Molecular Biology and Evolution* 39(12): msac261. https://doi.org/10.1093/molbev/msac261

Xie PL, Guo YL, Teng Y, Zhou WB, Yu Y. 2024. **GeneMiner**: A tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 24(3): e13924. https://doi.org/10.1111/1755-0998.13924

If `--alignment-filter alifilter` is used, please also cite:

Bianchini G, Zhu R, Cicconardi F, Moody ERR. 2026. **AliFilter: a machine learning approach to alignment filtering.** *Molecular Biology and Evolution* 43(4): msag097. https://doi.org/10.1093/molbev/msag097
