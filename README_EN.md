# GeneMiner2-UCE

**[中文主页](README.md)**

GeneMiner2-UCE is the UCE extension of GeneMiner2 for target-enrichment, ultraconserved-element (UCE), and related short-read data. It provides reference-guided marker recovery, UCE flank recovery, marker profiling, ordinary animal mitochondrial assembly, and UCE population-genetic analysis.

Please cite the [GeneMiner2-UCE GitHub repository](https://github.com/GUIBA-EX/GeneMiner2-UCE) when using this software.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Features

- Recover markers from genome-skimming or target-capture reads.
- Preserve UCE cores and read-supported flanks.
- Assess support for reference sequences from WGS or metagenomic reads.
- Build UCE cohort pseudo-references and generate joint VCF, PCA, and ADMIXTURE inputs.
- Recover ordinary circular animal mitochondria from annotated GenBank references.

## Choose a mode

| Mode | Suitable data | Main result |
| --- | --- | --- |
| `--assembly-mode original` | Exons, SCOs, and nuclear or mitochondrial markers | Reference-guided contigs with default reference trimming |
| `--assembly-mode uce` | UCEs from genome skimming or target capture | UCE cores and read-supported flanks |
| `profiling` | Any marker in WGS or metagenomic data | Per-reference support evidence |
| `mito` | Ordinary circular animal mitochondria | Read-supported circular or partial assembly |
| `population` | Multiple samples with completed UCE assemblies | Cohort pseudo-reference, VCF, PCA, and ADMIXTURE inputs |

## Installation

Build in an environment containing Cython, PyInstaller, and Rust:

```bash
conda activate geneminer2uce
make build
```

The entry point is `cli/geneminer2`. Run `make build` again after pulling a version with source changes. See the [command-line guide](manual/EN_US/command_line.md) for all dependencies.

## Quick start

Use a tab-delimited sample table: sample name, R1, and optional R2. Store each locus in a separate FASTA file; its file stem is the locus name.

```text
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz
Sample_B\t/path/to/Sample_B_R1.fq.gz
```

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

Start with UCE recovery:

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

Inspect `output/uce_assembly_summary.csv`, `output/uce_rescue_summary.csv`, and `output/uce_contigs/` first. See the [output guide](manual/EN_US/output.md) for the complete layout.

## Original mode

**Purpose.** For reference-guided recovery of exons, SCOs, and nuclear or mitochondrial markers; it is the default when `--assembly-mode` is omitted. Its default backend is deterministic, single-threaded `original-rust`, and results are trimmed to reference coordinates.

**Run.**

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

**Notes.** `original` is a workflow mode. Add `--assembler-implementation original` only for strict comparison with the fixed upstream Python implementation. See the [Assembler chapter](docs/assembler_EN.md).

## UCE mode

**Purpose.** For UCE recovery from genome skimming or target capture. It reduces short-probe boundary effects, skips reference-guided `trim` by default, and preserves read-supported flanks; the full paired-end fragment is retained when either mate passes filtering.

**Run.**

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

**Notes.** UCE uses `uce-rust` only. Rescue recruits again with first-round contigs plus original references and reverts on worse quality. See the [Assembler chapter](docs/assembler_EN.md).

## Mito mode

**Purpose.** Limited to **ordinary single circular animal mitochondrial genomes** with an annotated GenBank reference. It is not intended for multipartite genomes, major rearrangements, strong heteroplasmy, or complex plant and fungal mitochondria.

**Run.**

```bash
cli/geneminer2 mito -f samples.tsv -o mito_output -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

**Notes.** Circularity requires agreement among overlap, local read-graph paths, and junction-spanning reads; reference coordinates never fill gaps. See the [Mitochondrial chapter](docs/5.mito.md).

## Profiling mode

**Purpose.** A read-evidence workflow, not an assembler: it recruits marker-related reads and pseudoaligns them with Themisto to a reference library. It supports ITS2, 16S, and other markers with reference libraries.

**Run.**

```bash
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o output -p 8
```

**Notes.** The main result, `<output>/<sample>/marker_profile/marker_reference_support.tsv`, reports hit, fractional shared, and singleton support per reference. A query with N candidates contributes `1/N` to each. It is reference-compatibility evidence, not unique identification or organismal abundance. Optional `--profile-group-map` only adds a `group` annotation. See the [Profiling chapter](docs/profiling_EN.md).

## Population mode

**Purpose.** Uses completed UCE assemblies and original reads from multiple samples to generate a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs.

**Run.**

```bash
cli/geneminer2 population -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --engine panrefv2
```

**Notes.** `panrefv2` freezes cores, rescues flanks under constraints, and writes QC-passing loci by default. It requires minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional. See the [Population chapter](docs/population_EN.md).

## Documentation

| Topic | Document |
| --- | --- |
| Commands and dependencies | [Command-line guide](manual/EN_US/command_line.md) |
| Output files | [Output-file guide](manual/EN_US/output.md) |
| Read filtering | [Filter chapter](docs/filter_EN.md) |
| Assemblers | [Assembler chapter](docs/assembler_EN.md) |
| Mitochondria | [Mitochondrial chapter](docs/5.mito.md) |
| Marker profiling | [Profiling chapter](docs/profiling_EN.md) |
| UCE population genetics | [Population chapter](docs/population_EN.md) |
| Release history | [CHANGELOG](CHANGELOG.md) |

## Citation and contact

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Related tools include [Easy353](https://doi.org/10.1093/molbev/msac261) and [GeneMiner](https://doi.org/10.1111/1755-0998.13924). When using `--alignment-filter alifilter`, also cite [AliFilter](https://doi.org/10.1093/molbev/msag097). Contact the maintainer at [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).
