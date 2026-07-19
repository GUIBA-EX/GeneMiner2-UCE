
**[中文主页](README.md)**

GeneMiner2-UCE is the UCE extension of GeneMiner2 for target-enrichment, ultraconserved-element (UCE), and related short-read data. It retains reference-guided read recruitment while extending UCE flank recovery, marker profiling, and population-genetic analysis of UCE data.

Please cite the [GeneMiner2-UCE GitHub repository](https://github.com/GUIBA-EX/GeneMiner2-UCE) when using this software. A formal publication will be added when available.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Main features

- Recover molecular markers from genome-skimming or target-capture reads.
- Preserve UCE cores and read-supported flanking sequences.
- Recruit marker-associated reads from WGS or other shotgun data and report reference-level support by k-mer pseudoalignment.
- Build a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs from UCE samples.
- Export PHYLUCE-compatible contigs and summarize recovery quality by sample and locus.

## Choose a mode

| Mode | Suitable data | Main result |
| --- | --- | --- |
| `--assembly-mode original` | Exons, SCOs, and nuclear or mitochondrial markers | Reference-guided contigs with reference trimming in the default workflow |
| `--assembly-mode uce` | UCE recovery from genome skimming or target capture | UCE cores and read-supported flanking sequences |
| `profiling` subcommand | Any amplicon marker in WGS or metagenomic data | Per-reference hit, fractional, and singleton support |
| `mito` subcommand | Ordinary circular animal mitochondria with an annotated GenBank reference | Read-supported circular or partial mitochondrial assembly |
| `population` subcommand | Multiple samples with completed UCE assemblies | Cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs |

## Installation

GeneMiner2-UCE is built from source. Activate an environment containing Cython, PyInstaller, and Rust, then run in the repository root:

```bash
conda activate geneminer2uce
make build
```

The resulting entry point is:

```bash
cli/geneminer2
```

Run `make build` again after pulling a release that changes the source. See the [command-line guide](manual/EN_US/command_line.md) for complete build and runtime dependencies.

## Quick start

Prepare a tab-delimited sample table. Each row contains a sample name, R1, and an optional R2:

```text
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz
Sample_B\t/path/to/Sample_B_R1.fq.gz
```

Store each locus in a separate FASTA file. The file stem becomes the locus name:

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

Run UCE recovery with one controlled raw-read rescue round:

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --uce-rescue-reads
```

Inspect these results first:

- `output/uce_assembly_summary.csv`: assembly status and support metrics by sample and locus;
- `output/uce_rescue_summary.csv`: before/after rescue comparison and fallback status;
- `output/uce_contigs/`: sample contigs formatted for PHYLUCE.

See the [output guide](manual/EN_US/output.md) for the complete directory layout.

## Original mode

`original` is the default assembly mode for reference-guided recovery of exons, SCOs, and nuclear or mitochondrial markers; it is selected when `--assembly-mode` is omitted. Its default backend is `original-rust`. Use `--assembler-implementation original` for a strict comparison with the fixed upstream GeneMiner2 Python original. `uce-rust` is the UCE-oriented backend; use `--assembly-mode uce` for routine UCE recovery.

```bash
cli/geneminer2 \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original` names the workflow mode, while `--assembler-implementation original` selects the optional Python backend within that mode.

## UCE mode

UCE mode reduces the influence of short-probe boundaries, skips reference-guided `trim` in the default workflow, and favors longer candidates that retain read support. During refiltering, the complete paired-end fragment is retained whenever either mate passes the locus filter.

The default Rust assembler follows a backbone strategy without repeated backtracking. `--uce-rescue-reads` recruits reads again using first-round contigs plus the original references, and restores the first-round result when rescue quality deteriorates. See the [Assembler chapter](docs/assembler_EN.md) for parameters, guardrails, reference caching, and fallback rules.

## Mito mode

`mito` is limited to **ordinary single circular animal mitochondrial genomes**. It recruits reads from an annotated GenBank reference and reports circularity only when overlap, local read-graph paths, and junction-spanning reads support it. It does not fill gaps from reference coordinates and is not intended for multipartite genomes, major rearrangements, strong heteroplasmy, or complex plant and fungal mitochondria.

```bash
cli/geneminer2 mito \
  -f samples.tsv \
  -o mito_output \
  -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

See the [Mitochondrial chapter](docs/5.mito.md) for the workflow, success criteria, and expert parameters.

## Profiling mode

`profiling` is a **read-level evidence workflow, not an assembler**. GeneMiner2 recruits marker-related reads once, then Themisto pseudoaligns them to the reference library. The primary result is support for each reference sequence; it does not run `refilter`, `assemble`, `combine`, or `tree`.

```bash
cli/geneminer2 profiling \
  -f samples.tsv \
  -r marker_reference.fasta \
  -o output -p 8
```

`<output>/<sample>/marker_profile/marker_reference_support.tsv` reports hit queries, fractional shared-query support, and singleton support for every hit reference. A query with N candidates contributes `1/N` to each candidate, so shared reads are not counted N times. This is compatibility evidence, not a unique identification or organismal abundance.

`--profile-group-map` is optional and only fills the `group` annotation column with a species, clade, or curated DIV label. It does not calculate group abundance. See the [Profiling chapter](docs/profiling_EN.md) for cache and QC details.

## Population mode

`population` uses multiple samples with completed UCE assemblies and their original reads to build a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs. It creates or reuses the pseudo-reference, maps every sample to it, performs joint variant calling, and emits all-SNP, one-SNP-per-UCE, and LD-pruned panels.

After UCE assembly is complete, run:

```bash
cli/geneminer2 population \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --engine panref \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

Runtime dependencies are minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional. See the [Population chapter](docs/population_EN.md) for pseudo-reference strategies, staged restarts, SNP panels, and required QC checks.

## Documentation

The README is an entry point; parameter definitions, QC rules, and output fields live in the chapters below. `original` uses `original-rust` by default; select the fixed [upstream GeneMiner2 Python assembler](https://github.com/sculab/GeneMiner2/blob/36e06feeb99654bdb87f45d4cde225d8c3e311d0/scripts/main_assembler.py) with `--assembler-implementation original` only for strict comparison. `uce` uses `uce-rust` only.

- [Command-line guide](manual/EN_US/command_line.md)
- [Output-file guide](manual/EN_US/output.md)
- [Filter chapter](docs/filter_EN.md)
- [Assembler chapter](docs/assembler_EN.md)
- [Mitochondrial chapter](docs/5.mito.md)
- [Profiling chapter](docs/profiling_EN.md)
- [Population chapter](docs/population_EN.md)
- [Release history](CHANGELOG.md)

## Citation and contact

Primary GeneMiner2 reference:

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Related tools include [Easy353](https://doi.org/10.1093/molbev/msac261) and [GeneMiner](https://doi.org/10.1111/1755-0998.13924). When using `--alignment-filter alifilter`, also cite [AliFilter](https://doi.org/10.1093/molbev/msag097).

Contact the maintainer at [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).
