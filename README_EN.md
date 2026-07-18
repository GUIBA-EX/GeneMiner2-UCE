# GeneMiner2-UCE

**[中文主页](README.md)**

GeneMiner2-UCE is the UCE extension of GeneMiner2 for target-enrichment, ultraconserved-element (UCE), and related short-read data. It retains reference-guided read recruitment while extending UCE flank recovery, multi-candidate ITS2 assembly, and population-genetic analysis of UCE data.

Please cite the [GeneMiner2-UCE GitHub repository](https://github.com/GUIBA-EX/GeneMiner2-UCE) when using this software. A formal publication will be added when available.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Main features

- Recover molecular markers from genome-skimming or target-capture reads.
- Preserve UCE cores and read-supported flanking sequences.
- Retain multiple ITS2 candidates with paired, diagnostic, and EM-abundance evidence.
- Build a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs from UCE samples.
- Export PHYLUCE-compatible contigs and summarize recovery quality by sample and locus.

## Choose a mode

| Mode | Suitable data | Main result |
| --- | --- | --- |
| `--assembly-mode reference` | Genome skimming and conventional gene recovery | Reference-guided contigs with reference trimming in the default workflow |
| `--assembly-mode uce` | UCE target capture | UCE cores and read-supported flanking sequences |
| `--assembly-mode its2` | Loci with several genuine ITS2 variants | Multiple candidates with fragment-level support |
| `population` subcommand | Diploid UCE resequencing or target capture | Cohort pseudo-reference, joint SNPs, PCA, and ADMIXTURE panels |

## Installation

GeneMiner2-UCE is currently built from source. Run in the repository root:

```bash
make
```

The resulting entry point is:

```bash
cli/geneminer2
```

Run `make` again after pulling a release that changes the source. See the [command-line guide](manual/EN_US/command_line.md) for complete build and runtime dependencies.

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

## UCE mode

UCE mode reduces the influence of short-probe boundaries, skips reference-guided `trim` in the default workflow, and favors longer candidates that retain read support. During refiltering, the complete paired-end fragment is retained whenever either mate passes the locus filter.

The default Rust assembler follows a backbone strategy without repeated backtracking. `--uce-rescue-reads` recruits reads again using first-round contigs plus the original references, and restores the first-round result when rescue quality deteriorates. See the [UCE workflow guide](docs/uce-workflow_EN.md) for parameters, guardrails, reference caching, and fallback rules.

## ITS2 mode

ITS2 mode requires the Rust assembler and uses a fixed 21-mer. It can emit several candidates per locus and reports fragment, paired-fragment, diagnostic-fragment, and EM-abundance evidence. Candidates that reads cannot distinguish remain in an equivalence group.

The standard `combine` stage expects one sequence per locus, so run ITS2 explicitly as:

```bash
cli/geneminer2 filter refilter assemble \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode its2
```

ITS2 cannot use the original Python assembler and does not fall back when the Rust assembler fails. Results are written under each sample's `results/`, `its2_assembly_summary.csv`, and `<locus>.its2_support.tsv`.

## Population mode

`population` builds a consistent unphased diploid SNP matrix from accepted UCE contigs and the original reads. It creates or reuses a cohort pseudo-reference, maps every sample to that reference, performs joint variant calling, and emits all-SNP, one-SNP-per-UCE, and LD-pruned panels.

After UCE assembly is complete, run:

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

Runtime dependencies are minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional. See the [Population workflow guide](docs/population_EN.md) for pseudo-reference strategies, staged restarts, SNP panels, and required QC checks.

## Implementation and documentation

The default build includes the Rust MainFilter, Refilter, Assembler, Population workflow, the optional `main_assembler-original-rust` compatibility binary, and Rust helper tools. In `reference` mode, `--assembler-implementation auto` uses `original-rust`, the deterministic single-thread Rust port of the original reference algorithm. Choose `original` for the byte-identical [upstream GeneMiner2 assembler](https://github.com/sculab/GeneMiner2/blob/36e06feeb99654bdb87f45d4cde225d8c3e311d0/scripts/main_assembler.py), or `uce-rust` to test the UCE-oriented Rust assembler. With `--reuse-reference-cache`, `original-rust` reuses a binary k-mer cache validated by format version, implementation version, k, and reference-file identity; stale or corrupt files are rebuilt automatically. UCE and ITS2 are Rust-only: an unavailable or failed Rust assembler is reported as an error with no Python fallback. `original` and `original-rust` are reference-only. The main CLI orchestrator and consensus program remain in Python.

- [Command-line guide](manual/EN_US/command_line.md)
- [Output-file guide](manual/EN_US/output.md)
- [UCE workflow guide](docs/uce-workflow_EN.md)
- [Population workflow guide](docs/population_EN.md)
- [MainFilter performance and compatibility note](docs/mainfilter-performance.md)
- [Release history](CHANGELOG.md)

## Citation and contact

Primary GeneMiner2 reference:

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Related tools include [Easy353](https://doi.org/10.1093/molbev/msac261) and [GeneMiner](https://doi.org/10.1111/1755-0998.13924). When using `--alignment-filter alifilter`, also cite [AliFilter](https://doi.org/10.1093/molbev/msag097).

Contact the maintainer at [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).
