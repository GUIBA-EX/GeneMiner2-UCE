# GeneMiner2-UCE

**[中文主页](README.md)**

GeneMiner2-UCE is the UCE extension of GeneMiner2 for target-enrichment, ultraconserved-element (UCE), and related short-read data. It retains reference-guided read recruitment while extending UCE flank recovery, marker profiling, and population-genetic analysis of UCE data.

Please cite the [GeneMiner2-UCE GitHub repository](https://github.com/GUIBA-EX/GeneMiner2-UCE) when using this software. A formal publication will be added when available.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Main features

- Recover molecular markers from genome-skimming or target-capture reads.
- Preserve UCE cores and read-supported flanking sequences.
- Recruit marker-associated reads from WGS or other shotgun data and quantify them directly by k-mer pseudoalignment.
- Build a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs from UCE samples.
- Export PHYLUCE-compatible contigs and summarize recovery quality by sample and locus.

## Choose a mode

| Mode | Suitable data | Main result |
| --- | --- | --- |
| `--assembly-mode original` | Exons, SCOs, and nuclear or mitochondrial markers | Reference-guided contigs with reference trimming in the default workflow |
| `--assembly-mode uce` | UCE recovery from genome skimming or target capture | UCE cores and read-supported flanking sequences |
| `profiling` subcommand | Any amplicon marker in WGS or metagenomic data | Relative marker-group signal, detection state, and QC |
| `population` subcommand | Multiple samples with completed UCE assemblies | Cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs |

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

## Profiling mode

`profiling` is a **read-level quantification workflow, not an assembler**. It supports any amplicon marker: GeneMiner2 performs one k-mer recruitment, Themisto pseudoaligns the recruited reads to the marker library, and mSWEEP estimates relative signal among reference groups that share reads. It does not run `refilter`, `assemble`, `combine`, or `tree`.

Pass one `.fasta` or `.fa` marker library directly with `-r`, plus a required two-column `--profile-group-map` TSV: `reference_id<TAB>group`. The reference ID is the first whitespace-delimited FASTA-header field, and every reference must map to exactly one group. `themisto` and `mSWEEP` must be available, either on `PATH` or supplied via `--profile-themisto` and `--profile-msweep`.

```bash
cli/geneminer2 profiling \
  -f samples.tsv \
  -r marker_reference.fasta \
  --profile-group-map marker_groups.tsv \
  -o output \
  -p 8 \
  --profile-decoy non_target_sequences.fasta
```

The primary per-sample result is `<output>/<sample>/marker_profile/marker_group_abundance.tsv`; `marker_qc.tsv` records pseudoalignment and mSWEEP statistics, while `marker_reference_metadata.tsv` records the reference-color/group mapping. `evidence_queries` and `exclusive_queries` count individual FASTA/FASTQ query records, not paired fragments. `relative_proportion` is an uncalibrated marker-signal proportion renormalized after the minimum-exclusive-evidence rule; it is not a cell or organism proportion.

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
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

Runtime dependencies are minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional. See the [Population chapter](docs/population_EN.md) for pseudo-reference strategies, staged restarts, SNP panels, and required QC checks.

## Implementation and documentation

The default build includes Rust MainFilter, Refilter, Assembler, Population, marker-profiling helpers, and other Rust utilities. `original` uses `original-rust` by default and can select the fixed [upstream GeneMiner2 Python assembler](https://github.com/sculab/GeneMiner2/blob/36e06feeb99654bdb87f45d4cde225d8c3e311d0/scripts/main_assembler.py) with `--assembler-implementation original` for strict comparison; `uce` uses only `uce-rust`. With `--reuse-reference-cache`, `original-rust` reuses a binary cache validated by format, implementation version, k, and reference-file identity; stale or corrupt files are rebuilt automatically. The main CLI orchestrator and consensus program remain in Python.

- [Command-line guide](manual/EN_US/command_line.md)
- [Output-file guide](manual/EN_US/output.md)
- [Filter chapter](docs/filter_EN.md)
- [Assembler chapter](docs/assembler_EN.md)
- [Profiling chapter](docs/profiling_EN.md)
- [Population chapter](docs/population_EN.md)
- [Release history](CHANGELOG.md)

## Citation and contact

Primary GeneMiner2 reference:

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Related tools include [Easy353](https://doi.org/10.1093/molbev/msac261) and [GeneMiner](https://doi.org/10.1111/1755-0998.13924). When using `--alignment-filter alifilter`, also cite [AliFilter](https://doi.org/10.1093/molbev/msag097).

Contact the maintainer at [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).
