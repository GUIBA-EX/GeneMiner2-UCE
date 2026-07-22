# GeneMiner2-UCE

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[中文](README.md)** · [Changelog](CHANGELOG.md) · [Report an issue](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

GeneMiner2-UCE is a native Rust toolkit for short-read recovery: recruit reads with references, then assemble or quantify them for the selected task. It supports genome skimming, target capture, UCEs, RAD augmentation, and repeatome analysis without a Python runtime.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Install and first run

Install Rust/Cargo and the required bioinformatics tools, then build from the repository root:

```bash
make build
```

The common entry point is `cli/geneminer2`. Sample manifests are tab-delimited `sample_id  R1  [R2]`; each FASTA in a reference directory is one locus or bait.

```bash
# UCE: paired reads → selective recruitment → UCE assembly
cli/geneminer2 filter assemble \
  -f samples.tsv -r uce_references -o uce_out -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

Start with `uce_out/uce_assembly_summary.csv` and `uce_out/uce_contigs/`. `--uce-rescue-reads` is an optional, evidence-constrained extension of at most two rounds; it never reference-fills a gap.

## Choose a command

| Goal | Command | Main result |
| --- | --- | --- |
| Conventional exons, SCOs, or nuclear markers | `filter assemble` | Reference-guided contigs |
| UCE cores and read-supported flanks | `filter assemble --assembly-mode uce` | UCE contigs and QC summary |
| Mitochondria | `mito` | Closed, linear, or explicitly ambiguous structure call |
| Read support for markers | `profiling` | Per-reference support evidence |
| UCE population genetics | `population` | Cohort pseudo-reference, VCF, PCA, and more |
| Nuclear gene families | `gene` | Family candidates, copy states, and downstream inputs |
| Add WGS samples to a RAD matrix | `rad-probe` → `rad` → `rad-validate` | Independent-arm recovery and strict matrix |
| Reference-free repeatome | `te` | Conservative repeat library, annotation, and RPM |

Conventional assembly defaults to deterministic `original-rust`; UCE uses dedicated `ucefilter` and `uce-rust`. Both routes are Rust implementations.

## Common commands

```bash
# Conventional markers; original is the default mode
cli/geneminer2 filter assemble -f samples.tsv -r references -o marker_out -p 8

# Ordinary circular animal mitochondrion; requires annotated GenBank reference
cli/geneminer2 mito -f samples.tsv -o mito_out -p 8 \
  --mito-genbank mitochondrial_reference.gb

# Nuclear gene families
cli/geneminer2 gene -f samples.tsv -r family_baits -o gene_out -p 8

# Read support only, without assembly
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o profile_out -p 8

# Population workflow after UCE recovery
cli/geneminer2 population -f samples.tsv -r uce_references -o population_out -p 8 \
  --assembly-mode uce --engine panrefv2

# Reference-free repeatome; use a separate manifest
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

## Interpretation boundaries

- Default UCE fuses broad recruitment and locus-level read selection into one FASTQ scan. `--legacy-uce-filter` is only for comparison and diagnosis.
- `mito` targets ordinary single circular animal mitochondrial genomes. Short reads cannot reliably count identical tandem repeats beyond the insert size; the result remains linear or ambiguous.
- `profiling` reports reference compatibility, not taxonomic identification or abundance.
- RAD R1 and R2 are independent restriction-site arms. WGS recovery is not direct evidence of allele dropout; use the two-arm checks in `rad-validate`.
- `--cleanup-intermediates` removes reproducible filtered reads only after a successful complete invocation; final contigs, summaries, raw reads, and references remain.

## Documentation

| Topic | Document |
| --- | --- |
| Installation, options, external tools | [Command-line guide](manual/EN_US/command_line.md) |
| Output directories and tables | [Output guide](manual/EN_US/output.md) |
| Filtering and cache | [Filter](docs/filter_EN.md) |
| Conventional/UCE assembly | [Assembler](docs/assembler_EN.md) |
| Mito, Gene, RAD, and TE | [Mito](docs/mitochondria_EN.md) · [Gene](docs/gene_EN.md) · [RAD](docs/rad_EN.md) · [TE](docs/te_EN.md) |
| Population and profiling | [Population](docs/population_EN.md) · [Profiling](docs/profiling_EN.md) |

`--workflow-profile` writes `workflow_profile.tsv` with timing and I/O only; it does not change results.

## Citation

Please cite: Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

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

Released under [GPL-3.0-or-later](LICENSE). See [NOTICE](NOTICE) for provenance boundaries of third-party and ported code.
