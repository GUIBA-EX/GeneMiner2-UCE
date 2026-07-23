# TStools (formerly GeneMiner2-UCE)

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Fuzz smoke](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/fuzz-smoke.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/fuzz-smoke.yml)
[![MSRV: 1.87](https://img.shields.io/badge/MSRV-1.87-orange)](rust-toolchain.toml)
[![Rust edition: 2021](https://img.shields.io/badge/Rust%20edition-2021-orange)](Cargo.toml)
[![Dependency policy: cargo-deny](https://img.shields.io/badge/dependency%20policy-cargo--deny-blue)](deny.toml)
[![SBOM: SPDX 2.3](https://img.shields.io/badge/SBOM-SPDX%202.3-blueviolet)](rust/xtask/src/main.rs)
[![Release integrity: SHA-256](https://img.shields.io/badge/release%20integrity-SHA--256-blueviolet)](rust/xtask/src/main.rs)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[中文](README.md)** · [Changelog](CHANGELOG.md) · [Report an issue](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

TStools (formerly GeneMiner2-UCE) is a reference-guided short-read recovery toolkit: it recruits reads with references, then assembles, quantifies evidence, or performs cohort analysis by workflow. It covers genome skimming, target capture, UCEs, mitochondria, nuclear gene families, RAD augmentation, and reference-free repeatomes. Production workflows are Rust-native and require no Python runtime.

> **Relationship to GeneMiner2:** GeneMiner2 is this project's algorithmic provenance and result-compatibility baseline. TStools is not a maintenance fork: it is an independently evolving workflow toolkit. `cli/geneminer2`, some legacy output names, the repository URL, and historical citation names are compatibility identifiers; they do not mean that the two projects have identical features or implementations.

![TStools workflow](docs/images/summary_EN.png)

## Main differences from GeneMiner2

| Dimension | GeneMiner2 | Current TStools (v1.5.8) |
| --- | --- | --- |
| Positioning | Original algorithm and workflow for gene recovery from genome-skimming data | Independently evolving short-read recovery and analysis toolkit built on a compatibility baseline |
| Production implementation | Upstream implementation | Rust-native production paths; no Python runtime required |
| Read recruitment | Original recruitment semantics | Canonical bidirectional 2-bit k-mers, content-validated reference caches, and bounded streaming I/O preserve recruitment semantics and legacy output compatibility while reducing CPU, memory, and I/O cost |
| Conventional assembly | Upstream algorithmic baseline | Deterministic `original-rust` by default; upstream `original` remains available for strict comparison and reproduction |
| UCE assembly | General, non-specialized assembly route | `ucefilter → uce-rust` combines recruitment, paired-fragment retention, orientation/exact-match evidence, and per-locus selection in one FASTQ scan. Optional rescue accepts only read-supported extension and never reference-fills a gap |
| Workflow scope | Conventional gene recovery | Also includes mitochondria, marker profiling, UCE population analysis, nuclear gene families, RAD matrix augmentation, and reference-free repeatomes |
| Interpretation | Primarily recovered sequences | Mitochondrial closure, RAD strict matrices, and population graph paths require explicit evidence and retain QC, provenance, and audit output |

For an upstream-baseline comparison, select `original` explicitly and retain the complete inputs and parameters. For UCE, population, and extended workflows, interpret results using TStools-specific modes and QC rules rather than treating them as interchangeable with GeneMiner2 output.

See [Filter](docs/filter_EN.md), [Assembler](docs/assembler_EN.md), and [Population](docs/population_EN.md) for algorithms, measurements, and scope.

## Install and first run

Install Rust/Cargo and the required bioinformatics tools, then build from the repository root:

```bash
cargo run -p xtask -- build
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
- `--cleanup-intermediates` removes reproducible filtered reads only after a successful complete invocation; final contigs, summaries, raw reads, and references remain. Use `--cleanup-dry-run` first to write `cleanup_preview.tsv` for review without deleting files.

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

Each standard workflow atomically writes `workflow_manifest.tsv` in its output root, recording the CLI version, command, reference and sample-sheet SHA-256 values, key options, and input-read paths/sizes/mtimes for reproduction and audit.

`--resume` is conservative whole-workflow recovery: it returns a no-op success only when the existing manifest exactly matches the current inputs and options and `workflow_status.tsv` records success. It refuses failed or mismatched runs and never overwrites prior status or skips partial stages.

Once an output directory exists, the CLI atomically writes `workflow_status.tsv` on completion. Its `state` is `succeeded` or `failed`; failures include an error summary so batch systems can identify incomplete outputs.

`cargo run -p xtask -- build` also creates `SHA256SUMS` and `SBOM.spdx.json` in `cli/` for release integrity and the software bill of materials. FASTX fuzz targets live in `fuzz/` and run only in manually triggered or weekly, constrained CI smoke jobs (1,000 inputs, up to 60 seconds).

## Citation

Please cite: Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

```bibtex
@software{TStools,
  author  = {XIA, Fei and TANG, Zizhen and XU, Yan},
  title   = {TStools (formerly GeneMiner2-UCE): Reference-Guided Short-Read Recovery for UCE, Mitochondrial, Gene-Family, and RAD Workflows},
  year    = {2026},
  version = {1.5.8},
  url     = {https://github.com/GUIBA-EX/GeneMiner2-UCE},
  publisher = {GitHub},
  note    = {GPL-3.0-or-later licensed software}
}
```

Released under [GPL-3.0-or-later](LICENSE). See [NOTICE](NOTICE) for provenance boundaries of third-party and ported code.
