# GeneMiner2-UCE

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[中文](README.md)** · [Changelog](CHANGELOG.md) · [Report an issue](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

GeneMiner2-UCE is a reference-guided short-read recovery toolkit for genome skimming, target capture, and related data. Rust is the default compute backend. The project supports UCEs, conventional markers, mitochondria, marker profiling, UCE population genetics, nuclear gene families, and repeatomes.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Choose a workflow

| Goal | Command or mode | Main result |
| --- | --- | --- |
| Recover exons, SCOs, or conventional nuclear/mitochondrial markers | `--assembly-mode original` (default) | Reference-guided contigs, trimmed to the reference by default |
| Recover UCE cores and read-supported flanks | `--assembly-mode uce` | UCE contigs plus assembly and rescue summaries |
| Assess read support without assembly | `profiling` | Per-reference support evidence |
| Assemble ordinary circular animal mitochondria | `mito` | Read-supported closed or partial result |
| Analyse variation across UCE samples | `population` | Cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs |
| Recover nuclear gene-family candidates | `gene` | Family status, candidate FASTA, and cohort summaries |
| Build a repeatome from short reads | `te` | Conservative repeat library, annotation, and sample RPM |

`original` uses deterministic `original-rust` by default. UCE always uses Rust `ucefilter` and `uce-rust`; it does not fall back to Python.

## Install

Rust/Cargo, Python 3.11, Cython, PyInstaller, and bioinformatics dependencies are required. Follow the [command-line guide](manual/EN_US/command_line.md) to create the environment, then run this from the repository root:

```bash
make build
```

The entry point is `cli/geneminer2`. Run `make build` again after updating source code.

## Quick start: UCE

Use a tab-delimited sample manifest with `sample_id  R1  [R2]`. UCE requires paired reads. Put one FASTA per locus in the reference directory; the filename stem is the locus name.

```text
samples.tsv
Sample_A\t/path/to/Sample_A_R1.fq.gz\t/path/to/Sample_A_R2.fq.gz

references/
  uce-0001.fasta
  uce-0002.fasta
```

```bash
cli/geneminer2 filter assemble \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce --uce-rescue-reads
```

Inspect these outputs first:

- `output/uce_assembly_summary.csv`: locus-level assembly choice and QC;
- `output/uce_rescue_summary.csv`: optional rescue gains and rollbacks;
- `output/uce_contigs/`: final UCE contigs.

Without `--uce-rescue-reads`, only the primary UCE recovery runs. The option enables at most two constrained extension rounds.

## Other common entry points

```bash
# Conventional markers; original is also the default when --assembly-mode is omitted
cli/geneminer2 filter assemble -f samples.tsv -r references -o original_out -p 8 \
  --assembly-mode original

# Nuclear gene families
cli/geneminer2 gene -f samples.tsv -r family_references -o gene_out -p 8

# Ordinary circular animal mitochondria; requires an annotated GenBank reference
cli/geneminer2 mito -f samples.tsv -o mito_out -p 8 \
  --mito-genbank mitochondrial_reference.gb

# Marker read evidence, without assembly
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o profile_out -p 8

# Population workflow after UCE recovery
cli/geneminer2 population -f samples.tsv -r references -o population_out -p 8 \
  --assembly-mode uce --engine panrefv2

# Reference-free repeatome; uses a separate TE sample manifest
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

`mito` is limited to ordinary single circular animal mitochondrial genomes; it is not intended for multipartite, strongly rearranged, or complex plant/fungal mitochondria. `profiling` measures reference compatibility, not unique taxonomic identification or abundance.

## Documentation

| Task | Document |
| --- | --- |
| Installation, all options, and external dependencies | [Command-line guide](manual/EN_US/command_line.md) |
| Output directories and tables | [Output-file guide](manual/EN_US/output.md) |
| Read filtering and cache | [Filter](docs/filter_EN.md) |
| Original/UCE assembly and rescue | [Assembler](docs/assembler_EN.md) |
| Marker read profiling | [Profiling](docs/profiling_EN.md) |
| UCE population genetics | [Population](docs/population_EN.md) |
| Nuclear gene families | [Gene](docs/gene_EN.md) |
| Repeatome | [TE / repeatome](docs/te_EN.md) |
| Mitochondria | [Mito](docs/mitochondria_EN.md) |
| Performance design and scope | [MainFilter performance note](docs/development/mainfilter-performance.md) |

Add `--workflow-profile` to write `workflow_profile.tsv` at the output root. Supported Rust assemblers also write per-sample `assembly_profile.tsv`. Both record timing and I/O only; neither changes results.

## Citation and support

Please cite: Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Report bugs and feature requests through [GitHub Issues](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues). This project is released under [GPL-3.0-or-later](LICENSE); see [NOTICE](NOTICE) for provenance boundaries of third-party and ported code.
