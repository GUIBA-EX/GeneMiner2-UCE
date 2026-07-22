# GeneMiner2-UCE

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)

**[中文主页](README.md)**

GeneMiner2-UCE is the UCE extension of GeneMiner2 for target-enrichment, ultraconserved-element (UCE), and related short-read data. It provides reference-guided marker recovery, UCE flank recovery, marker profiling, ordinary animal mitochondrial assembly, UCE population genetics, nuclear gene-family analysis, and short-read repeatome analysis.

Please cite the [GeneMiner2-UCE GitHub repository](https://github.com/GUIBA-EX/GeneMiner2-UCE) when using this software.

> **Testing stage.** The project is still under active validation. Please report bugs, unexpected results, or documentation issues through [GitHub Issues](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues) or by email to [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).

> **v1.3.** MainFilter now uses canonical bidirectional k-mers, compact locus postings, and reusable dictionaries validated by a reference-content SHA-256; filtering output is unchanged for identical inputs and parameters.

> **v1.3.1.** MainFilter further reduces temporary allocation in long-k-mer and GM2/text encoding paths. At the global memory budget, distributed locus output flushes only its largest buffer. A synthetic long-k-mer stress test was about 9% faster; see the [performance notes](docs/development/mainfilter-performance.md).

> **v1.3.2.** The consensus stage, including SAM/CIGAR handling, IUPAC sequence generation, and the optional mutation-density image, now runs in Rust; the main workflow no longer depends on matplotlib or SciPy.

> **v1.4.** MainFilter output now keeps file handles open across flushes, FASTQ/FASTA parsing is byte-level, and the build automatically links native zlib-ng through `pkg-config` when available; runtime detection remains as a portable fallback to zlib-ng or system zlib. filtering output is unchanged, byte-for-byte, for identical inputs and parameters. `mito` adds a joint rescue round: contigs retained from the first UCE-style pass become sample-specific seeds combined with the GenBank baits into a rescue reference, then recruited and reassembled together with the original paired reads; the adaptive-stop check is now an exact, cut- and strand-independent circular sequence comparison.

![GeneMiner2-UCE workflow](docs/images/summary_EN.png)

## Performance profiling

Add `--workflow-profile` to write `workflow_profile.tsv` at the output root. It records wall-clock time and input/output bytes by sample and rescue round for UCEFilter recruitment/selection, the Rust assembler, and terminal QC/rollback; it does not change filtering or assembly results. When the selected Rust assembler supports it, this option also writes each sample's `assembly_profile.tsv`, splitting read decoding, k-mer counting, graph processing, path search, and output.

```bash
geneminer2 filter assemble -f samples.txt -r references -o run --assembly-mode uce --workflow-profile
```

## Features

- Recover markers from genome-skimming or target-capture reads.
- Preserve UCE cores and read-supported flanks.
- Assess support for reference sequences from WGS or metagenomic reads.
- Build UCE cohort pseudo-references and generate joint VCF, PCA, and ADMIXTURE inputs.
- Recover ordinary circular animal mitochondria from annotated GenBank references.
- Recover nuclear gene-family candidates and resolve conservative one-to-one clades.
- Build and quantify conservative repeatome units from short reads.

## Choose a mode

| Mode | Suitable data | Main result |
| --- | --- | --- |
| `--assembly-mode original` | Exons, SCOs, and nuclear or mitochondrial markers | Reference-guided contigs with default reference trimming |
| `--assembly-mode uce` | UCEs from genome skimming or target capture | UCE cores and read-supported flanks |
| `profiling` | Any marker in WGS or metagenomic data | Per-reference support evidence |
| `mito` | Ordinary circular animal mitochondria | Read-supported circular or partial assembly |
| `population` | Multiple samples with completed UCE assemblies | Cohort pseudo-reference, VCF, PCA, and ADMIXTURE inputs |
| `te` | Genome-skimming or WGS short reads | Conservative repeatome library, sample RPM, and calls |
| `gene` | Nuclear gene families defined by multi-species baits | Candidate-contig status and cohort FASTA/matrices |

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

**Notes.** UCE now uses the fused Rust `ucefilter` by default. One pass over the original PE reads performs rolling-k-mer recruitment, run-k orientation verification, maximum-exact-match evidence, and adaptive per-locus selection, then writes `filtered/` directly. Low-depth or incompletely covered loci pass through; only saturated loci spanning the reference have redundant core evidence reduced, while PE reads crossing bait or contig edges are retained in overhang-length strata. Users do not need to declare target capture versus genome skimming. No GM2 file, candidate FASTQ, or separate Refilter process is produced; each fragment is stored once and both mates remain atomic. Candidate data beyond 256 MiB per sample are written sequentially to an internal spool that is removed automatically. `refilter` remains accepted in commands for compatibility but is unnecessary in UCE mode. Every sample remains single-threaded from recruitment through rescue. With `--uce-rescue-reads`, a bounded whole-contig round is followed by a terminal round for loci that are still growing. The old core is frozen and unsupported added sides roll back. See the [Assembler chapter](docs/assembler_EN.md).

Experimental `--uce-alignment-shadow` is off by default. It samples at most 64 selected fragments per locus and records internal affine-gap alignment identity, overlap, linked-mate, terminal, and 64-bin breadth evidence without changing read selection or assembly input. Target-capture and genome-skimming evidence still require separate interpretation; a short bait boundary is not a biological contig terminal.

## TE / repeatome mode

**Purpose.** A reference-free repeatome workflow for genome-skimming or WGS short reads. `discover`, `curate`, `annotate`, and `quantify` create conservative EQ units, annotation evidence, and sample RPM in sequence. It is not a complete-TE annotation, insertion-site calling, or TE-phylogeny workflow.

**Run.**

```bash
cli/geneminer2 te -f te_samples.tsv -o te_output -p 32
```

**Notes.** TE uses its own manifest: `taxon_id sample_id read1 read2` for paired reads, or `taxon_id sample_id read1` for single-end reads; no `-r` is required. Optional `--te-library` supplies conservative library annotation and never merges EQs. See the [TE / repeatome chapter](docs/te_EN.md) for stages, thresholds, outputs, and interpretation.

## Gene mode

**Purpose.** Recover within-sample candidate contigs for nuclear gene families defined by multi-species baits; one reference FASTA defines one family. Candidate count is assembly evidence, not an allele or biological copy-number call.

**Run.**

```bash
cli/geneminer2 gene -f samples.tsv -r family_reference -o gene_output -p 8
```

**Notes.** Main results are written to `gene_output/gene/`, including candidate status, family matrices, `pseudo_sco/`, and `multiple_candidate_families/`. Optionally add `--gene-protein-reference family_proteins` for automatic protein-guided annotation. The subsequent `gene-annotate`, `gene-resolve`, and `gene-tree` subcommands prepare strict one-to-one or multicopy species-tree inputs; they require miniprot, MAFFT/IQ-TREE, and ASTER2 `astral`, respectively. See the [Gene chapter](docs/gene_EN.md) for the complete workflow, QC, and output interpretation.

## Mito mode

**Purpose.** Limited to **ordinary single circular animal mitochondrial genomes** with an annotated GenBank reference. It is not intended for multipartite genomes, major rearrangements, strong heteroplasmy, or complex plant and fungal mitochondria.

**Run.**

```bash
cli/geneminer2 mito -f samples.tsv -o mito_output -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

**Notes.** Circularity requires agreement among overlap, local read-graph paths, and junction-spanning reads; reference coordinates never fill gaps. See the [Mitochondrial chapter](docs/mitochondria_EN.md).

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

### Navigate by task

- **Recover target sequences:** start with [Filter](docs/filter_EN.md) and [Assembler](docs/assembler_EN.md).
- **Interpret read-level or population variation:** see [Profiling](docs/profiling_EN.md) and [Population](docs/population_EN.md).
- **Specialised workflows:** mitochondrial, nuclear gene-family, and repeatome routes are listed below.

| Topic | Document |
| --- | --- |
| Commands and dependencies | [Command-line guide](manual/EN_US/command_line.md) |
| Output files | [Output-file guide](manual/EN_US/output.md) |
| Read filtering | [Filter chapter](docs/filter_EN.md) |
| Assemblers | [Assembler chapter](docs/assembler_EN.md) |
| Mitochondria | [Mitochondrial chapter](docs/mitochondria_EN.md) |
| Marker profiling | [Profiling chapter](docs/profiling_EN.md) |
| UCE population genetics | [Population chapter](docs/population_EN.md) |
| Gene subcommands | [Gene chapter](docs/gene_EN.md) |
| TE / repeatome | [TE / repeatome chapter](docs/te_EN.md) |
| Release history | [CHANGELOG](CHANGELOG.md) |

## Citation and contact

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [https://doi.org/10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

Related tools include [Easy353](https://doi.org/10.1093/molbev/msac261) and [GeneMiner](https://doi.org/10.1111/1755-0998.13924). When using `--alignment-filter alifilter`, also cite [AliFilter](https://doi.org/10.1093/molbev/msag097). Contact the maintainer at [xf@g.ecc.u-tokyo.ac.jp](mailto:xf@g.ecc.u-tokyo.ac.jp).

## Copyright and citation

Unless a file states otherwise, copyright for the original Rust implementation, documentation, tests, and workflow extensions in the current and later releases belongs to Fei Xia and is released under [GPL-3.0-or-later](LICENSE). Earlier releases made under the MIT License remain available under their original terms. This project retains the original copyright and license notices of third-party or directly ported code; see [NOTICE](NOTICE) for provenance boundaries. When publishing work that uses this software or its results, please cite this project and the original papers for relevant methods and software.
