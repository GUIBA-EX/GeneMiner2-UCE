# GeneMiner2-UCE Command-Line Guide

GeneMiner2-UCE is a command-line workflow for UCEs and other phylogenetic markers. This repository does not include the former GUI, bundled demonstration data, or legacy graphical documentation.

## 1. Building from source

### 1.1 Build dependencies

A complete build requires:

- a C/C++ compiler and [zlib](https://zlib.net/);
- [Rust and Cargo](https://www.rust-lang.org/tools/install/);
- Python 3.11, Cython, PyInstaller, Biopython, NumPy, SciPy, pandas, matplotlib, and seaborn.

Rust/Cargo is required for the complete current build, including the primary read filter, secondary filter, and population workflow. The Haxe source remains available only as a compatible implementation of the primary filter; it does not replace the complete Rust build.

On Ubuntu, first install the system dependencies:

```bash
sudo apt install build-essential zlib1g zlib1g-dev
```

Then create the Python environment:

```bash
conda create -c conda-forge -n geneminer \
  python=3.11 numpy=2.1.3 biopython cython matplotlib \
  pandas seaborn pyinstaller scipy setuptools wheel
conda activate geneminer
```

### 1.2 Runtime dependencies

The phylogenetic workflow calls a subset of the following external programs, depending on the selected stages:

```bash
conda install -c bioconda \
  aster blast clustalo fasttree iqtree mafft magicblast minimap2 \
  muscle raxml-ng trimal veryfasttree
```

Optional and population-specific tools:

- `--alignment-filter alifilter` requires a separate AliFilter installation with `AliFilter` in `PATH`.
- `population` requires `minibwa`, samtools, bcftools, and PLINK 1.9.
- ADMIXTURE performs ancestry analysis. If it is absent, population analysis still completes cohort-reference construction, VCF generation, PLINK export, and PCA, while recording its status as `unavailable`.

Population executables are located through `PATH` by default. Override them individually with `--population-minibwa`, `--population-samtools`, `--population-bcftools`, `--population-plink`, and `--population-admixture`.

### 1.3 Download and build

```bash
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2-UCE.git
cd GeneMiner2-UCE
make
```

The generated entry point is:

```bash
cli/geneminer2 -h
```

Run `make` again after updating the source.

## 2. Input files

Every command requires:

- `-f FILE`: tab-delimited sample table;
- `-r DIR`: reference-sequence directory;
- `-o DIR`: output directory.

### 2.1 Sample table

Each row represents one sample. Single-end data use two columns; paired-end data use three:

```text
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

Absolute paths are recommended. Sample names must be unique. Population analysis writes a manifest linking the original names to internal directory and VCF sample IDs.

### 2.2 Reference directory

Use one FASTA file per locus; the file name becomes the locus name. A file may contain one or several reference sequences:

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

## 3. Subcommands and default workflows

One or more subcommands can be listed in execution order:

| Subcommand | Function |
| --- | --- |
| `filter` | Recruit reads using reference k-mers |
| `refilter` | Refine per-locus read assignment and filtering |
| `assemble` | Assemble target sequences with the wDBG assembler |
| `population` | Build a cohort UCE reference and generate SNP, PCA, and ADMIXTURE results |
| `consensus` | Generate consensus sequences at heterozygous sites |
| `trim` | Remove flanks relative to reference sequences |
| `combine` | Merge samples, align loci, clean sequences, and filter alignment columns |
| `tree` | Infer coalescent or concatenated trees |
| `stats` | Summarize UCE recovery and optionally generate heatmaps |

When no subcommand is given:

- `--assembly-mode reference` (default) runs `filter refilter assemble trim combine tree`;
- `--assembly-mode uce` runs `filter refilter assemble combine tree`, omitting `trim` so newly recovered UCE flanks are not cut back to the reference interval.

Default reference-mode example:

```bash
cli/geneminer2 \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8
```

## 4. UCE assembly

UCE mode relaxes boundaries imposed by short probes and favors longer flanking sequences that retain read support. During refiltering, a paired-end read pair is retained whenever either mate passes the locus filter.

Basic example:

```bash
cli/geneminer2 \
  -f samples.tsv \
  -r references \
  -o output \
  -p 8 \
  --assembly-mode uce \
  --uce-rescue-reads
```

`--uce-rescue-reads` uses preliminary contigs plus the original references to recruit raw reads again, followed by one additional re-filtering and assembly round. Rescue processes at most four samples concurrently, with up to four threads per sample and an overall limit set by `-p`.

After relaxing `-sb`, `-e`, or the assembly k-mer range, inspect `uce_assembly_summary.csv`, `uce_rescue_summary.csv`, and downstream alignments. See the [output guide](output.md) for details.

## 5. Population-genetic analysis

### 5.1 Scope

`population` derives an unphased diploid SNP matrix from UCE assemblies and the original reads of multiple samples. Its principal uses are PCA, ADMIXTURE, ancestry comparison, and species delimitation.

Before running it, each sample must have completed UCE assembly, and the following must remain available:

- `uce_assembly_summary.csv`;
- accepted UCE contigs under `results/`;
- the original reads listed in the sample table.

This mode reports unphased genotypes, not two complete haplotype sequences. It does not replace phasing when haplotype sequences, recombination information, or per-locus gene trees are required.

### 5.2 Workflow

1. Pool accepted contigs by locus and build a cohort UCE reference.
2. Map every sample uniformly to the same reference with minibwa.
3. Jointly call variants with bcftools and apply genotype- and site-level filters.
4. Generate all-SNP, one-SNP-per-UCE, and LD-pruned panels.

PLINK runs PCA on all three panels. ADMIXTURE uses the one-SNP-per-UCE primary panel by default and evaluates cross-validation error across the requested K range.

### 5.3 Example

After UCE assembly is complete, run population analysis independently:

```bash
cli/geneminer2 population \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8 \
  --assembly-mode uce \
  --population-admixture-k-min 2 \
  --population-admixture-k-max 6
```

Before interpreting ancestry results, inspect:

- mapping rate, coverage breadth, and depth in `population/mapping/mapping_qc.tsv`;
- whether cohort-reference contributions are concentrated in a few samples in `population/reference/reference_contribution.tsv`;
- sample and site missingness;
- agreement among PCA results from the three SNP panels.

## 6. Staged-run examples

Rebuild only a concatenated tree from existing results:

```bash
cli/geneminer2 tree \
  -f samples.tsv -r references -o output \
  -m concatenation --phylo-program iqtree
```

Generate consensus sequences, trim them to references, and combine them:

```bash
cli/geneminer2 consensus trim combine \
  -f samples.tsv -r references -o output \
  -c 0.75 -ts consensus -tm all -tr 0.5 \
  -cs trimmed --msa-program muscle
```

Generate tables from existing UCE results without plotting heatmaps:

```bash
cli/geneminer2 stats \
  -f samples.tsv -r references -o output \
  --stats-no-heatmap
```

## 7. Parameter reference

The tables below list the main public options and current defaults. Run `cli/geneminer2 -h` for the complete help associated with the checked-out source.

### 7.1 General input and parallelism

| Option | Description |
| --- | --- |
| `-f FILE` | Sample table; required |
| `-r DIR` | Reference-sequence directory; required |
| `-o DIR` | Output directory; required |
| `-p INT` | Total parallel-process limit; default `1` |

### 7.2 Read filtering and refiltering

| Option | Description |
| --- | --- |
| `-kf INT` | Filter k-mer size; default `31` |
| `-s, --step-size INT` | Read-scanning step; default `4` |
| `--max-reads INT` | Maximum million reads processed per file; `0` means unlimited |
| `--reuse-reference-cache` | Reuse a fingerprinted reference k-mer index |
| `--reference-cache-dir DIR` | Reference-cache directory; default `output/.gm2_reference_cache`; requires the preceding option |
| `--depth-low-water-mark INT` | Below this depth, attempt relaxed read recruitment; default `50` |
| `--depth-limit INT` | Maximum depth processed during refiltering; default `768` |
| `--file-size-limit INT` | Refiltering file-size limit; default `6` |

### 7.3 Assembly and UCE options

| Option | Description |
| --- | --- |
| `-ka INT` | Assembly k-mer size; default `0` for automatic estimation |
| `--min-ka INT` / `--max-ka INT` | Automatic-estimation range; defaults `21` / `51` |
| `-e, --error-threshold INT` | k-mer error threshold; default `2` |
| `-sb, --soft-boundary VALUE` | Integer, `auto`, or `unlimited`; default `auto` |
| `-i, --search-depth INT` | Search depth; default `4096` |
| `--min-coverage INT` | Minimum contig read depth; default `0` |
| `--assembly-mode MODE` | `reference` or `uce`; default `reference` |
| `--uce-side-candidates INT` | One-sided branch candidates; default `8`, minimum `3` |
| `--uce-max-contig-length INT` | Maximum UCE contig length before scoring; default `5000`; `0` disables |
| `--uce-min-read-density FLOAT` | Minimum unique-read/length ratio for long contigs; default `0.003` |
| `--uce-density-check-min-length INT` | Minimum contig length for the density guardrail; default `1000` |
| `--uce-max-depth-cv FLOAT` | Maximum k-mer-depth CV; default `0` disables |
| `--uce-max-depth-ratio FLOAT` | Maximum/median k-mer-depth ratio; default `0` disables |
| `--uce-rescue-reads` | Perform one UCE raw-read rescue round |
| `--uce-rescue-min-contig-length INT` | Minimum rescue-reference contig; default `60` and never below `-kf` |
| `--uce-rescue-min-density-ratio FLOAT` | Minimum rescue/first-round density retained; default `0.5` |

### 7.4 Population options

| Option | Description |
| --- | --- |
| `--population-reference-strategy MODE` | `sqcl-longest` (default) or `supported` |
| `--population-min-mapq INT` / `--population-min-baseq INT` | Minimum MAPQ / base quality; defaults `20` / `20` |
| `--population-min-dp INT` / `--population-min-gq INT` | Set lower-quality genotypes to missing; defaults `5` / `20` |
| `--population-min-qual FLOAT` | Minimum site QUAL; default `20` |
| `--population-min-call-rate FLOAT` | Minimum non-missing genotype fraction; default `0.8` |
| `--population-min-mac INT` | Minimum minor allele count; default `2` |
| `--population-ld-window INT` / `--population-ld-step INT` | LD-pruning window and step; defaults `50` / `5` SNPs |
| `--population-ld-r2 FLOAT` | LD-pruning r² threshold; default `0.2` |
| `--population-admixture-k-min INT` / `--population-admixture-k-max INT` | ADMIXTURE K range; defaults `2` / `6`; maximum K is capped at sample count |
| `--population-admixture-cv INT` | ADMIXTURE CV folds; default `10`, capped at sample count |
| `--population-stop-after STAGE` | Stop after `reference`, `mapping`, `calling`, or `selection`; default `selection` |
| `--population-skip-mark-duplicates` | Skip samtools duplicate marking |
| `--population-skip-plink` | Omit PLINK, PCA, LD-pruned, and ADMIXTURE outputs |
| `--population-skip-admixture` | Generate PLINK and PCA without running ADMIXTURE |
| `--population-minibwa PATH` | minibwa executable; default `minibwa` |
| `--population-samtools PATH` | samtools executable; default `samtools` |
| `--population-bcftools PATH` | bcftools executable; default `bcftools` |
| `--population-plink PATH` | PLINK 1.9 executable; default `plink` |
| `--population-admixture PATH` | ADMIXTURE executable; default `admixture` |

### 7.5 Consensus, trimming, and combining

| Option | Description |
| --- | --- |
| `-c, --consensus-threshold FLOAT` | Consensus threshold; default `0.75` |
| `-ts, --trim-source SOURCE` | `assembly` or `consensus`; default is the preceding stage's output |
| `-tm, --trim-mode MODE` | `all`, `longest`, `terminal`, or `isoform`; default `terminal` |
| `-tr, --trim-retention FLOAT` | Reference-trimming retention fraction; default `0` |
| `-cs, --combine-source SOURCE` | `assembly`, `consensus`, or `trimmed`; default is the preceding stage's output |
| `-cd, --clean-difference FLOAT` | Maximum acceptable pairwise alignment difference; default `1.0` |
| `-cn, --clean-sequences INT` | Minimum sequences per alignment; default `0` |
| `--msa-program PROGRAM` | `clustalo`, `mafft`, or `muscle`; default `mafft` |
| `--msa-threads INT` | Threads per MSA job; default `1` and cannot exceed `-p` |
| `--alignment-filter PROGRAM` | `trimal`, `alifilter`, or `none`; default `trimal` |
| `--filter-processes INT` | Concurrent alignment-filter jobs; default equals `-p` |
| `--alifilter-model MODEL` | AliFilter model specification or `model.json` path |
| `--strict-combine-errors` | Stop when any locus fails alignment, cleanup, or filtering |
| `--no-alignment` | Skip multiple-sequence alignment |
| `--no-trimal` | Deprecated alias for `--alignment-filter none` |

### 7.6 Tree inference and statistics

| Option | Description |
| --- | --- |
| `-m, --tree-method METHOD` | `coalescent` or `concatenation`; default `coalescent` |
| `-b, --bootstrap INT` | Bootstrap replicates; default `1000` |
| `--phylo-program PROGRAM` | `raxmlng`, `iqtree`, `fasttree`, or `veryfasttree`; default `fasttree` |
| `--stats-no-heatmap` | Do not generate UCE statistics heatmaps |
| `--stats-count-input-reads` | Count input FASTQ reads for `InputReads` and `PctFiltered`; slow on large datasets |

Fractional thresholds use decimals from `0` to `1`. Before changing filters on a full dataset, run a small test with defaults and inspect intermediate VCFs, mapping QC, and locus occupancy.
