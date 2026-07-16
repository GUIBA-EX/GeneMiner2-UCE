# GeneMiner2 UCE CLI

**[查看中文版说明](README_zh_cn.md)**

This repository is a command-line focused fork of GeneMiner2 for target-enrichment and UCE workflows. GUI project files, graphical interface documentation, screenshots, and bundled demo datasets have been removed so the repository only contains the CLI source, build files, and command-line documentation.

This fork is not distributed as an installable Python package. Build the standalone CLI with `make`, then run `cli/geneminer2`. After every `git pull`, rebuild with `make`; otherwise the generated executable may still contain older code.

## Main Features

- Reference-guided marker recovery from next-generation sequencing reads.
- UCE-oriented assembly mode through `--assembly-mode uce`, designed to retain longer read-supported flanking contigs.
- Optional one-round UCE raw-read rescue through `--uce-rescue-reads`.
- phyluce-compatible UCE contig export under `uce_contigs/`.
- HybPiper-style UCE statistics through the `stats` subcommand.
- Optional AliFilter alignment filtering through `--alignment-filter alifilter`.
- Controlled combine-stage parallelism through `--msa-threads` and `--filter-processes`.

![GeneMiner2 UCE CLI workflow](docs/images/summary_EN.png)

## Changes in This Fork

This fork keeps the original GeneMiner2 reference-guided recovery model, but adjusts the command-line workflow for UCE loci, where the bait/probe sequence can be short and the useful phylogenetic signal often comes from flanking sequence.

### CLI-only layout

The GUI project, screenshots, bundled demo data, and large historical files have been removed from the active repository. The project is intentionally kept as a small source tree plus `Makefile` build. There is no Python `console_scripts` entry point and no `pyproject.toml`; use the generated `cli/geneminer2` executable after building.

### Rust re-filter implementation

The secondary read filter now has a Rust implementation under `rust/main_refilter_new/`. It is intended as a drop-in replacement for `scripts/main_refilter_new.py`: the command-line options and output layout are kept compatible, including `--keep-linked-mates` for UCE workflows.

The Python implementation is still kept in `scripts/main_refilter_new.py` as a readable reference implementation and fallback. During `make`, GeneMiner2 builds the Rust binary when `cargo` is available; if Cargo is not installed, the build falls back to packaging the Python implementation with PyInstaller.

### UCE assembly mode

`--assembly-mode uce` changes the assembly behavior so GeneMiner2 is less likely to trim contigs back to the short reference/probe interval. In UCE mode, the assembler prefers longer candidates that still have read support, penalizes weakly supported over-extension using read-density and k-mer-depth continuity, and the default command set skips the reference-based `trim` step unless `trim` is requested explicitly.

Recommended UCE-oriented assembly options are:

```bash
--assembly-mode uce \
-sb unlimited \
-ka 0 \
--min-ka 17 \
--max-ka 31 \
-e 1
```

These settings keep boundary trimming permissive, allow automatic assembly k-mer selection over a lower range, and reduce the k-mer count threshold. They are designed for short UCE baits and divergent samples, but they can also admit noisier candidates, so the rescue summary and downstream alignments should still be inspected.

### Reference/index cache

Repeated runs against the same reference directory can reuse reference k-mer indexes:

```bash
--reuse-reference-cache
```

The cache is fingerprinted by reference file names, sizes, modification times, `-kf`, and `-s`. By default it is written under `output/.gm2_reference_cache`; use `--reference-cache-dir` to place it in a shared project or scratch directory. This speeds up repeated filter and assembly runs, but it does not change contig selection or improve assembly quality. UCE rescue references are rebuilt per sample because they include sample-specific preliminary contigs.

UCE mode also applies conservative contig guardrails before candidate selection. By default, contigs longer than 5000 bp are rejected, and contigs at least 1000 bp long must have `uniquely_placed_read_count / contig_length >= 0.003`. Repetitively placed read slices remain visible in the summary but do not provide positional support. These defaults are intended to suppress very long, weakly supported rescue artifacts while leaving rejected candidates available for inspection. They can be adjusted with:

```bash
--uce-max-contig-length 5000 \
--uce-min-read-density 0.003 \
--uce-density-check-min-length 1000
```

Two optional advanced guardrails are available but disabled by default:

```bash
--uce-max-depth-cv 0 \
--uce-max-depth-ratio 0
```

Set these to positive values only when you want to reject candidates with highly uneven k-mer depth or strong repeat-like depth spikes.

### Paired-end mate retention

In UCE mode, the re-filtering step keeps a paired-end read pair when either mate passes the locus filter. This is important for short probes: one mate may overlap the conserved UCE core while the other extends into the flanking region. Keeping the pair gives the assembler more information for contig extension.

### One-round raw-read rescue

`--uce-rescue-reads` runs one additional recruitment round after the first assembly:

1. Build temporary rescue references from the original locus reference plus the first-round contig.
2. Re-filter raw reads against these rescue references.
3. Re-run re-filtering and assembly using the rescue references.
4. Compare the rescue result with the first-round result.

The rescue stage uses controlled parallelism: up to four samples are rescued at the same time, with up to four threads per sample. This avoids launching too many independent read-filtering jobs when many samples are present.

### Density-ratio rollback

Raw-read rescue can occasionally create very long but weakly supported contigs. To avoid accepting those artifacts, this fork compares read density before and after rescue:

```text
before_density = before_read_count / before_contig_length
rescue_density = rescue_read_count / rescue_contig_length
density_ratio = rescue_density / before_density
```

By default, a rescue result is rejected only when:

```text
density_ratio < 0.5
```

Rejected or missing rescue results are restored to the first-round contig and marked as `reverted_failed_rescue`; accepted rescue results whose unique-read density drops below the threshold are restored and marked as `reverted_density_drop`. The threshold can be changed with:

```bash
--uce-rescue-min-density-ratio 0.5
```

The rescue summary records `before_read_density`, `after_read_density`, and `density_ratio` using uniquely placed reads when available, so the reason for rollback is visible. The final accepted sequence may be the first-round contig if rollback occurred.

### UCE and phyluce outputs

When `--assembly-mode uce` is used, the workflow writes:

- `uce_assembly_summary.csv`: per-sample and per-locus acceptance status and rejection reason, selected contig length, legacy left-to-right support span, union of uniquely placed slice-supported bases, support breadth and maximum unsupported gap, total/unique/multi-mapping read counts, density metrics, k-mer depth metrics, candidate count, and low-quality flag.
- `uce_rescue_summary.csv`: rescue before/after comparison, density ratio, rollback status, and errors.
- `uce_contigs/`: phyluce-compatible per-sample contig FASTA files.
- `contigs_all_low/`: rejected UCE candidates retained for inspection but never promoted to primary results, rescue references, combined matrices, or phyluce exports.

After a UCE run, the `stats` subcommand summarizes recovery across samples and loci:

```bash
cli/geneminer2 stats \
  -f samples.tsv \
  -r references \
  -o output \
  --stats-no-heatmap
```

It writes `uce_stats.tsv`, `uce_locus_stats.tsv`, `uce_seq_lengths.tsv`, `uce_read_counts.tsv`, and `uce_filtered_read_counts.tsv`. If `pandas`, `seaborn`, and `matplotlib` are available and `--stats-no-heatmap` is not supplied, it also writes `uce_recovery_heatmap.png` and `uce_read_counts_heatmap.png`.

### AliFilter integration

`--alignment-filter alifilter` can be used during the combine stage as an alternative to trimAl. This is useful when many UCE alignments include noisy or sparsely occupied columns. AliFilter must be available in `PATH`; it is not bundled in this repository. Omit `--alifilter-model`, or set it to `default`, to use AliFilter's built-in default model; pass a real `model.json` path only when using a custom model.

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
