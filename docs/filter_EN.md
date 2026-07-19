# 1. Filter

[中文版本](filter_ZH.md)

Filtering is the shared read-recruitment layer for the other three workflows. It scans reads for reference k-mers and writes locus-assigned reads; it does not itself decide contig validity, marker abundance, or population structure.

## Inputs and reference layout

The sample table is tab-delimited: sample name, R1, and optional R2. For `original` and `uce`, `-r` is a directory containing one FASTA per locus; the file stem is the locus name. Profiling is the exception: it accepts one `.fa` or `.fasta` marker library directly.

## First recruitment and refiltering

`filter` recruits reads with `-kf` and scan step `-s`, writing temporary assignments to `filtered_pe/`. `refilter` resolves per-locus assignment and applies depth and file-size limits, writing retained reads to `filtered/`.

In `uce` mode, a paired fragment is retained when either mate passes, so a core-mapping mate can retain flank evidence from its partner. Profiling performs only the first recruitment; it does not run `refilter`.

## Cache and reproducibility

Use `--reuse-reference-cache` for unchanged references. The cache fingerprint includes reference files, filter k-mer, and scan step. `--reference-cache-dir` selects a shared location; otherwise the cache is under `output/.gm2_reference_cache`.

## Outputs and checks

- `ref_reads_count_dict.txt` records first-pass read counts per locus.
- `filtered_pe/` contains first-pass assignments and is normally removed after successful refiltering.
- `filtered/` contains reads passed to assembly.
- `large_files/` appears only when refilter depth or size limits are exceeded.

A locus with no recruited reads is not evidence of absence by itself: check reference divergence, `-kf`, sequencing depth, and the read-count table. File definitions are in the [output guide](../manual/EN_US/output.md); option definitions are in the [command-line guide](../manual/EN_US/command_line.md).
