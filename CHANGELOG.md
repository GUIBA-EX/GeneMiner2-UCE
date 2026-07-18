# Changelog

## v0.7 — ITS2 assembly and Rust utility migration

- Added ITS2 multi-candidate assembly with paired-fragment compatibility, equivalence groups, diagnostic support, and EM abundance estimates; ITS2 now remains strictly Rust-only on failure.
- Reimplemented alignment cleanup, sequence merging, reference trimming, and UCE statistics as readable Rust utilities while preserving their command-line contracts.
- Removed the unused MUSCLE integration and obsolete validation-only helper scripts.
- Fixed single-end statistics, rescue scheduling after sample failures, deterministic gene-tree ordering, and top-level CLI error handling.
- Synchronized Chinese and English README, command-line, and output documentation with the current CLI; removed obsolete GUI-era console output and local debug artifacts.

## v0.65 — MainFilter deterministic lookup optimization

- Optimized the Rust primary filter's short-k-mer scan with a DNA lookup table, modulo-free probe scheduling, and `AHashMap` k-mer lookup.
- Kept filtering semantics, command-line options, and cache/output formats unchanged; documented byte-level output verification and the decisions not to add threads, LRU output handles, or low-gain hash alternatives.

## v0.6 — Scalable UCE graph assembly

- Stream filtered reads in bounded batches and count k-mers through parallel, sorted per-batch aggregation.
- Compress non-branching UCE backbone paths into unitigs, retaining bounded decisions only at graph junctions.
- Add optional compact GFA and DOT assembly-graph output via `--assembler-graph-format`.
- Add `--assembler-read-chunk-size` and `--assembler-kmer-count-threads`, while preserving the unmodified Python fallback.

## v0.5 — Rust UCE assembly and reusable population analysis

- Added the high-performance Rust UCE assembler with compact rolling k-mers, a bounded non-backtracking backbone path strategy, reference caching, and parallel per-locus assembly.
- Made Rust assembly the default via `--assembler-implementation auto`; failed or unavailable Rust runs now clean incomplete outputs and retry the unmodified Git-baseline Python assembler.
- Retained strict `rust` and direct `original` assembler modes for reproducibility and diagnosis.
- Added fixed external cohort-reference support, checked resume stages (`mapping`, `calling`, and `selection`), and per-stage variant-count QC to the population workflow.
- Added reusable tools for summarizing UCE validation runs and comparing a population pseudo-reference with an external reference.
- Updated Chinese and English command-line/output documentation and regression coverage.

## v0.4 — Population analysis

- Added the Rust `population` workflow: cohort-reference construction, uniform minibwa mapping, joint bcftools variant calling, and one representative SNP per UCE.
- Added SqCL-inspired longest-eligible-contig reference selection, with a read-support-first alternative and per-sample reference-contribution diagnostics.
- Added all-SNP, one-SNP-per-UCE, and LD-pruned VCF/PLINK panels with PCA for each panel.
- Added automated ADMIXTURE K-range analysis, cross-validation summaries, status reporting, and retained logs.
- Added mapping-rate, coverage-breadth, depth, sample-name, and reference-provenance quality-control reports.
- Added real-tool integration tests for minibwa, samtools, bcftools, PLINK, and ADMIXTURE.
- Updated Chinese and English usage and output documentation.

## v0.3 — Rust primary filter

- Reimplemented the primary read filter in Rust while retaining the original command-line and cache compatibility behavior.

## v0.2 — UCE assembly validation and rescue

- Added UCE assembly guardrails, read-support validation, and controlled raw-read rescue.

## v0.1 — UCE workflow foundations

- Added sequence-integrity fixes and the initial UCE-focused command-line workflow.
