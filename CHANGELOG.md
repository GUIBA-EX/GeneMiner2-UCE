# Changelog

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
