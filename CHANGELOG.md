# Changelog

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
