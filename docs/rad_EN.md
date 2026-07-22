# RAD: supplement an existing ipyrad loci matrix with WGS

`rad` supplements new WGS samples missing from an ipyrad `.loci` matrix. It does not invent RAD alleles from absent RAD reads.

| Start with | Get | Boundary |
| --- | --- | --- |
| Completed `.loci` or demultiplexed paired RAD reads + new WGS reads | Independent R1/R2 arms, recovery state, and a validated strict matrix | Does not infer the R1/R2 insert or directly prove allele dropout |

## Recommended route

```bash
# Build a reusable probe from a completed ipyrad matrix
cli/geneminer2 rad-probe --ipyrad-loci assembly.loci -o rad_probe

# Recover independent arms for new WGS samples
cli/geneminer2 rad --rad-probe rad_probe/rad_reference \
  -f wgs_samples.tsv -o rad_out -p 8

# Write a separate validated strict matrix
cli/geneminer2 rad-validate --rad-probe rad_probe/rad_reference \
  --rad-recovery rad_out/rad_recovery -o rad_validate_out
```

## Inputs and probe construction

`rad-probe` supports these input routes:

- `--ipyrad-loci FILE`: a completed `.loci` file; this is the recommended route.
- `--ipyrad-params FILE`: run ipyrad first. The default is steps `1--7`, then `<project_dir>/<assembly>_outfiles/<assembly>.loci`; use `--ipyrad-executable` or `--ipyrad-steps` when needed. `--ipyrad-loci` may additionally name a relocated ipyrad output.
- `--rad-denovo -f paired_rad_samples.tsv`: already-demultiplexed paired RAD reads.

`rad-probe` does not replace demultiplexing, restriction-enzyme recognition, adapter handling, or full within/between-sample clustering. Those belong in the ipyrad setup. The native de novo mode is a conservative candidate-probe builder, not an ipyrad-equivalent RAD matrix.

```text
paired_rad_samples.tsv
sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz
```

The native de novo builder counts per-sample canonical k-mers, retains solid k-mers, groups paired arms by minimizer, and emits only supported paired arms. It retains at most 128 deterministically sampled read pairs per sample-stack for consensus while recording total supporting pairs. Defaults are `k=31`, minimum solid depth `3`, support from `2` samples, and arm length `60 bp`. `--rad-overhang` checks R1 only; pass `--rad-overhang-r2` only when R2 starts at a known second restriction end.

## Recovery and validation

Both probe types use MainFilter, refilter, and `original-rust` assembly. `rad` accepts new WGS sample IDs only: a sample already present in the RAD reference, mismatched R1/R2 sample sets, duplicate normalized names, or malformed arm FASTA are rejected rather than silently merged.

`--rad-min-arm-breadth` defaults to `0.80`. `rad` creates a recovery matrix; run `rad-validate` afterwards. It compares each candidate arm against its own locus's multi-allelic bait and against every other locus's bait on the same arm. Defaults require query and reference breadth of at least `0.80`, identity of at least `0.90`, and an own-locus score at least 5% above the best foreign locus. A WGS sample enters the strict matrix only when both R1 and R2 pass. Tune the checks with `--rad-validate-min-breadth`, `--rad-validate-min-identity`, and `--rad-validate-min-delta`.

## Outputs and interpretation

- `rad_reference/arms/`: multi-allelic R1/R2 bait FASTA per locus.
- `rad_recovery/`: per-WGS-sample recruitment, refiltering, and assembly outputs.
- `rad_matrix/rad_sample_locus.tsv`: arm and joint recovery status per sample × locus.
- `rad_matrix/recovered_arms/`: unaligned FASTA retaining individually supported arms.
- `rad_matrix/paired_arms/`: unaligned FASTA for WGS samples with both recovered arms; not yet validated.
- `rad_validated/rad_validation.tsv`: per sample × locus × arm validation metrics and status.
- `rad_validated/strict_arms/`: original RAD baits plus only WGS samples whose two arms validate.

R1 and R2 are always independent observations. The workflow never bridges their unsequenced insert. All arm FASTA are unaligned and should be aligned only after choosing an appropriate missing-data strategy. `rad-validate` does not recruit or assemble reads again, and it leaves `rad_recovery/` and `rad_matrix/` unchanged.

`rad_missing_wgs_recovered` means that a locus absent from the input RAD matrix was recovered from WGS. It is not, by itself, evidence for restriction-site allele dropout; that interpretation additionally needs restriction-site and cross-locus WGS evidence.
