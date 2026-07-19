# 3. Profiling

[中文版本](profiling_ZH.md)

Profiling is an assembly-free workflow for recovering and quantifying any amplicon marker from WGS or metagenomic reads. It reports relative marker-group signal, not contigs and not organism or cell proportions.

## Workflow

```text
one GeneMiner2 k-mer recruitment → Themisto pseudoalignment → mSWEEP group estimation
```

It runs neither `refilter`, `assemble`, `combine`, nor `tree`. `--profile-kmer-size` sets the same odd k-mer size, from 15 to 31, for both recruitment and Themisto.

## Required inputs

- A single `.fa` or `.fasta` marker library passed directly to `-r`.
- A two-column `--profile-group-map` TSV: `reference_id<TAB>group`.
- A sample table with WGS or metagenomic reads.

The reference ID is the first whitespace-delimited FASTA-header field. Every reference maps to exactly one group; duplicate rows are permitted only when they assign the same group. `themisto` and `mSWEEP` must be on `PATH` or supplied with `--profile-themisto` and `--profile-msweep`. An optional `--profile-decoy` FASTA represents plausible non-target sequence.

## Run and cache

```bash
cli/geneminer2 profiling \
  -f samples.tsv -r marker_reference.fasta \
  --profile-group-map marker_groups.tsv \
  -o output -p 8 --profile-decoy non_target.fasta
```

Themisto indices use a content-addressed cache. Use `--profile-index-dir` to share it and `--profile-force-rebuild` only when rebuilding is intended.

## Outputs and interpretation

Each sample writes `marker_profile/`:

- `marker_group_abundance.tsv`: groups, evidence, detection state, and relative signal.
- `marker_qc.tsv`: pseudoalignment, target/decoy, evidence, mSWEEP, and parameter counts.
- `marker_reference_metadata.tsv`: reference IDs, Themisto colors, and groups.

`evidence_queries` and `exclusive_queries` count individual FASTA/FASTQ query records, not paired fragments. `relative_proportion` is renormalized after the exclusive-evidence rule and remains an uncalibrated marker-signal proportion.

## QC and calibration

Check pseudoaligned queries, target-decoy shared evidence, decoy-exclusive evidence, and each group’s exclusive support before calling a detection. Set `--profile-min-evidence` using negative controls, mixtures, and downsampling matched to the marker and sequencing depth. Keep the reference library, group map, k-mer size, thresholds, and decoy strategy fixed for cross-sample comparisons.

See [Filter](filter_EN.md), the [output guide](../manual/EN_US/output.md), and the [command-line guide](../manual/EN_US/command_line.md).
