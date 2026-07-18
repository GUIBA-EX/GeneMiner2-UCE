# UCE workflow

[中文版本](uce-workflow_ZH.md)

This guide describes `--assembly-mode uce`, including path selection, read support, raw-read rescue, quality guardrails, and reference caching. See the [command-line guide](../manual/EN_US/command_line.md) for all options and the [output guide](../manual/EN_US/output.md) for file fields.

## Default workflow

Without explicit subcommands, UCE mode runs:

```text
filter → refilter → assemble → combine → tree
```

It skips `trim` so newly recovered flanks are not cut back to the probe interval. Add `trim` explicitly when reference-based trimming is required.

During refiltering, a paired-end fragment is retained whenever either mate passes the locus filter. One mate can therefore anchor in the conserved core while the other contributes flank information.

## Backbone assembly

The default `--uce-path-strategy backbone` constructs one path on each side of a seed. At a bubble, it performs a look-ahead bounded by `--uce-backbone-lookahead`, prefers the branch that remains extendable for longer, and breaks ties with cumulative k-mer support.

Once a branch is selected, sibling edges are discarded. The algorithm keeps no branch stack and does not repeatedly backtrack. Visited k-mers cannot be re-entered, so cycles terminate. Non-branching chains are consumed as unitigs to reduce per-node decisions.

For comparison with legacy branch enumeration, use:

```bash
--uce-path-strategy search \
  --uce-side-candidates 8
```

`--uce-side-candidates` affects only the `search` strategy.

## Candidate support and guardrails

Candidates are not ranked by length alone. The assembler records:

- total, uniquely placed, and multi-mapping reads;
- the union and breadth of read-supported bases;
- the largest unsupported gap;
- unique-read density;
- median k-mer depth, depth CV, and max/median depth ratio.

UCE mode has no default contig-length limit. Candidates at least 1000 bp long must satisfy:

```text
unique_read_density = uniquely_placed_read_count / contig_length
unique_read_density >= 0.003
```

Relevant options are:

```bash
--uce-max-contig-length 0 \
  --uce-min-read-density 0.003 \
  --uce-density-check-min-length 1000 \
  --uce-max-depth-cv 0 \
  --uce-max-depth-ratio 0
```

The length, depth-CV, and depth-ratio limits are disabled by default. Enable them only when abnormal extensions, repeat peaks, or strongly uneven depth are evident in the data.

## Relaxed settings for short probes

For short baits or moderately divergent samples, test:

```bash
-sb unlimited \
  -ka 0 \
  --min-ka 17 \
  --max-ka 31 \
  -e 1
```

These settings relax boundaries and allow lower assembly k-mers, but may also admit weak or erroneous candidates. Inspect `uce_assembly_summary.csv`, downstream alignments, and stability across parameter choices.

## Raw-read rescue

`--uce-rescue-reads` performs one additional recruitment round:

1. Combine the original locus reference with an accepted first-round contig.
2. Recruit matching reads again from the raw inputs.
3. Repeat refiltering and assembly.
4. Compare the first-round and rescue results and retain or revert.

First-round contigs shorter than `--uce-rescue-min-contig-length` (default 60) are not used as rescue references. Rescue processes at most four samples concurrently, with at most four threads per sample and an overall limit set by `-p`.

If the first round was accepted but rescue is missing or rejected, the first-round result is restored as `reverted_failed_rescue`. When both rounds pass, unique-read densities are compared:

```text
before_density = before_unique_read_count / before_contig_length
rescue_density = rescue_unique_read_count / rescue_contig_length
density_ratio = rescue_density / before_density
```

By default, `density_ratio < 0.5` restores the first-round result as `reverted_density_drop`. Change the threshold with `--uce-rescue-min-density-ratio`. The `after_*` columns in `uce_rescue_summary.csv` always describe the rescue attempt; after a reversion, the final FASTA still contains the first-round contig.

## Reference cache

Use `--reuse-reference-cache` when repeatedly processing the same references. The default cache is `output/.gm2_reference_cache`; `--reference-cache-dir` selects a shared directory.

The fingerprint includes reference names, sizes, modification times, the filter k-mer, and the scan step. Caching reduces index construction time without changing filtering or assembly. Rescue references contain sample-specific first-round contigs, so rescue indices are still built per sample.

## Recommended checks

Before using UCE results, inspect:

- acceptance, rejection reasons, and read support in `uce_assembly_summary.csv`;
- `uce_rescue_summary.csv` when rescue is enabled;
- rejected candidates in `contigs_all_low/` for abnormal extension or weak support;
- missingness, unusually long branches, and possible paralogs in merged alignments;
- sample- and locus-level recovery reported by `stats`.
