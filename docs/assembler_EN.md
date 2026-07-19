# 2. Assembler

[中文版本](assembler_ZH.md)

Assembly turns recruited locus reads into contigs. Choose the workflow by biological target, not by sequencer or library label.

| Workflow | Appropriate target | Default backend | Result |
|---|---|---|---|
| `--assembly-mode original` | Exons, SCOs, nuclear markers, and mitochondrial markers | `original-rust` | Reference-guided contigs, trimmed by default |
| `--assembly-mode uce` | UCEs from genome-skimming or target-capture data | `uce-rust` | UCE cores plus read-supported flanks |

## Backends

`original-rust` is the deterministic Rust compatibility backend for `original`. Use `--assembler-implementation original` only for strict comparison with the fixed upstream GeneMiner2 Python implementation. `uce-rust` is the only backend for `uce`; it does not fall back to Python.

All backends use reference-positioned seeds, read-k-mer support, bidirectional extension, and read-slice validation.

## Original workflow

The original algorithm extends the highest-weight edge first and keeps alternatives on a stack for backtracking. It retains up to three candidates per side before combining and scoring them. This is appropriate for shorter reference-guided targets such as exons and organellar markers.

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original-rust` can reuse a versioned binary reference-k-mer cache with `--reuse-reference-cache`.

## UCE workflow and QC

UCE mode keeps linked mates during refiltering and skips default trimming. The `backbone` strategy makes one bounded look-ahead decision at each bubble, commits the winning branch, and does not backtrack. `--uce-rescue-reads` performs one additional raw-read recruitment and assembly round using accepted first-round contigs; failed or lower-density results revert.

Inspect `uce_assembly_summary.csv`, `uce_rescue_summary.csv`, `contigs_all_low/`, and downstream alignments. Important guardrails include unique-read density, supported breadth, largest unsupported gap, k-mer depth CV, and max/median depth ratio. Use `--uce-path-strategy search` only for sensitivity comparison with legacy branch enumeration.

File definitions are in the [output guide](../manual/EN_US/output.md); options are in the [command-line guide](../manual/EN_US/command_line.md).
