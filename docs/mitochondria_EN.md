# Mitochondrial workflow: ordinary circular animal mitochondria

`mito` uses the existing MainFilter, refilter, and Rust UCE assembler to recruit reads from an annotated GenBank reference. It reports a circular sequence only when the sample reads support it.

## Scope

This workflow is for **ordinary single circular animal mitochondrial genomes**: gene order should remain broadly comparable with the reference, without multipartite structure, major rearrangement, or heteroplasmy that needs phased reporting. It is not intended for complex plant or fungal mitochondria, fragmented or multichromosomal mitochondria, major rearrangements, or studies requiring specialised heteroplasmy or NUMT handling.

## Minimal command

```bash
cli/geneminer2 mito \
  -f samples.tsv \
  -o mito_output \
  -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

`--mito-max-reads 320` caps adaptive input; `--no-mito-adaptive-stop` switches to one-pass filtering.

## Workflow

```text
GenBank gene/rRNA/tRNA + genome/tile baits
→ MainFilter paired-read recruitment
→ one mitochondrial read pool
→ refilter → Rust UCE assembler
→ overlap, mate-link, and local read-graph joins
→ junction-spanning read validation of circularity
```

All baits are written as one mitochondrial locus. The reference is used only for recruitment and seeding: final sequence is neither coordinate-stitched nor reference-filled. Mate links propose adjacency and orientation only; gap bases must be recovered from a unique path in the same filtered read pool. Unresolved gaps remain broken and are never filled with `N`.

## Success criteria

A circular result must have one component, no `N`, every join supported by an overlap or unique read-graph bridge, a closed terminus, and at least `--mito-min-junction-support` true reads spanning the final junction. Otherwise partial output is retained but the command fails.

## Common and expert parameters

- `--mito-genbank`: required annotated mitochondrial GenBank reference.
- `--mito-max-reads 320`: at most approximately 1.05M paired-read blocks per adaptive stage; the workflow stops early when two successive stages return the same circular sequence.
- `--no-mito-adaptive-stop`: disable staged early stopping and use the normal one-pass `--max-reads` behaviour.

The following hidden expert overrides should be changed only to diagnose a known recruitment, graph-joining, or circularity problem: `--mito-min-overlap`, `--mito-min-overlap-identity`, `--mito-terminal-window`, `--mito-link-kmer`, `--mito-min-link-hits`, `--mito-min-pair-support`, `--mito-bridge-kmer`, `--mito-bridge-min-depth`, `--mito-max-bridge`, and `--mito-min-junction-support`.

## Outputs

- `<sample>/mito/mitochondrial_assembly.fasta`: circular sequence or partial components.
- `<sample>/mito/mitochondrial_assembly_summary.tsv`: status, joins, and junction support.
- `<sample>/mito/mitochondrial_mate_links.tsv`: accepted read-supported links.
- `.gm2_mito_reference/metadata/mitochondrial_genes.tsv`: bait metadata in 0-based half-open coordinates; `segments_0_half_open` preserves every segment of cross-origin or `join(...)` features.
