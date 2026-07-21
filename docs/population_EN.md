# 4. Population

[中文版本](population_ZH.md)

`population` creates one cohort coordinate system for multiple diploid UCE samples, then produces a joint SNP set, PCA, and optional ADMIXTURE results. It reports unphased genotypes; it is not a phasing or per-locus gene-tree workflow.

## Input and quick start

Each sample needs its completed UCE assembly (`uce_assembly_summary.csv` and accepted contigs in `results/`) plus the original R1/R2 files listed in the sample table. For `population`, optional `population` and `batch` fields may follow R2; leave the R2 field empty for single-end samples with metadata. The workflow requires minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional.

```bash
cli/geneminer2 population \
  -f samples.tsv -r baits_by_locus -o output -p 32 \
  --assembly-mode uce --engine panrefv2
```

`baits_by_locus` contains one FASTA per locus. Use `pseudoref` for accepted-contig representative references, or `--population-reference-fasta FILE` with `pseudoref` for a fixed external reference.

## Reference engines

| Engine | Purpose |
| --- | --- |
| `pseudoref` | Default accepted-contig representative reference; supports `sqcl-longest` and `supported`. |
| `panref` | Legacy experimental implementation. |
| `panrefv2` | Recommended streaming local-graph UCE reference builder. |

PanRefV2 builds a frozen paired-read core per locus, then admits only bounded single-mate rescue reads connected to that core. It uses bait minimizers, writes no per-locus recruited FASTQ, and never scaffolds across UCE loci.

Backbone order: bait agreement > supporting samples > current graph-edge read/PE support > per-sample depth stability > length.

## Stages and outputs

```text
reference → mapping → calling → selection
```

The stages build the reference, map/QC reads, jointly call/filter variants, then create all-SNP, one-SNP-per-UCE, and LD-pruned panels with PCA and optional ADMIXTURE. Use `--population-start-at` and `--population-stop-after` only with validated outputs.

PanRefV2 writes `population/reference/panrefv2/`:

- `index_metadata.tsv` — bait minimizer summary;
- `recruitment_summary.tsv` — strong, candidate rescue, accepted rescue, and ambiguous pairs;
- `population_graph.gfa` — locus-partitioned local graphs;
- `locus_summary.tsv` — locus status and backbone evidence.

Only `pass` loci enter `population_reference.fasta` by default. `short`, `low_sample_support`, `low_coverage`, `complex`, and `no_core` remain in the QC report. Use `--population-panrefv2-include-low-confidence` only for diagnostics.

## QC

Default filters are MAPQ 20, base quality 20, genotype DP 5, GQ 20, site QUAL 20, call rate 0.8, and MAC 2. Check PanRefV2 locus status, mapping QC, variant QC, panel counts, PCA agreement across the three panels, and ADMIXTURE CV error before interpretation. Individual QC is written to `structure/qc/individuals.imiss`, `individuals.het`, and `individuals.genome` (pairwise PI_HAT); use `--skip-relatedness-qc` to omit the quadratic relatedness calculation for large cohorts.
