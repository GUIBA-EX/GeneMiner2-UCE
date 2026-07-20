# TE / repeatome mode

`geneminer2 te` is a reference-free repeatome workflow for genome-skimming or WGS short reads. It produces conservative repeat units, annotation evidence, and sample abundance; it is **not** complete-TE annotation, insertion-site calling, or a TE-phylogeny workflow.

| Start with | Get | Main boundary |
| --- | --- | --- |
| Genome-skimming or WGS short reads | EQ units, conservative annotation, and sample RPM | Does not report complete TEs or insertion sites |

## Quick start

```bash
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

TE uses its own sample manifest and does not require `-r`:

```text
taxon_id  sample_id  read1  read2
Taxon_A   A01        /data/A01_R1.fq.gz  /data/A01_R2.fq.gz
Taxon_B   B01        /data/B01.fq.gz
```

The fourth field, `read2`, is optional. `--te-read-ledger FILE` is optional and accepts `sample_id read_id` rows to exclude reads already assigned with high confidence to target sequences.

## Workflow

```text
discover → curate → annotate → quantify
```

- **discover** taxon-balances sampling and finds atomic seed groups.
- **curate** runs MainFilter once, retains candidate reads and paired-end linkage evidence, and builds an exact-equivalence (EQ) library.
- **annotate** builds up to three bounded, read-supported local fragments per EQ. Ambiguous local bubbles are retained as alternatives rather than merged; it then detects tandem/satellite features and performs conservative library-homology annotation.
- **quantify** maps every eligible input read pair back to the annotated fragment set (not only MainFilter recruits), reports RPM, mapped bases, and fragment mean depth; it never reruns MainFilter.

The default `--te-stage all` runs all four stages. `--te-stage discover|curate|annotate|quantify` reruns one stage. Downstream stages validate upstream manifests, so changed inputs or parameters require rerunning the appropriate upstream stage.

## Optional library annotation

```bash
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32 \
  --te-library curated_te.fa
```

`--te-library` is an optional classified TE FASTA. Headers use `name#Class/Subclass`, for example `RTE1#DNA/TcMar`. A high-confidence class requires enough fragment length and read-pair support, identity, coverage, and a score margin over the second class. Tune these conditions with `--te-annotate-min-*`.

Annotation never merges or removes an EQ. PE bridges and read transitions remain `linked_not_merged` evidence only; short, weak, or conflicting fragments remain `unknown_repeat` or `unknown_interspersed_repeat`. The built-in comparison is a conservative rolling-k-mer/diagonal screen, not a substitute for assembly-dependent EDTA or RepeatModeler.

## Outputs and interpretation

```text
te_out/
├── 01_discover/                 atomic catalog, seeds, manifest
├── 02_curate/                   EQ library, candidate reads, linkage, manifest
├── 03_annotate/
│   ├── fragments/EQ00001.fasta  read-supported fragment (when available)
│   ├── annotation_evidence.tsv  fragment and homology evidence
│   ├── fragment_metrics.tsv    local-fragment support and state
│   ├── annotated_catalog.tsv    final class, confidence, and decision
│   └── manifest.tsv
└── 04_quantify/
    ├── repeat_signal.tsv        per-sample × EQ RPM, coverage, and annotation
    ├── fragment_coverage.tsv    EQ-level coverage summary from all eligible reads
    └── taxon_repeat_matrix.tsv  taxon median RPM and call
```

`signal_rpm = 1,000,000 × specific_pairs / effective_pairs`. `PRESENT` requires at least 100 effective pairs, 3 specific pairs, and at least 70% specific support. For UCE capture data, interpret results only as off-target repeat signal, not unbiased genome-wide TE content.

## Local assembly and mapping

The default local assembler retains at most three non-redundant fragment hypotheses per EQ. It stops at unsupported sequence and does not join different EQs from PE links. `--te-assemble-min-kmer-count`, `--te-assemble-branch-ratio`, and `--te-assemble-max-fragments` tune this bounded walk. Quantification scans all reads remaining after the optional ledger exclusion. `mean_depth` is mapped bases divided by retained fragment length and `kmer_breadth` is the fraction of unique fragment k-mers observed in specific reads; `anchor_identity` is the mean identity of accepted k-mer-anchored ungapped local alignments, not a gapped alignment identity or genome-wide copy-number estimate.
