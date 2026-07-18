# Population workflow

[中文版本](population_ZH.md)

The `population` workflow targets multiple diploid UCE samples and builds a cohort pseudo-reference, joint SNPs, PCA, and ADMIXTURE panels from accepted contigs and the original reads. It reports unphased genotypes and does not replace phasing when complete haplotypes or per-locus gene trees are required.

See the [command-line guide](../manual/EN_US/command_line.md) for all options and the [output guide](../manual/EN_US/output.md) for file descriptions.

## Required inputs

Every sample must have completed UCE assembly, retaining:

- `uce_assembly_summary.csv` in the sample directory;
- accepted UCE contigs under `results/`;
- the original R1 and optional R2 listed in the sample table.

Only contigs that pass UCE acceptance criteria and have unique-read support are eligible. All samples are then mapped to one cohort reference instead of calling variants independently against sample-specific assemblies.

## Four stages

```text
reference → mapping → calling → selection
```

1. `reference`: build or copy the cohort pseudo-reference;
2. `mapping`: map with minibwa, then sort, mark duplicates, and collect QC with samtools;
3. `calling`: perform joint calling and SNP filtering with bcftools;
4. `selection`: build three SNP panels, run PLINK PCA, and optionally run ADMIXTURE.

Use `--population-stop-after` to stop after a stage. `--population-start-at mapping|calling|selection` reuses only validated existing outputs. Missing references, BAMs, or filtered VCFs cause an error rather than silent file mixing.

## Cohort pseudo-reference

The default `--population-reference-strategy sqcl-longest` selects the longest accepted contig for each locus, following the representative-sequence idea of SqCL `make_PRG.py`. The `supported` strategy prioritizes read support before length and other metrics.

Internal construction writes:

- `population_reference.fasta`: the common mapping reference;
- `population_reference_provenance.tsv`: source sample, candidate count, and support metrics for each locus;
- `reference_contribution.tsv`: the number and fraction of reference loci contributed by each sample;
- `locus_name_map.tsv`: original locus names and VCF-safe names.

To validate the pseudo-reference or repeat an analysis with a fixed reference, use:

```bash
--population-reference-fasta cohort_reference.fasta
```

The external FASTA is copied into the Population output and its source is recorded in `reference_source.tsv`. The internal pseudo-reference selects one real accepted contig per locus rather than synthesizing a forced multi-individual consensus. Reference contribution and mapping QC must still be checked for bias.

## Mapping and variant calling

The default external tools are minibwa, samtools, and bcftools. Calling defaults are:

```text
minimum MAPQ = 20
minimum base quality = 20
minimum genotype DP = 5
minimum genotype GQ = 20
minimum site QUAL = 20
minimum call rate = 0.8
minimum minor allele count = 2
```

Low-DP or low-GQ genotypes are first set to missing; QUAL, call rate, and MAC are then applied at the site level. Thresholds should follow sequencing depth and missingness patterns rather than being relaxed only to increase the final SNP count.

## SNP panels

The workflow produces three complementary panels:

- `all_snps`: retains multiple SNPs per UCE and maximizes information, but gives linked and long loci more weight;
- `one_snp_per_uce`: selects one representative SNP per UCE and is the default primary panel for PCA and ADMIXTURE;
- `ld_pruned`: prunes the complete SNP set for linkage and serves as a sensitivity panel.

PLINK runs PCA on all three panels. ADMIXTURE uses `one_snp_per_uce` by default and evaluates cross-validation error from `--population-admixture-k-min` through `--population-admixture-k-max`. If ADMIXTURE is unavailable, pseudo-reference construction, VCFs, PLINK files, and PCA still complete, and the status is recorded as `unavailable`.

## Required QC

Before interpreting PCA or ancestry components, inspect:

- `mapping/mapping_qc.tsv`: mapping rate, coverage breadth, and depth per sample;
- `variants/variant_qc.tsv`: variant counts after every filtering stage;
- `reference/reference_contribution.tsv`: whether a few samples dominate an internally built reference;
- `structure/panel_summary.tsv`: SNP counts in all three panels;
- sample and site missingness;
- agreement among the principal structures from all three PCA panels;
- stability of ADMIXTURE CV error across K.

Samples with low mapping or breadth should be checked for reference distance, library quality, sequencing depth, and contamination. Reference contribution is relevant only when the cohort reference was built internally.

## Recommended pseudo-reference validation

For validation, compare:

1. the default `sqcl-longest` strategy;
2. the `supported` strategy;
3. a fixed external reference, when a reliable genome or validated UCE reference is available.

Compare mapping rate, coverage breadth, SNP missingness, PCA axes, and per-locus variant counts across the three analyses. Population structure that appears only under one pseudo-reference should first be treated as reference sensitivity rather than direct biological evidence.
