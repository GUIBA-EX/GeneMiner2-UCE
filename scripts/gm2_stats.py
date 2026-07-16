from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import csv
import gzip
import os
import statistics


FASTA_EXTENSIONS = ('.fa', '.fas', '.fasta')


def read_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def assembly_row_is_accepted(row):
    if not row:
        return False

    accepted = str(row.get('accepted', '')).strip().lower()
    if accepted:
        return accepted in {'1', 'true', 'yes'}

    low_quality = str(row.get('low_quality', '')).strip().lower()
    return row.get('status') == 'success' and low_quality not in {'1', 'true', 'yes'}


def fmt(value, digits=3):
    if value == '':
        return ''
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f'{value:.{digits}f}'.rstrip('0').rstrip('.')
    return str(value)


def safe_mean(values):
    return statistics.mean(values) if values else ''


def safe_median(values):
    return statistics.median(values) if values else ''


def fasta_lengths(path):
    lengths = []

    with open(path) as handle:
        for _, seq in SimpleFastaParser(handle):
            lengths.append(len(seq.replace('-', '').replace('N', '').replace('n', '')))

    return lengths


def reference_lengths(ref_dir):
    lengths = {}

    with os.scandir(ref_dir) as entries:
        for entry in entries:
            if not entry.is_file() or not entry.name.endswith(FASTA_EXTENSIONS):
                continue

            locus = os.path.splitext(entry.name)[0]
            values = fasta_lengths(entry.path)
            lengths[locus] = round(statistics.mean(values)) if values else 0

    return dict(sorted(lengths.items()))


def read_csv_rows(path):
    if not os.path.isfile(path):
        return []

    with open(path, newline='') as handle:
        return list(csv.DictReader(handle))


def read_filtered_counts(path):
    counts = {}

    if not os.path.isfile(path):
        return counts

    with open(path) as handle:
        for line in handle:
            parts = line.strip().split(',')
            if len(parts) >= 2 and parts[0]:
                counts[parts[0]] = read_int(parts[1])

    return counts


def count_fastq_reads(path):
    opener = gzip.open if path.endswith('.gz') else open

    with opener(path, 'rt') as handle:
        return sum(1 for _ in handle) // 4


def maybe_count_input_reads(read_paths, enabled):
    if not enabled:
        return ''

    total = 0

    for path in dict.fromkeys(read_paths):
        if path and os.path.isfile(path):
            total += count_fastq_reads(path)

    return total


def read_sample_assembly_rows(out_dir, sample):
    rows = {}
    path = os.path.join(out_dir, sample, 'uce_assembly_summary.csv')

    for row in read_csv_rows(path):
        locus = row.get('locus', '')
        if locus:
            rows[locus] = row

    return rows


def collect_rescue_rows(out_dir, samples):
    rows = []
    global_path = os.path.join(out_dir, 'uce_rescue_summary.csv')

    if os.path.isfile(global_path):
        return read_csv_rows(global_path)

    for sample in samples:
        sample_path = os.path.join(out_dir, sample, 'uce_rescue_summary.csv')
        rows.extend(read_csv_rows(sample_path))

    return rows


def build_stats(out_dir, ref_dir, samples, count_input_reads=False):
    refs = reference_lengths(ref_dir)
    loci = list(refs)
    sample_names = list(samples)
    assembly = {
        sample: read_sample_assembly_rows(out_dir, sample)
        for sample in sample_names
    }
    filtered_counts = {
        sample: read_filtered_counts(os.path.join(out_dir, sample, 'ref_reads_count_dict.txt'))
        for sample in sample_names
    }
    rescue_rows = collect_rescue_rows(out_dir, sample_names)

    rescue_by_sample = defaultdict(lambda: defaultdict(int))
    for row in rescue_rows:
        sample = row.get('sample', '')
        status = row.get('rescue_status', '')
        if sample and status:
            rescue_by_sample[sample][status] += 1

    matrices = {
        'lengths': {},
        'read_counts': {},
        'filtered_counts': {},
    }
    sample_stats = []

    for sample in sample_names:
        rows = assembly[sample]
        lengths = {}
        read_counts = {}
        spans = []
        densities = []
        statuses = defaultdict(int)
        total_bases = 0

        for locus in loci:
            row = rows.get(locus, {})
            status = row.get('status', 'missing')
            accepted = assembly_row_is_accepted(row)
            length = read_int(row.get('selected_contig_length')) if accepted else 0
            read_count = read_int(row.get('read_count')) if accepted else 0
            span = read_int(row.get('read_supported_span'))

            if row:
                statuses[status] += 1
            else:
                statuses['missing'] += 1

            lengths[locus] = length
            read_counts[locus] = read_count
            total_bases += length

            if accepted and span:
                spans.append(span)
            if length > 0:
                densities.append(read_count / length)

        matrices['lengths'][sample] = lengths
        matrices['read_counts'][sample] = read_counts
        matrices['filtered_counts'][sample] = {
            locus: filtered_counts[sample].get(locus, 0)
            for locus in loci
        }

        recovered_lengths = [length for length in lengths.values() if length > 0]
        input_reads = maybe_count_input_reads(samples[sample], count_input_reads)
        reads_filtered = sum(filtered_counts[sample].values())
        pct_filtered = '' if input_reads in ('', 0) else reads_filtered / input_reads * 100

        sample_stats.append({
            'Name': sample,
            'InputReads': input_reads,
            'ReadsFiltered': reads_filtered,
            'PctFiltered': pct_filtered,
            'LociWithFilteredReads': sum(1 for value in filtered_counts[sample].values() if value > 0),
            'LociWithContigs': sum(1 for value in recovered_lengths if value > 0),
            'LociSuccess': statuses['success'],
            'LociLowQuality': statuses['low quality'],
            'LociMissing': statuses['missing'],
            'LociNoFilteredFile': statuses['no filtered file'],
            'LociNoSeed': statuses['no seed'],
            'LociNoContigs': statuses['no contigs'],
            'LociInsufficientGenomicKmers': statuses['insufficient genomic kmers'],
            'LociAt25pct': count_loci_at_threshold(lengths, refs, 0.25),
            'LociAt50pct': count_loci_at_threshold(lengths, refs, 0.50),
            'LociAt75pct': count_loci_at_threshold(lengths, refs, 0.75),
            'LociAt150pct': count_loci_at_threshold(lengths, refs, 1.50),
            'RescueSuccess': rescue_by_sample[sample]['success'],
            'RescueFailedRolledBack': rescue_by_sample[sample]['failed_rolled_back'],
            'RescueRevertedDensityDrop': rescue_by_sample[sample]['reverted_density_drop'],
            'RescueRevertedFailed': rescue_by_sample[sample]['reverted_failed_rescue'],
            'TotalBasesRecovered': total_bases,
            'MeanContigLength': safe_mean(recovered_lengths),
            'MedianContigLength': safe_median(recovered_lengths),
            'MeanReadSupportedSpan': safe_mean(spans),
            'MeanReadDensity': safe_mean(densities),
        })

    return refs, loci, sample_names, matrices, sample_stats, build_locus_stats(loci, refs, sample_names, assembly)


def count_loci_at_threshold(lengths, refs, threshold):
    return sum(
        1 for locus, ref_len in refs.items()
        if ref_len > 0 and lengths.get(locus, 0) >= ref_len * threshold
    )


def build_locus_stats(loci, refs, sample_names, assembly):
    rows = []

    for locus in loci:
        lengths = []
        read_counts = []
        spans = []
        balances = []
        candidates = []
        success = 0
        low_quality = 0
        filtered_failures = 0

        for sample in sample_names:
            row = assembly[sample].get(locus, {})
            status = row.get('status', 'missing')
            accepted = assembly_row_is_accepted(row)
            length = read_int(row.get('selected_contig_length')) if accepted else 0

            if status == 'success':
                success += 1
            elif status == 'low quality':
                low_quality += 1
            elif status == 'no filtered file':
                filtered_failures += 1

            if accepted and length > 0:
                lengths.append(length)
                read_counts.append(read_int(row.get('read_count')))
                spans.append(read_int(row.get('read_supported_span')))
                balances.append(read_float(row.get('flank_balance')))
                candidates.append(read_int(row.get('candidate_count')))

        total = len(sample_names)
        rows.append({
            'Locus': locus,
            'MeanReferenceLength': refs.get(locus, 0),
            'Samples': total,
            'SuccessSamples': success,
            'LowQualitySamples': low_quality,
            'NoFilteredFileSamples': filtered_failures,
            'Occupancy': success / total if total else '',
            'MeanLength': safe_mean(lengths),
            'MedianLength': safe_median(lengths),
            'MaxLength': max(lengths) if lengths else '',
            'MeanReadCount': safe_mean(read_counts),
            'MeanReadSupportedSpan': safe_mean(spans),
            'MeanFlankBalance': safe_mean(balances),
            'MeanCandidateCount': safe_mean(candidates),
        })

    return rows


def write_matrix(path, loci, sample_names, matrix, mean_lengths=None):
    with open(path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['Species', *loci])

        if mean_lengths is not None:
            writer.writerow(['MeanLength', *[mean_lengths.get(locus, 0) for locus in loci]])

        for sample in sample_names:
            writer.writerow([sample, *[matrix[sample].get(locus, 0) for locus in loci]])


def write_table(path, rows, fieldnames):
    with open(path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in rows:
            writer.writerow({field: fmt(row.get(field, '')) for field in fieldnames})


def write_rescue_stats(out_dir, rescue_rows):
    if not rescue_rows:
        return

    fieldnames = sorted({field for row in rescue_rows for field in row})
    write_table(os.path.join(out_dir, 'uce_rescue_stats.tsv'), rescue_rows, fieldnames)


def maybe_write_heatmap(matrix_path, out_path, title, relative_to_mean=False):
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
    except ImportError:
        print(f'Warning: pandas, seaborn or matplotlib is unavailable; skipped {out_path}')
        return

    df = pd.read_csv(matrix_path, sep='\t')

    if relative_to_mean:
        df = df.astype('object')
        df.loc[:, df.columns[1]:] = df.loc[:, df.columns[1]:].div(df.iloc[0][df.columns[1]:])
        df.where(df.loc[:, df.columns[1]:] < 1, 1, inplace=True)
        df.drop(labels=0, axis=0, inplace=True)

    df = df.melt(id_vars=['Species'], var_name='locus', value_name='value')
    df['value'] = pd.to_numeric(df['value'])
    df = df.pivot(index='Species', columns='locus', values='value')

    width = max(8, min(400, len(df.columns) / 3))
    height = max(4, min(400, len(df.index) / 3))
    dpi = 100

    plt.figure(figsize=(width, height))
    sns.heatmap(df, vmin=0, cmap='bone_r', xticklabels=True, yticklabels=True)
    plt.title(title)
    plt.xlabel('Locus')
    plt.ylabel('Sample')
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def run(args, samples):
    out_dir = args.o.strip()
    refs, loci, sample_names, matrices, sample_stats, locus_stats = build_stats(
        out_dir,
        args.r,
        samples,
        count_input_reads=args.stats_count_input_reads,
    )

    seq_lengths_path = os.path.join(out_dir, 'uce_seq_lengths.tsv')
    read_counts_path = os.path.join(out_dir, 'uce_read_counts.tsv')

    write_matrix(seq_lengths_path, loci, sample_names, matrices['lengths'], refs)
    write_matrix(read_counts_path, loci, sample_names, matrices['read_counts'])
    write_matrix(os.path.join(out_dir, 'uce_filtered_read_counts.tsv'), loci, sample_names, matrices['filtered_counts'])
    write_table(os.path.join(out_dir, 'uce_stats.tsv'), sample_stats, list(sample_stats[0]) if sample_stats else ['Name'])
    write_table(os.path.join(out_dir, 'uce_locus_stats.tsv'), locus_stats, list(locus_stats[0]) if locus_stats else ['Locus'])
    write_rescue_stats(out_dir, collect_rescue_rows(out_dir, sample_names))

    if not args.stats_no_heatmap:
        maybe_write_heatmap(seq_lengths_path, os.path.join(out_dir, 'uce_recovery_heatmap.png'),
                            'UCE length recovery relative to mean reference length', relative_to_mean=True)
        maybe_write_heatmap(read_counts_path, os.path.join(out_dir, 'uce_read_counts_heatmap.png'),
                            'UCE selected-contig read counts')

    print(f'Wrote UCE statistics to {out_dir}')
