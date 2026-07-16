from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
import argparse
import csv
import hashlib
import math
import os
import shutil
import statistics
import subprocess
import sys
import threading

import build_trimed
import fix_alignment
import gm2_stats
import muscle_wrapper

COMMAND_HELP = '''
filter    Reference-based filtering of raw reads
refilter  Refinement of filtered reads
assemble  Gene assembly using wDBG
consensus Consensus generation on heterozygous sites
trim      Flank sequence removal
combine   Gene alignment, concatenation and cleanup
tree      Phylogenetic tree reconstruction
stats     UCE recovery statistics and heatmaps
'''

DEPTH_DEPRECATION_EXPLAINER = '''
  Gene assembly involves two different types of depth measurements: k-mer frequencies and read depths. In early versions of GeneMiner2, the two measurements were confounded, causing slight deviation in contig depth calculation. Several new options have since been added to avoid ambiguity.

  --depth-limit INT
    This option limits the highest depth of a gene to prevent pathological assemblies. It should not significantly affect quality unless set to a low value.

  --depth-low-water-mark INT
    This option corresponds to the depth at which GeneMiner2 begins to relax its filtering criteria, so highly diverged sequences can still be recovered. Its value should be sufficiently large, since most irrelevant reads will be removed at a later stage. This option is not for strict quality control!

  --error-threshold INT
    This option exists in previous versions and corresponds to the minimum frequency of a k-mer. Increasing its value will improve quality at the cost of contiguity.

  --min-coverage INT
    This option specifically targets the minimum read depth of contigs. Any contig with a lower read depth will be removed. Note that GeneMiner2 can recover certain types of sequence without a full supporting read. These reads are not counted towards the read depth, making this option very stringent. Use this option if the goal is to ensure read coverage of assembled contigs.
'''

HELP_EPILOG = 'quality control of assembled genes:' + DEPTH_DEPRECATION_EXPLAINER

SCRIPT_ROOT = os.path.join(sys._MEIPASS, os.pardir) if hasattr(sys, '_MEIPASS') else os.path.dirname(__file__)
REFERENCE_EXTENSIONS = ('.fa', '.fas', '.fasta')

def is_reference_file_name(name):
    return os.path.splitext(name)[1].lower() in REFERENCE_EXTENSIONS

def find_executable(prog, internal=False):
    bin_path = os.path.join(SCRIPT_ROOT, prog)

    if not shutil.which(bin_path):
        if internal:
            raise RuntimeError(f"A GeneMiner component is missing from '{bin_path}'")
        else:
            bin_path = shutil.which(prog)

    if not bin_path:
        raise RuntimeError(f"Unable to find {prog} executable")

    return bin_path

def get_ref_genes(ref_dir):
    genes = set()

    for entry in iter_reference_files(ref_dir):
        genes.add(os.path.splitext(entry.name))

    return genes

def get_sample_ext(data_path):
    data_name, data_ext = os.path.splitext(data_path)

    if data_ext == '.gz':
        data_name, data_ext = os.path.splitext(data_name)

    if data_ext == '.fq' or data_ext == '.fastq':
        return '.fq'
    else:
        return '.fasta'

def iter_reference_files(ref_dir):
    with os.scandir(ref_dir) as entries:
        for entry in sorted(entries, key=lambda x: x.name):
            if not entry.is_file():
                continue

            if is_reference_file_name(entry.name):
                yield entry

def reference_cache_key(ref_dir, kmer_size, step_size):
    digest = hashlib.sha256()
    digest.update(os.path.abspath(ref_dir).encode())
    digest.update(b'\0')
    digest.update(str(kmer_size).encode())
    digest.update(b'\0')
    digest.update(str(step_size).encode())

    for entry in iter_reference_files(ref_dir):
        stat = entry.stat()
        digest.update(b'\0')
        digest.update(entry.name.encode())
        digest.update(b'\0')
        digest.update(str(stat.st_size).encode())
        digest.update(b'\0')
        digest.update(str(stat.st_mtime_ns).encode())

    return digest.hexdigest()[:16]

def get_reference_kmer_dict_path(args, out_loc):
    if not args.reuse_reference_cache:
        return os.path.join(out_loc, f'kmer_dict_k{args.kf}.dict')

    cache_dir = args.reference_cache_dir or os.path.join(out_loc, '.gm2_reference_cache')
    cache_name = f'reference_kmer_k{args.kf}_s{args.step_size}_{reference_cache_key(args.r, args.kf, args.step_size)}.dict'
    return os.path.join(cache_dir, cache_name)

def get_assembler_reference_cache_dir(args, out_loc):
    if not args.reuse_reference_cache:
        return None

    return os.path.join(args.reference_cache_dir or os.path.join(out_loc, '.gm2_reference_cache'), 'assembler')

def prepare_workdir(args):
    samples = {}
    tsv_loc = args.f

    try:
        sp_id = 0
        tsv_loc = os.path.realpath(tsv_loc, strict=True)

        with open(tsv_loc, 'r') as f:
            for row in csv.reader(f, delimiter="\t", quotechar='"'):
                if not row:
                    continue

                sp_id += 1
                sp_name = "".join('_' if c in ' -' else c for c in row[0].strip() if c.isalnum() or c in " -_.").capitalize()

                if not sp_name:
                    print(f"Invalid sample name '{row[0]}'")
                    return {}

                if len(row) == 1:
                    print(f"Sample '{row[0]}' has no data files")
                    return {}

                samples[f'{sp_id}_{sp_name}'] = (row[1], row[1] if len(row) == 2 else row[2])

    except OSError as e:
        print(f"Unable to read sample list '{tsv_loc}': {e}")
        return {}

    out_loc = args.o.strip()

    print(f"Preparing working directory '{out_loc}'")

    try:
        os.makedirs(out_loc, exist_ok=True)
    except OSError as e:
        print(f"Unable to create working directory '{out_loc}': {e}")
        return {}

    for name in samples.keys():
        sp_path = os.path.join(out_loc, name)

        try:
            os.makedirs(sp_path, exist_ok=True)
        except OSError as e:
            print(f"Unable to create directory '{sp_path}': {e}")
            return {}

    return samples

def write_fasta_record(out, header, sequence, line_width=80):
    sequence = ''.join(sequence.split()).upper()

    if not sequence:
        return False

    out.write(f'>{header}\n')

    for i in range(0, len(sequence), line_width):
        out.write(sequence[i:i + line_width] + '\n')

    return True

def build_uce_rescue_refs(ref_dir, sample_dir, rescue_ref_dir, min_contig_len):
    results_dir = os.path.join(sample_dir, 'results')
    summary_rows = read_uce_summary(os.path.join(sample_dir, 'uce_assembly_summary.csv'))
    added_contigs = 0

    if not os.path.isdir(results_dir):
        return 0

    if os.path.isdir(rescue_ref_dir):
        shutil.rmtree(rescue_ref_dir, ignore_errors=True)

    os.makedirs(rescue_ref_dir, exist_ok=True)

    with os.scandir(ref_dir) as entries:
        for entry in entries:
            if not entry.is_file():
                continue

            if not is_reference_file_name(entry.name):
                continue

            gene = os.path.splitext(entry.name)[0]
            contig_path = os.path.join(results_dir, f'{gene}.fasta')
            rescue_path = os.path.join(rescue_ref_dir, f'{gene}.fasta')
            contig_index = 0
            wrote_any = False

            with open(rescue_path, 'w') as out:
                with open(entry.path, 'r') as ref_in:
                    for title, seq in SimpleFastaParser(ref_in):
                        wrote_any |= write_fasta_record(out, title, seq)

                if (os.path.isfile(contig_path)
                        and uce_summary_row_is_accepted(summary_rows.get(gene))):
                    with open(contig_path, 'r') as contig_in:
                        for _, seq in SimpleFastaParser(contig_in):
                            if len(seq) < min_contig_len:
                                continue

                            contig_index += 1
                            added_contigs += 1
                            wrote_any |= write_fasta_record(out, f'{gene}_gm2_rescue_contig_{contig_index}', seq)

            if not wrote_any:
                os.remove(rescue_path)

    return added_contigs

def read_uce_summary(summary_path):
    rows = {}

    if not os.path.isfile(summary_path):
        return rows

    with open(summary_path, newline='') as f:
        for row in csv.DictReader(f):
            locus = row.get('locus')

            if locus:
                rows[locus] = row

    return rows


def uce_summary_row_is_accepted(row):
    """Read new acceptance decisions while remaining compatible with old summaries."""
    if not row:
        return False

    accepted = str(row.get('accepted', '')).strip().lower()
    if accepted:
        return accepted in {'1', 'true', 'yes'}

    low_quality = str(row.get('low_quality', '')).strip().lower()
    return row.get('status') == 'success' and low_quality not in {'1', 'true', 'yes'}

def int_or_blank(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return ''


def float_or_blank(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return ''

def delta_or_blank(after, before):
    after_value = int_or_blank(after)
    before_value = int_or_blank(before)

    if after_value == '' or before_value == '':
        return ''

    return after_value - before_value

def read_density_or_blank(row):
    length = int_or_blank(row.get('selected_contig_length'))
    read_count = int_or_blank(row.get('unique_read_count'))
    if length != '' and length > 0 and read_count != '':
        return read_count / length

    unique_density = float_or_blank(row.get('unique_read_density'))
    if unique_density != '':
        return unique_density

    read_count = int_or_blank(row.get('read_count'))
    if length == '' or read_count == '' or length <= 0:
        return ''

    return read_count / length

def density_ratio_or_blank(before, after):
    before_density = read_density_or_blank(before)
    after_density = read_density_or_blank(after)

    if before_density == '' or after_density == '':
        return ''

    if before_density <= 0:
        return ''

    return after_density / before_density

def rescue_density_below_ratio(before, after, min_density_ratio):
    density_ratio = density_ratio_or_blank(before, after)

    if density_ratio == '':
        return False

    return density_ratio < min_density_ratio

UCE_ASSEMBLY_SUMMARY_FIELDS = [
    'locus',
    'status',
    'accepted',
    'rejection_reason',
    'selected_contig_length',
    'read_supported_span',
    'slice_supported_bases',
    'slice_support_breadth',
    'max_slice_support_gap',
    'read_count',
    'unique_read_count',
    'multi_mapping_read_count',
    'read_density',
    'unique_read_density',
    'support_fraction',
    'flank_balance',
    'kmer_median_depth',
    'kmer_depth_cv',
    'kmer_max_depth_ratio',
    'candidate_count',
    'low_quality',
]

UCE_RESCUE_SUMMARY_FIELDS = [
    'sample',
    'locus',
    'rescue_status',
    'before_status',
    'after_status',
    'before_length',
    'after_length',
    'length_delta',
    'before_read_count',
    'after_read_count',
    'read_count_delta',
    'before_read_density',
    'after_read_density',
    'density_ratio',
    'before_read_supported_span',
    'after_read_supported_span',
    'error',
]

SAMPLE_STATE_BACKUP_ITEMS = [
    'results',
    'result_dict.txt',
    'uce_assembly_summary.csv',
    'contigs_all',
    'contigs_all_low',
    'filtered',
    'filtered_pe',
    'ref_reads_count_dict.txt',
]

def write_uce_assembly_summary_rows(summary_path, rows):
    with open(summary_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_ASSEMBLY_SUMMARY_FIELDS)
        writer.writeheader()

        for locus in sorted(rows):
            row = rows[locus]
            writer.writerow({field: row.get(field, '') for field in UCE_ASSEMBLY_SUMMARY_FIELDS})

def write_result_dict_from_uce_summary(sample_dir, rows):
    result_path = os.path.join(sample_dir, 'result_dict.txt')

    with open(result_path, 'w') as out:
        for locus in sorted(rows):
            row = rows[locus]

            if row.get('status') == 'skipped':
                continue

            out.write(f"{locus},{row.get('status', '')},{row.get('read_count', '')},\n")

def restore_locus_file(sample_dir, backup_dir, subdir, locus):
    rel_path = os.path.join(subdir, f'{locus}.fasta')
    src = os.path.join(backup_dir, rel_path)
    dest = os.path.join(sample_dir, rel_path)

    if os.path.isfile(src):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest)
    elif os.path.isfile(dest):
        os.remove(dest)

def locus_file_name_matches(name, locus, paired=False):
    stem = os.path.splitext(name)[0]
    if paired:
        return stem in (f'{locus}_1', f'{locus}_2')
    return stem == locus

def restore_locus_directory_files(sample_dir, backup_dir, subdir, locus):
    """Restore only one locus's read files while keeping accepted rescue loci."""
    source_dir = os.path.join(backup_dir, subdir)
    destination_dir = os.path.join(sample_dir, subdir)
    names = set()

    for directory in (source_dir, destination_dir):
        if os.path.isdir(directory):
            names.update(entry.name for entry in os.scandir(directory) if entry.is_file())

    for name in names:
        if not locus_file_name_matches(name, locus, paired=subdir == 'filtered_pe'):
            continue

        source = os.path.join(source_dir, name)
        destination = os.path.join(destination_dir, name)
        if os.path.isfile(source):
            os.makedirs(destination_dir, exist_ok=True)
            shutil.copy2(source, destination)
        elif os.path.isfile(destination):
            os.remove(destination)

def restore_locus_read_count(sample_dir, backup_dir, locus):
    filename = 'ref_reads_count_dict.txt'
    source = os.path.join(backup_dir, filename)
    destination = os.path.join(sample_dir, filename)

    def read_rows(path):
        if not os.path.isfile(path):
            return []
        with open(path) as handle:
            return [line for line in handle if line.strip()]

    backup_rows = read_rows(source)
    current_rows = read_rows(destination)
    backup_locus_rows = [line for line in backup_rows if line.split(',', 1)[0] == locus]
    merged_rows = [line for line in current_rows if line.split(',', 1)[0] != locus]
    merged_rows.extend(backup_locus_rows)

    if merged_rows:
        with open(destination, 'w') as handle:
            handle.writelines(merged_rows)
    elif os.path.isfile(destination):
        os.remove(destination)

def format_float_or_blank(value, digits=6):
    if value == '':
        return ''

    return f'{value:.{digits}f}'.rstrip('0').rstrip('.')

def revert_invalid_rescue_loci(sample_dir, backup_dir, before_rows, rescue_rows, min_density_ratio):
    reverted = {}
    final_rows = {locus: row.copy() for locus, row in rescue_rows.items()}

    for locus, before in before_rows.items():
        after = rescue_rows.get(locus)
        if not uce_summary_row_is_accepted(before):
            continue

        if not uce_summary_row_is_accepted(after):
            status = 'reverted_failed_rescue'
        elif rescue_density_below_ratio(before, after, min_density_ratio):
            status = 'reverted_density_drop'
        else:
            continue

        for subdir in ('results', 'contigs_all', 'contigs_all_low'):
            restore_locus_file(sample_dir, backup_dir, subdir, locus)
        for subdir in ('filtered', 'filtered_pe'):
            restore_locus_directory_files(sample_dir, backup_dir, subdir, locus)
        restore_locus_read_count(sample_dir, backup_dir, locus)

        final_rows[locus] = before.copy()
        reverted[locus] = status

    if reverted:
        write_uce_assembly_summary_rows(os.path.join(sample_dir, 'uce_assembly_summary.csv'), final_rows)
        write_result_dict_from_uce_summary(sample_dir, final_rows)

    return reverted

def write_sample_uce_rescue_summary(sample_dir, sample, before_rows, after_rows, rescue_status, error='', status_by_locus=None, error_by_locus=None):
    out_path = os.path.join(sample_dir, 'uce_rescue_summary.csv')
    loci = sorted(set(before_rows) | set(after_rows))
    status_by_locus = {} if status_by_locus is None else status_by_locus
    error_by_locus = {} if error_by_locus is None else error_by_locus

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_RESCUE_SUMMARY_FIELDS)
        writer.writeheader()

        for locus in loci:
            before = before_rows.get(locus, {})
            after = after_rows.get(locus, {})
            before_length = before.get('selected_contig_length', '')
            after_length = after.get('selected_contig_length', '')
            before_count = before.get('read_count', '')
            after_count = after.get('read_count', '')
            before_density = read_density_or_blank(before)
            after_density = read_density_or_blank(after)

            writer.writerow({
                'sample': sample,
                'locus': locus,
                'rescue_status': status_by_locus.get(locus, rescue_status),
                'before_status': before.get('status', ''),
                'after_status': after.get('status', ''),
                'before_length': before_length,
                'after_length': after_length,
                'length_delta': delta_or_blank(after_length, before_length),
                'before_read_count': before_count,
                'after_read_count': after_count,
                'read_count_delta': delta_or_blank(after_count, before_count),
                'before_read_density': format_float_or_blank(before_density),
                'after_read_density': format_float_or_blank(after_density),
                'density_ratio': format_float_or_blank(density_ratio_or_blank(before, after)),
                'before_read_supported_span': before.get('read_supported_span', ''),
                'after_read_supported_span': after.get('read_supported_span', ''),
                'error': error_by_locus.get(locus, error),
            })

def backup_sample_state(sample_dir):
    backup_dir = os.path.join(sample_dir, '.uce_rescue_backup')

    if os.path.isdir(backup_dir):
        shutil.rmtree(backup_dir, ignore_errors=True)

    os.makedirs(backup_dir, exist_ok=True)

    for item in SAMPLE_STATE_BACKUP_ITEMS:
        src = os.path.join(sample_dir, item)

        if os.path.exists(src):
            shutil.move(src, os.path.join(backup_dir, item))

    return backup_dir

def restore_sample_state(sample_dir, backup_dir):
    if not os.path.isdir(backup_dir):
        return

    for item in SAMPLE_STATE_BACKUP_ITEMS:
        dest = os.path.join(sample_dir, item)

        if os.path.isdir(dest):
            shutil.rmtree(dest, ignore_errors=True)
        elif os.path.isfile(dest):
            os.remove(dest)

    with os.scandir(backup_dir) as entries:
        for entry in entries:
            dest = os.path.join(sample_dir, entry.name)
            shutil.move(entry.path, dest)

    shutil.rmtree(backup_dir, ignore_errors=True)

def discard_sample_state_backup(backup_dir):
    if os.path.isdir(backup_dir):
        shutil.rmtree(backup_dir, ignore_errors=True)

def write_failed_samples(out_loc, failures):
    out_path = os.path.join(out_loc, 'failed_samples.tsv')

    if not failures:
        if os.path.isfile(out_path):
            os.remove(out_path)

        return

    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample', 'stage', 'error'])
        writer.writerows(failures)

def get_uce_rescue_parallelism(total_threads, sample_count):
    rescue_threads = max(1, min(4, total_threads))
    rescue_workers = max(1, min(4, sample_count, total_threads // rescue_threads))
    return rescue_workers, rescue_threads

def build_uce_rescue_filter_commands(filter_bin, rescue_ref_dir, sample_dir, q1, q2, args, rescue_kmer_dict_path):
    dict_cmd = [filter_bin, '-r', rescue_ref_dir, '-o', sample_dir, '-kf', str(args.kf),
                '-s', str(args.step_size), '-gr', '-lkd', rescue_kmer_dict_path, '-m', '2']
    reads_cmd = [filter_bin, '-r', rescue_ref_dir, '-q1', q1, '-q2', q2, '-o', sample_dir,
                 '-kf', str(args.kf), '-s', str(args.step_size), '-gr', '-subdir', 'filtered_pe',
                 '-m', '5', '-lb', '-lkd', rescue_kmer_dict_path]

    if args.max_reads > 0:
        reads_cmd.extend(['-m_reads', str(args.max_reads)])

    return dict_cmd, reads_cmd

def build_assembler_command(assembler_bin, args, sample_dir, ref_dir, soft_boundary, thr):
    command = [
        assembler_bin, '-r', ref_dir, '-o', sample_dir, '-ka', str(args.ka),
        '-k_min', str(args.min_ka), '-k_max', str(args.max_ka),
        '-limit_count', str(args.error_threshold), '-iteration', str(args.search_depth),
        '-sb', soft_boundary, '-cov_min', str(args.min_coverage), '-p', str(thr),
        '--assembly-mode', args.assembly_mode,
        '--uce-side-candidates', str(args.uce_side_candidates),
        '--uce-max-contig-length', str(args.uce_max_contig_length),
        '--uce-min-read-density', str(args.uce_min_read_density),
        '--uce-density-check-min-length', str(args.uce_density_check_min_length),
        '--uce-max-depth-cv', str(args.uce_max_depth_cv),
        '--uce-max-depth-ratio', str(args.uce_max_depth_ratio),
    ]

    assembler_cache_dir = getattr(args, 'assembler_reference_cache_dir', None)
    original_ref_dir = getattr(args, 'r', None)
    if assembler_cache_dir and original_ref_dir and os.path.abspath(ref_dir) != os.path.abspath(original_ref_dir):
        assembler_cache_dir = None

    if assembler_cache_dir:
        command.extend(['--assembler-reference-cache-dir', assembler_cache_dir])

    return command

def do_filter_assemble(args, samples, do_filter, do_refilter, do_assemble, ignore_hook=lambda *_, **__: None):
    out_loc = args.o.strip()
    kmer_dict_path = get_reference_kmer_dict_path(args, out_loc)
    args.assembler_reference_cache_dir = get_assembler_reference_cache_dir(args, out_loc)
    rescue_enabled = args.uce_rescue_reads
    failed_samples = []
    rescue_workers, rescue_threads = get_uce_rescue_parallelism(args.p, len(samples))

    if rescue_enabled and (args.assembly_mode != 'uce' or not do_filter or not do_refilter or not do_assemble):
        raise RuntimeError('--uce-rescue-reads requires --assembly-mode uce and the filter, refilter and assemble steps')

    if args.soft_boundary == 'auto':
        soft_boundary = -1
    elif args.soft_boundary == 'unlimited':
        soft_boundary = 10000
    else:
        try:
            soft_boundary = int(args.soft_boundary)
        except ValueError:
            raise RuntimeError(f"Invalid soft boundary {args.soft_boundary} (must be an integer)")

        if soft_boundary < 0:
            raise RuntimeError(f"Invalid soft boundary {args.soft_boundary} (must be positive or zero)")

    soft_boundary = str(soft_boundary)

    if do_filter:
        filter_bin = find_executable('MainFilterNew', internal=True)

        os.makedirs(os.path.dirname(kmer_dict_path), exist_ok=True)

        if os.path.isfile(kmer_dict_path) and args.reuse_reference_cache:
            print(f'Reusing reference k-mer cache: {kmer_dict_path}')
        elif os.path.isfile(kmer_dict_path):
            os.remove(kmer_dict_path)

        if not os.path.isfile(kmer_dict_path):
            try:
                subprocess.run([filter_bin, '-r', args.r, '-o', out_loc, '-kf', str(args.kf), '-s', str(args.step_size),
                                '-gr', '-lkd', kmer_dict_path, '-m', '2'], check=True)
            except subprocess.SubprocessError as e:
                raise RuntimeError(f"Unable to build k-mer dictionary: {e}")

        def run_filter(name):
            q1, q2 = samples[name]
            read_count_path = os.path.join(out_loc, name, 'ref_reads_count_dict.txt')
            out_dir = os.path.join(out_loc, name, 'filtered_pe')

            if os.path.isfile(read_count_path):
                os.remove(read_count_path)

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            params = [filter_bin, '-r', args.r, '-q1', q1, '-q2', q2, '-o', os.path.join(out_loc, name),
                      '-kf', str(args.kf), '-s', str(args.step_size), '-gr', '-subdir', 'filtered_pe',
                      '-m', '5', '-lb', '-lkd', kmer_dict_path]

            if args.max_reads > 0:
                params.extend(['-m_reads', str(args.max_reads)])

            subprocess.run(params, check=True)

            if not os.path.isfile(read_count_path):
                raise RuntimeError('Filter failed')

            if not do_refilter and os.path.isdir(out_dir):
                merge_dir = os.path.join(out_loc, name, 'filtered')
                sample_ext = get_sample_ext(q1)

                if os.path.isdir(merge_dir):
                    shutil.rmtree(merge_dir, ignore_errors=True)

                os.makedirs(merge_dir, exist_ok=True)

                genes = set()

                with open(read_count_path, 'r') as f:
                    for line in f:
                        line = line.strip()

                        if not line:
                            continue

                        genes.add(line.split(',')[0])

                for gene in genes:
                    read_1 = os.path.join(out_dir, f'{gene}_1{sample_ext}')
                    read_2 = os.path.join(out_dir, f'{gene}_2{sample_ext}')

                    if not os.path.isfile(read_1):
                        continue

                    with open(os.path.join(merge_dir, gene + sample_ext), 'wb') as f:
                        with open(read_1, 'rb') as r:
                            shutil.copyfileobj(r, f)

                        if not os.path.isfile(read_2):
                            continue

                        with open(read_2, 'rb') as r:
                            shutil.copyfileobj(r, f)

    else:
        run_filter = ignore_hook

    if do_refilter:
        refilter_bin = find_executable('main_refilter_new', internal=True)

        def run_refilter(name, thr=1, ref_dir=None):
            in_dir  = os.path.join(out_loc, name, 'filtered_pe')
            out_dir = os.path.join(out_loc, name, 'filtered')
            ref_dir = args.r if ref_dir is None else ref_dir

            if not os.path.isdir(in_dir):
                raise RuntimeError('No successful filter run, cannot re-filter')

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            params = [refilter_bin, '-r', ref_dir, '-qd', in_dir, '-o', out_dir, '-kf', str(args.kf),
                      '-p', str(thr), '--log-file', os.path.join(out_loc, name, 'log.txt'),
                      '--min-depth', str(args.depth_low_water_mark), '--max-depth', str(args.depth_limit),
                      '--max-size', str(args.file_size_limit), '--use-gm2-format']

            if args.assembly_mode == 'uce':
                params.append('--keep-linked-mates')

            subprocess.run(params, check=True)

            if do_filter and os.path.isdir(in_dir) and os.path.isdir(out_dir):
                shutil.rmtree(in_dir, ignore_errors=True)

    else:
        run_refilter = ignore_hook

    if do_assemble:
        assembler_bin = find_executable('main_assembler', internal=True)

        def run_assembler(name, thr=1, ref_dir=None):
            sample_dir = os.path.join(out_loc, name)
            in_dir = os.path.join(sample_dir, 'filtered')
            out_dir = os.path.join(sample_dir, 'results')
            result_path = os.path.join(sample_dir, 'result_dict.txt')
            ref_dir = args.r if ref_dir is None else ref_dir

            if not os.path.isdir(in_dir):
                raise RuntimeError('No successful filter run, cannot assemble')

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir, ignore_errors=True)

            if os.path.isfile(result_path):
                os.remove(result_path)

            uce_summary_path = os.path.join(sample_dir, 'uce_assembly_summary.csv')

            if os.path.isfile(uce_summary_path):
                os.remove(uce_summary_path)

            subprocess.run(build_assembler_command(assembler_bin, args, sample_dir, ref_dir, soft_boundary, thr), check=True)

            if not os.path.isfile(result_path):
                raise RuntimeError('Assembly failed')

    else:
        run_assembler = ignore_hook

    if rescue_enabled:
        def run_uce_rescue(name, thr=1):
            sample_dir = os.path.join(out_loc, name)
            rescue_ref_dir = os.path.join(sample_dir, 'uce_rescue_refs')
            rescue_kmer_dict_path = os.path.join(sample_dir, f'uce_rescue_kmer_dict_k{args.kf}.dict')
            summary_path = os.path.join(sample_dir, 'uce_assembly_summary.csv')
            read_count_path = os.path.join(sample_dir, 'ref_reads_count_dict.txt')
            filtered_pe_dir = os.path.join(sample_dir, 'filtered_pe')
            filtered_dir = os.path.join(sample_dir, 'filtered')

            before_rows = read_uce_summary(summary_path)
            added_contigs = build_uce_rescue_refs(args.r, sample_dir, rescue_ref_dir, args.uce_rescue_min_contig_length)

            if added_contigs == 0:
                print(f'No preliminary UCE contigs for {name}; skipping raw-read rescue.')
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, before_rows, 'skipped')
                return

            print(f'Running one-round UCE raw-read rescue for {name} using {added_contigs} preliminary contig(s).')
            backup_dir = backup_sample_state(sample_dir)

            try:
                if os.path.isfile(rescue_kmer_dict_path):
                    os.remove(rescue_kmer_dict_path)

                q1, q2 = samples[name]
                dict_cmd, reads_cmd = build_uce_rescue_filter_commands(
                    filter_bin,
                    rescue_ref_dir,
                    sample_dir,
                    q1,
                    q2,
                    args,
                    rescue_kmer_dict_path,
                )

                subprocess.run(dict_cmd, check=True)

                if os.path.isfile(read_count_path):
                    os.remove(read_count_path)

                if os.path.isdir(filtered_pe_dir):
                    shutil.rmtree(filtered_pe_dir, ignore_errors=True)

                if os.path.isdir(filtered_dir):
                    shutil.rmtree(filtered_dir, ignore_errors=True)

                subprocess.run(reads_cmd, check=True)

                if not os.path.isfile(read_count_path):
                    raise RuntimeError('UCE rescue filter failed')

                run_refilter(name, thr=thr, ref_dir=rescue_ref_dir)
                run_assembler(name, thr=thr, ref_dir=rescue_ref_dir)

            except Exception as e:
                restore_sample_state(sample_dir, backup_dir)
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, before_rows, 'failed_rolled_back', str(e))
                print(f'Warning: UCE raw-read rescue failed for {name}; first-round assembly was restored: {e}')
                return

            else:
                after_rows = read_uce_summary(summary_path)
                reverted_loci = revert_invalid_rescue_loci(
                    sample_dir,
                    backup_dir,
                    before_rows,
                    after_rows,
                    args.uce_rescue_min_density_ratio,
                )
                status_by_locus = reverted_loci
                error_by_locus = {}
                for locus, status in reverted_loci.items():
                    if status == 'reverted_density_drop':
                        error_by_locus[locus] = (
                            f'rescue unique read density ratio below '
                            f'{args.uce_rescue_min_density_ratio:g}; first-round contig restored')
                    else:
                        error_by_locus[locus] = (
                            'rescue result missing or rejected; first-round contig restored')
                final_rows = read_uce_summary(summary_path)
                write_sample_uce_rescue_summary(sample_dir, name, before_rows, final_rows, 'success', status_by_locus=status_by_locus, error_by_locus=error_by_locus)
                discard_sample_state_backup(backup_dir)
    else:
        run_uce_rescue = ignore_hook

    if args.p > 1:
        avail_cpu = args.p
        asm_thr   = max(min(args.p // 2, 6), 2)
        filt_thr  = 1 if args.p < 4 else 2

        def calc_task_thr():
            min_thr = min(asm_thr, filt_thr) if filter_list else asm_thr
            return avail_cpu if avail_cpu - asm_thr < min_thr else asm_thr

        filter_list   = []
        refilter_list = []
        assemble_list = []

        executor      = ThreadPoolExecutor(max_workers=math.ceil(avail_cpu / filt_thr))
        running_tasks = {}
        task_metadata = {} # (stage, threads)

        if do_filter:
            filter_list.extend(reversed(samples.keys()))
        elif do_refilter:
            refilter_list.extend(reversed(samples.keys()))
        elif do_assemble:
            assemble_list.extend(reversed(samples.keys()))

        while True:
            while refilter_list and avail_cpu >= filt_thr:
                sample = refilter_list.pop()
                task_thr = calc_task_thr()
                avail_cpu -= task_thr
                running_tasks[sample] = executor.submit(run_refilter, sample, thr=task_thr)
                task_metadata[sample] = (2, task_thr)

            while assemble_list and avail_cpu >= asm_thr:
                sample = assemble_list.pop()
                task_thr = calc_task_thr()
                avail_cpu -= task_thr
                running_tasks[sample] = executor.submit(run_assembler, sample, thr=task_thr)
                task_metadata[sample] = (3, task_thr)

            while filter_list and avail_cpu >= filt_thr:
                sample = filter_list.pop()
                avail_cpu -= filt_thr
                running_tasks[sample] = executor.submit(run_filter, sample)
                task_metadata[sample] = (1, filt_thr)

            if not running_tasks:
                break

            wait(running_tasks.values(), return_when=FIRST_COMPLETED)

            processed_samples = set()

            for sample, task in running_tasks.items():
                if not task.done():
                    continue

                processed_samples.add(sample)

                stage, thr_cnt = task_metadata[sample]
                avail_cpu += thr_cnt

                try:
                    task.result()
                except Exception as e:
                    print(f'An error occurred while processing {sample}: {e}')
                    failed_samples.append((sample, {1: 'filter', 2: 'refilter', 3: 'assemble'}.get(stage, 'unknown'), str(e)))
                    continue

                if stage == 1:
                    if do_refilter:
                        refilter_list.append(sample)
                    elif do_assemble:
                        assemble_list.append(sample)
                elif stage == 2 and do_assemble:
                    assemble_list.append(sample)

            for sample in processed_samples:
                del running_tasks[sample]
                del task_metadata[sample]

        executor.shutdown()

    else:
        for name in samples.keys():
            try:
                if do_filter:
                    run_filter(name)
                if do_refilter:
                    run_refilter(name)
                if do_assemble:
                    run_assembler(name)
            except Exception as e:
                print(f'An error occurred while processing {name}: {e}')
                failed_samples.append((name, 'filter/refilter/assemble', str(e)))
                continue

    if rescue_enabled:
        print(f'Running UCE raw-read rescue with up to {rescue_workers} sample(s) in parallel and {rescue_threads} thread(s) per sample.')

        if rescue_workers > 1:
            with ThreadPoolExecutor(max_workers=rescue_workers) as executor:
                running_rescues = {
                    executor.submit(run_uce_rescue, name, thr=rescue_threads): name
                    for name in samples.keys()
                }

                for task in as_completed(running_rescues):
                    name = running_rescues[task]

                    try:
                        task.result()
                    except Exception as e:
                        print(f'An error occurred during UCE raw-read rescue for {name}: {e}')
                        failed_samples.append((name, 'uce_rescue', str(e)))
        else:
            for name in samples.keys():
                try:
                    run_uce_rescue(name, thr=rescue_threads)
                except Exception as e:
                    print(f'An error occurred during UCE raw-read rescue for {name}: {e}')
                    failed_samples.append((name, 'uce_rescue', str(e)))
                    continue

    write_failed_samples(out_loc, failed_samples)

    if failed_samples:
        raise RuntimeError(f'{len(failed_samples)} sample task(s) failed; see {os.path.join(out_loc, "failed_samples.tsv")}')

def make_phyluce_sample_name(sample):
    name = ''.join(c if ord(c) < 128 and (c.isalnum() or c == '_') else '_' for c in sample).strip('_')

    if not name:
        name = 'sample'

    if not name[0].isalpha():
        name = 'sample_' + name

    return name

def get_contig_read_count(header):
    parts = header.split('_')

    if len(parts) >= 6 and parts[0] == 'contig' and parts[5].isdigit():
        return parts[5]

    if parts and parts[-1].isdigit():
        return parts[-1]

    return '0'

def write_uce_contigs_for_phyluce(args, samples):
    out_loc = args.o.strip()
    uce_dir = os.path.join(out_loc, 'uce_contigs')

    if os.path.isdir(uce_dir):
        shutil.rmtree(uce_dir, ignore_errors=True)

    os.makedirs(uce_dir, exist_ok=True)

    name_map_path = os.path.join(uce_dir, 'sample_name_map.tsv')

    with open(name_map_path, 'w', newline='') as map_file:
        writer = csv.writer(map_file, delimiter='\t')
        writer.writerow(['geneminer_sample', 'phyluce_sample', 'contigs_file', 'contig_count'])
        used_names = set()

        for sample in samples.keys():
            phyluce_sample = make_phyluce_sample_name(sample)

            if phyluce_sample in used_names:
                suffix = 2
                base_name = phyluce_sample

                while f'{base_name}_{suffix}' in used_names:
                    suffix += 1

                phyluce_sample = f'{base_name}_{suffix}'

            used_names.add(phyluce_sample)
            results_dir = os.path.join(out_loc, sample, 'results')
            summary_rows = read_uce_summary(os.path.join(out_loc, sample, 'uce_assembly_summary.csv'))
            accepted_loci = {locus for locus, row in summary_rows.items()
                             if uce_summary_row_is_accepted(row)}
            out_path = os.path.join(uce_dir, phyluce_sample + '.contigs.fasta')
            contig_count = 0

            with open(out_path, 'w') as out:
                if os.path.isdir(results_dir):
                    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
                        if not entry.is_file() or not is_reference_file_name(entry.name):
                            continue

                        locus = os.path.splitext(entry.name)[0]
                        if locus not in accepted_loci:
                            continue

                        with open(entry.path) as f:
                            for header, seq in SimpleFastaParser(f):
                                contig_count += 1
                                read_count = get_contig_read_count(header)
                                out.write(f'>NODE_{contig_count}_length_{len(seq)}_cov_{read_count}.0_{locus}\n')
                                out.write(seq + '\n')

            writer.writerow([sample, phyluce_sample, out_path, contig_count])

def write_uce_assembly_summary(args, samples):
    out_loc = args.o.strip()
    out_path = os.path.join(out_loc, 'uce_assembly_summary.csv')
    fieldnames = ['sample', *UCE_ASSEMBLY_SUMMARY_FIELDS]

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples.keys():
            summary_path = os.path.join(out_loc, sample, 'uce_assembly_summary.csv')

            if not os.path.isfile(summary_path):
                continue

            with open(summary_path, newline='') as f:
                for row in csv.DictReader(f):
                    row['sample'] = sample
                    writer.writerow({name: row.get(name, '') for name in fieldnames})

def write_uce_rescue_summary(args, samples):
    out_loc = args.o.strip()
    out_path = os.path.join(out_loc, 'uce_rescue_summary.csv')
    wrote_any = False

    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=UCE_RESCUE_SUMMARY_FIELDS)
        writer.writeheader()

        for sample in samples.keys():
            summary_path = os.path.join(out_loc, sample, 'uce_rescue_summary.csv')

            if not os.path.isfile(summary_path):
                continue

            with open(summary_path, newline='') as f:
                for row in csv.DictReader(f):
                    writer.writerow({name: row.get(name, '') for name in UCE_RESCUE_SUMMARY_FIELDS})
                    wrote_any = True

    if not wrote_any:
        os.remove(out_path)

def write_uce_outputs(args, samples):
    write_uce_contigs_for_phyluce(args, samples)
    write_uce_assembly_summary(args, samples)

    if args.uce_rescue_reads:
        write_uce_rescue_summary(args, samples)

def generate_consensus(args, samples):
    out_loc = args.o.strip()

    consensus_bin = find_executable('build_consensus', internal=True)
    minimap2_bin = find_executable('minimap2')

    if args.consensus_threshold <= 0 or args.consensus_threshold > 1:
        raise RuntimeError(f"Invalid consensus threshold {args.consensus_threshold} (must be between 0.0 and 1.0)")

    genes = get_ref_genes(args.r)

    def iterate_gene(sample):
        in_dir = os.path.join(out_loc, sample, 'results')

        if not os.path.isdir(in_dir):
            print(f'Error: Sample {sample} has no assembled genes, cannot generate consensus')
            return

        cns_dir  = os.path.join(out_loc, sample, 'consensus')
        filt_dir = os.path.join(out_loc, sample, 'filtered')

        if os.path.isdir(cns_dir):
            shutil.rmtree(cns_dir, ignore_errors=True)

        os.makedirs(cns_dir, exist_ok=True)

        for name, ext in genes:
            asm_path = os.path.join(in_dir, name + ext)
            read_path = os.path.join(filt_dir, name + get_sample_ext(samples[sample][0]))

            if os.path.isfile(asm_path) and os.path.isfile(read_path):
                sam_path = os.path.join(cns_dir, name + '.sam')
                yield (name, asm_path, read_path, sam_path)

    def process_gene(task):
        gene, asm_path, read_path, sam_path = task

        subprocess.run([minimap2_bin, '-ax', 'sr', '-t', '1', '--sam-hit-only', '--secondary=no',
                        '-o', sam_path, asm_path, read_path],
                       check=True)

        if os.path.isfile(sam_path):
            subprocess.run([consensus_bin, '-i', sam_path, '-c', str(args.consensus_threshold),
                            '-o', os.path.dirname(sam_path), '-s', '0'],
                           check=True)

            os.remove(sam_path)

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for _ in executor.map(process_gene, (task
                                                for sample in samples.keys()
                                                for task in iterate_gene(sample))):
                pass

    else:
        for sample in samples.keys():
            for task in iterate_gene(sample):
                process_gene(task)

def blast_trim(args, samples):
    out_loc = args.o.strip()

    makeblastdb_bin = find_executable('makeblastdb')

    if args.trim_mode == 'isoform':
        blast_bin = find_executable('magicblast')
        blast_iter = build_trimed.execute_magicblast
    else:
        blast_bin = find_executable('blastn')
        blast_iter = build_trimed.execute_blastn

    if args.trim_retention < 0 or args.trim_retention > 1:
        raise RuntimeError(f"Invalid trim retention threshold {args.trim_retention} (must be between 0.0 and 1.0)")

    if args.trim_mode == 'longest' or args.trim_mode == 'isoform':
        criterion = 'longest'
    elif args.trim_mode == 'terminal':
        criterion = 'terminal'
    else:
        criterion = 'all'

    genes = get_ref_genes(args.r)

    os.makedirs(os.path.join(out_loc, 'blast_db'), exist_ok=True)

    def build_blast_db(name_tup):
        name, ext = name_tup
        subprocess.run([makeblastdb_bin, "-in", os.path.realpath(os.path.join(args.r, name + ext)),
                        "-dbtype", "nucl", "-out", name],
                       cwd=os.path.join(out_loc, 'blast_db'), check=True)

    def iterate_gene(sample):
        if args.trim_source == 'consensus':
            in_dir = os.path.join(out_loc, sample, 'consensus')
        else:
            in_dir = os.path.join(out_loc, sample, 'results')

        if not os.path.isdir(in_dir):
            print(f'Error: Sample {sample} has no {args.trim_source} sequences, cannot trim')
            return

        blast_dir = os.path.join(out_loc, sample, 'blast')

        if os.path.isdir(blast_dir):
            shutil.rmtree(blast_dir, ignore_errors=True)

        os.makedirs(blast_dir, exist_ok=True)

        for name, ext in genes:
            asm_path = os.path.join(in_dir, name + '.fasta')
            ref_path = os.path.join(args.r, name + ext)

            if os.path.isfile(asm_path):
                yield (name, asm_path, ref_path, os.path.join(blast_dir, name + '.fasta'))

    def process_gene(task):
        name, asm_path, ref_path, out_path = task
        blast_output = blast_iter(asm_path, os.path.join(out_loc, 'blast_db', name), executable_path=blast_bin)
        build_trimed.process_file(asm_path, ref_path, blast_output, out_path, args.trim_retention * 100, criterion)

    gene_count = len(genes) * len(samples)
    trimmed_count = 0

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for _ in executor.map(build_blast_db, genes):
                pass

            for _ in executor.map(process_gene, (task
                                                for sample in samples.keys()
                                                for task in iterate_gene(sample))):
                trimmed_count += 1

                if trimmed_count >= 2:
                    print(f'{trimmed_count}/{gene_count} genes trimmed\r', end='')

    else:
        for gene_tup in genes:
            build_blast_db(gene_tup)

        for sample in samples.keys():
            for task in iterate_gene(sample):
                process_gene(task)

                trimmed_count += 1

                if trimmed_count >= 2:
                    print(f'{trimmed_count}/{gene_count} genes trimmed\r', end='')

    print('\n')

def combine_genes(args, samples):
    out_loc = args.o.strip()
    alignment_filter = get_alignment_filter(args)

    if not args.no_alignment:
        msa_threads = get_msa_threads(args)
        filter_processes = get_filter_processes(args)

        alifilter_model = get_alifilter_model(args)

        if alifilter_model and alignment_filter != 'alifilter':
            raise RuntimeError("--alifilter-model requires --alignment-filter alifilter")

        if args.msa_program == 'clustalo':
            msa_bin = find_executable('clustalo')
        elif args.msa_program == 'muscle':
            msa_bin = find_executable('muscle')
        else:
            msa_bin = find_executable('mafft')

        if alignment_filter == 'trimal':
            trimal_bin = find_executable('trimal')
            alifilter_bin = None
        elif alignment_filter == 'alifilter':
            trimal_bin = None
            alifilter_bin = find_executable('AliFilter')
        else:
            trimal_bin = None
            alifilter_bin = None

        merge_seq_bin = find_executable('merge_seq', internal=True)

        if args.clean_difference < 0 or args.clean_difference > 1:
            raise RuntimeError(f"Invalid maximum difference {args.clean_difference} (must be between 0.0 and 1.0)")

        if args.clean_sequences < 0 or args.clean_sequences > len(samples):
            raise RuntimeError(f"Invalid required number of sequences {args.clean_sequences} (must be between 0 and {len(samples)})")

    combine_dir = os.path.join(out_loc, 'combined_results')

    if os.path.isdir(combine_dir):
        shutil.rmtree(combine_dir, ignore_errors=True)

    os.makedirs(combine_dir, exist_ok=True)

    if not args.no_alignment:
        alignment_dir = os.path.join(out_loc, 'combined_results', 'aligned')
        trim_dir = os.path.join(out_loc, 'combined_trimed')
        combined_fasta = os.path.join(out_loc, 'combined_results.fasta')
        trim_fasta = os.path.join(out_loc, 'combined_trimed.fasta')

        if os.path.isdir(trim_dir):
            shutil.rmtree(trim_dir, ignore_errors=True)

        if os.path.isfile(trim_fasta):
            os.remove(trim_fasta)

        if os.path.isfile(combined_fasta):
            os.remove(combined_fasta)

        os.makedirs(alignment_dir, exist_ok=True)

        if alignment_filter != 'none':
            os.makedirs(trim_dir, exist_ok=True)

    if args.combine_source == 'trimmed':
        in_name = 'blast'
    elif args.combine_source == 'consensus':
        in_name = 'consensus'
    else:
        in_name = 'results'

    genes = {t[0] for t in get_ref_genes(args.r)}
    accepted_loci_by_sample = {}
    if getattr(args, 'assembly_mode', 'reference') == 'uce':
        for name in samples.keys():
            summary_rows = read_uce_summary(
                os.path.join(out_loc, name, 'uce_assembly_summary.csv'))
            accepted_loci_by_sample[name] = {
                locus for locus, row in summary_rows.items()
                if uce_summary_row_is_accepted(row)
            }

    def merge_gene(gene):
        out_path = os.path.join(combine_dir, gene + '.fasta')
        written  = False

        with open(out_path, 'w+') as f:
            for name in samples.keys():
                if (accepted_loci_by_sample
                        and gene not in accepted_loci_by_sample.get(name, set())):
                    continue

                in_path = os.path.join(out_loc, name, in_name, gene + '.fasta')

                if not os.path.isfile(in_path):
                    continue

                with open(in_path, 'r') as r:
                    try:
                        _, seq = next(SimpleFastaParser(r))
                    except StopIteration:
                        continue

                    f.write(f'>{name}\n{seq}\n')

                written = True

        if not written:
            os.remove(out_path)

    def handle_locus_failure(gene, stage, exc=None, output_path=None):
        if output_path and os.path.isfile(output_path):
            os.remove(output_path)

        message = f"{stage} failed on {gene}"

        if args.strict_combine_errors:
            raise RuntimeError(message) from exc

        if exc:
            message = f"{message}: {exc}"

        print(f"Warning: {message}", file=sys.stderr)
        return False

    def align_gene(gene):
        in_path = os.path.join(combine_dir, gene + '.fasta')
        out_path = os.path.join(alignment_dir, gene + '.fasta')

        if not os.path.isfile(in_path):
            return False

        try:
            if args.msa_program == 'clustalo':
                subprocess.run([msa_bin, '-i', in_path, '-o', out_path, '--auto', '--force',
                                '--seqtype=DNA', f'--threads={msa_threads}'],
                               stderr=subprocess.DEVNULL, check=True)
            elif args.msa_program == 'muscle':
                subprocess.run([msa_bin, '-align', in_path, '-output', out_path, '-quiet',
                                '-nt', '-threads', str(msa_threads)],
                               stderr=subprocess.DEVNULL, check=True)
                muscle_wrapper.reorder_sequences(in_path, out_path)
            else:
                with open(out_path, 'w') as out:
                    subprocess.run([msa_bin, '--auto', '--quiet', '--nuc',
                                    '--thread', str(msa_threads), in_path],
                                   stdout=out, stderr=subprocess.DEVNULL, check=True)
        except (OSError, subprocess.CalledProcessError, RuntimeError, ValueError) as e:
            return handle_locus_failure(gene, args.msa_program, e, out_path)

        return os.path.isfile(out_path)

    def clean_gene(gene):
        gene_path = os.path.join(alignment_dir, gene + '.fasta')

        if os.path.isfile(gene_path):
            try:
                fix_alignment.clean_file(gene_path, args.clean_sequences, args.clean_difference)
            except (OSError, RuntimeError, ValueError) as e:
                return handle_locus_failure(gene, 'alignment cleanup', e, gene_path)

        return os.path.isfile(gene_path)

    def filter_gene(gene):
        in_path = os.path.join(alignment_dir, gene + '.fasta')
        out_path = os.path.join(trim_dir, gene + '.fasta')

        if not os.path.isfile(in_path):
            return False

        if alignment_filter == 'trimal':
            cmd = [trimal_bin, '-in', in_path, '-out', out_path, '-automated1']
        elif alignment_filter == 'alifilter':
            cmd = [alifilter_bin, '-i', in_path, '-o', out_path]

            if alifilter_model:
                cmd.extend(['-m', alifilter_model])
        else:
            return False

        try:
            subprocess.run(cmd, check=True)
        except (OSError, subprocess.CalledProcessError) as e:
            return handle_locus_failure(gene, alignment_filter, e, out_path)

        return os.path.isfile(out_path)

    alignment_count = 0
    gene_count = len(genes)

    if args.no_alignment:
        process_gene = merge_gene
    else:
        msa_slots = max(1, args.p // msa_threads)
        msa_semaphore = threading.Semaphore(msa_slots)
        filter_semaphore = threading.Semaphore(filter_processes)

        def process_gene(gene):
            merge_gene(gene)

            with msa_semaphore:
                aligned = align_gene(gene)

            if not aligned:
                return False

            cleaned = clean_gene(gene)

            if not cleaned:
                return False

            if alignment_filter != 'none':
                with filter_semaphore:
                    filtered = filter_gene(gene)

                if not filtered:
                    return False

            return True

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            for aligned in executor.map(process_gene, genes):
                if aligned:
                    alignment_count += 1

                    if alignment_count >= 2:
                        print(f'{alignment_count}/{gene_count} genes aligned\r', end='')

    else:
        for gene in genes:
            aligned = process_gene(gene)
            if aligned:
                alignment_count += 1

                if alignment_count >= 2:
                    print(f'{alignment_count}/{gene_count} genes aligned\r', end='')

    print('\n')

    if not args.no_alignment:
        subprocess.run([merge_seq_bin, '-input', alignment_dir, '-exts', '.fasta', '-missing', '-',
                        '-output', os.path.join(out_loc, 'combined_results.fasta')],
                       check=True)

        if alignment_filter != 'none':
            subprocess.run([merge_seq_bin, '-input', trim_dir, '-exts', '.fasta', '-missing', '-',
                            '-output', os.path.join(out_loc, 'combined_trimed.fasta')],
                           check=True)

def get_alignment_filter(args):
    if getattr(args, 'no_trimal', False):
        return 'none'

    return getattr(args, 'alignment_filter', None) or 'trimal'

def get_alifilter_model(args):
    model = getattr(args, 'alifilter_model', None)

    if not model:
        return None

    model = model.strip()

    if not model or model.lower() == 'default':
        return None

    return model

def get_msa_threads(args):
    return max(1, getattr(args, 'msa_threads', 1))

def get_filter_processes(args):
    filter_processes = getattr(args, 'filter_processes', None)

    if filter_processes is None:
        return max(1, args.p)

    return max(1, filter_processes)

def get_locus_alignment_dir(args):
    out_loc = args.o.strip()

    if get_alignment_filter(args) == 'none':
        return os.path.join(out_loc, 'combined_results', 'aligned')

    return os.path.join(out_loc, 'combined_trimed')

def get_concatenated_alignment_path(args):
    out_loc = args.o.strip()

    if get_alignment_filter(args) == 'none':
        return os.path.join(out_loc, 'combined_results.fasta')

    return os.path.join(out_loc, 'combined_trimed.fasta')

def build_single_tree(prog_name, prog_bin, in_path, bootstrap=0, quiet=False, threads=1):
    if prog_name == 'raxmlng':
        out_path = in_path + ".raxml.bestTree"
        params = [prog_bin, '--msa', in_path, '--msa-format', 'FASTA',
                  '--model', 'GTR+G', '--redo']

        if bootstrap:
            params.extend(['--all', '--bs-trees', str(bootstrap)])
        else:
            params.append('--search')

        if threads > 1:
            params.extend(['--threads', f'auto{{{threads}}}', '--workers', 'auto'])
        else:
            params.extend(['--threads', '1'])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stdout=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    elif prog_name == 'iqtree':
        out_path = in_path + ".treefile"
        params = [prog_bin, '-s', in_path, '-redo']

        if bootstrap:
            params.extend(['-B', str(bootstrap)])

        if quiet:
            params.append('-quiet')

        if threads > 1:
            params.extend(['-T', 'AUTO', '-ntmax', str(threads)])
        else:
            params.extend(['-T', '1'])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stdout=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    elif prog_name == 'veryfasttree':
        out_path = in_path + ".veryfasttree.tre"
        params = [prog_bin, '-out', out_path, '-gtr']

        if bootstrap:
            params.extend(['-boot', str(bootstrap)])
        else:
            params.append('-nosupport')

        if quiet:
            params.extend(['-quiet'])

        if threads > 1:
            params.extend(['-threads', str(threads)])

        params.extend(['-nt', in_path])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stderr=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

    else:
        out_path = in_path + ".fasttree.tre"
        params = [prog_bin, '-out', out_path, '-gtr']

        if bootstrap:
            params.extend(['-boot', str(bootstrap)])
        else:
            params.append('-nosupport')

        if quiet:
            params.append('-quiet')

        params.extend(['-nt', in_path])

        if os.path.isfile(out_path):
            os.remove(out_path)

        subprocess.run(params, stderr=subprocess.DEVNULL if quiet else None, check=True)

        return out_path

def run_gene_tree_job(gene, alignment_dir, make_gene_tree):
    in_path = os.path.join(alignment_dir, f'{gene}.fasta')

    try:
        return gene, in_path, make_gene_tree(gene), ''
    except Exception as e:
        return gene, in_path, None, str(e)

def write_failed_gene_trees(out_loc, failures):
    out_path = os.path.join(out_loc, 'failed_gene_trees.tsv')

    if not failures:
        if os.path.isfile(out_path):
            os.remove(out_path)

        return

    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['locus', 'alignment', 'error'])
        writer.writerows(failures)

def build_coalescent_tree(args):
    out_loc = args.o.strip()

    if args.phylo_program == 'raxmlng':
        phylo_bin = find_executable('raxml-ng')
    elif args.phylo_program == 'iqtree':
        phylo_bin = find_executable('iqtree')
    elif args.phylo_program == 'veryfasttree':
        phylo_bin = find_executable('VeryFastTree')
    else:
        phylo_bin = find_executable('FastTree')

    astral_bin = find_executable('astral')

    def find_genes(path):
        try:
            with os.scandir(path) as it:
                return {os.path.splitext(entry.name)[0] for entry in it if entry.is_file() and is_reference_file_name(entry.name)}
        except OSError:
            return set()

    alignment_dir = get_locus_alignment_dir(args)

    genes = {t[0] for t in get_ref_genes(args.r)} & find_genes(alignment_dir)
    gene_count = len(genes)

    if not genes:
        raise RuntimeError(f"No gene alignments found under '{alignment_dir}'")

    def make_gene_tree(gene):
        return build_single_tree(args.phylo_program, phylo_bin, os.path.join(alignment_dir, f'{gene}.fasta'), quiet=True)

    tree_files = set()

    failed_gene_trees = []

    def handle_gene_tree_result(result):
        gene, alignment_path, tree_path, error = result

        if error:
            print(f'Warning: gene tree failed on {gene}: {error}')
            failed_gene_trees.append((gene, alignment_path, error))
            return

        if tree_path and os.path.isfile(tree_path):
            tree_files.add(tree_path)
            tree_count = len(tree_files)

            if tree_count >= 2:
                print(f'{tree_count}/{gene_count} trees built\r', end='')

    if args.p > 1:
        with ThreadPoolExecutor(max_workers=args.p) as executor:
            futures = [executor.submit(run_gene_tree_job, gene, alignment_dir, make_gene_tree) for gene in genes]

            for task in as_completed(futures):
                handle_gene_tree_result(task.result())

    else:
        for gene in genes:
            handle_gene_tree_result(run_gene_tree_job(gene, alignment_dir, make_gene_tree))

    print('\n')
    write_failed_gene_trees(out_loc, failed_gene_trees)

    coal_trees_path = os.path.join(out_loc, 'combined_genes.trees')
    coal_out_path = os.path.join(out_loc, 'Coalescent.tree')
    written = False

    with open(coal_trees_path, 'w') as f:
        for path in tree_files:
            if os.path.getsize(path) <= 2: # Empty tree
                continue

            with open(path, 'r') as r:
                tree = next((line.strip() for line in r if line.strip()), '')

            if not tree:
                continue

            f.write(tree + '\n')

            written = True

    if not written:
        raise RuntimeError(f"Unable to reconstruct coalescent trees because no gene tree is available")

    if os.path.isfile(coal_out_path):
        os.remove(coal_out_path)

    subprocess.run([astral_bin, '-i', coal_trees_path, '-o', coal_out_path, '-t', str(args.p)], check=True)

def build_concatenation_tree(args):
    out_loc = args.o.strip()

    if args.phylo_program == 'raxmlng':
        phylo_bin = find_executable('raxml-ng')
    elif args.phylo_program == 'iqtree':
        phylo_bin = find_executable('iqtree')
    elif args.phylo_program == 'veryfasttree':
        phylo_bin = find_executable('VeryFastTree')
    else:
        phylo_bin = find_executable('FastTree')

    in_path = get_concatenated_alignment_path(args)

    if not os.path.isfile(in_path):
        raise RuntimeError(f"Unable to find the concatenated alignment at '{in_path}'")

    final_tree_path = os.path.join(out_loc, 'Concatenation.tree')

    if os.path.isfile(final_tree_path):
        os.remove(final_tree_path)

    out_path = build_single_tree(args.phylo_program, phylo_bin, in_path,
                                 bootstrap=args.bootstrap, threads=args.p)

    if not os.path.isfile(out_path):
        raise RuntimeError(f"Phylogenetic tree reconstruction failed")

    shutil.copyfile(out_path, final_tree_path)

def execute_tasks(args, samples):
    if not os.path.isdir(args.r):
        print(f"Reference directory '{args.r}' does not exist")
        return 2

    commands = frozenset(args.command)

    do_filter = 'filter' in commands
    do_refilter = 'refilter' in commands
    do_assemble = 'assemble' in commands
    do_consensus = 'consensus' in commands
    do_trim = 'trim' in commands
    do_combine = 'combine' in commands
    do_tree = 'tree' in commands
    do_stats = 'stats' in commands

    try:
        if do_filter or do_refilter or do_assemble:
            do_filter_assemble(args, samples, do_filter, do_refilter, do_assemble)

            if do_assemble and args.assembly_mode == 'uce':
                write_uce_outputs(args, samples)

        if do_consensus:
            generate_consensus(args, samples)

        if do_trim:
            if not args.trim_source:
                args.trim_source = 'consensus' if do_consensus else 'assembly'

            blast_trim(args, samples)

        if do_combine:
            if not args.combine_source:
                if do_trim:
                    args.combine_source = 'trimmed'
                elif do_consensus:
                    args.combine_source = 'consensus'
                else:
                    args.combine_source = 'assembly'

            combine_genes(args, samples)

        if do_tree:
            if args.tree_method == 'coalescent':
                build_coalescent_tree(args)
            else:
                build_concatenation_tree(args)

        if do_stats:
            gm2_stats.run(args, samples)

    except (RuntimeError, subprocess.SubprocessError) as e:
        print(f'Error: {e}')
        return 1

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='GeneMiner2-UCE extracts phylogenetic marker loci for UCE workflows.',
                                     epilog=HELP_EPILOG)
    parser.add_argument('command',
                        choices=('filter', 'refilter', 'assemble', 'consensus', 'trim', 'combine', 'tree', 'stats', []),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='command',
                        nargs='*')

    group_io = parser.add_argument_group('input/output parameters')
    group_io.add_argument('-f', help='Sample list file', metavar='FILE', required=True)
    group_io.add_argument('-r', help='Reference directory', metavar='DIR', required=True)
    group_io.add_argument('-o', help='Output directory', metavar='DIR', required=True)
    group_io.add_argument('-p', default=1, help='Number of parallel processes', metavar='INT', type=int)

    group_filter = parser.add_argument_group('arguments for filtering')
    group_filter.add_argument('-kf', default=31, help='Filter k-mer size', metavar='INT', type=int)
    group_filter.add_argument('-s', '--step-size', default=4, help='Filter step size', metavar='INT', type=int)
    group_filter.add_argument('--max-reads', default=0, help='Million reads to process per file', metavar='INT', type=int)
    group_filter.add_argument('--reuse-reference-cache', action='store_true', default=False, help='Reuse a fingerprinted reference k-mer cache instead of rebuilding it every run')
    group_filter.add_argument('--reference-cache-dir', default=None, help='Directory for --reuse-reference-cache files (default = output/.gm2_reference_cache)', metavar='DIR')

    group_refilter = parser.add_argument_group('arguments for futher filtering')
    group_refilter.add_argument('--depth-low-water-mark', default=50, help='If depth is lower than this value, try to find more reads with relaxed criteria', metavar='INT', type=int)
    group_refilter.add_argument('--depth-limit', default=768, help='Maximum depth processed during re-filtering', metavar='INT', type=int)
    group_refilter.add_argument('--file-size-limit', default=6, help='Maximum file size during re-filtering', metavar='INT', type=int)

    group_assembly = parser.add_argument_group('arguments for assembling')
    group_assembly.add_argument('-ka', default=0, help='Assembly k-mer size (default = auto)', metavar='INT', type=int)
    group_assembly.add_argument('--min-ka', default=21, help='Minimum auto-estimated assembly k-mer size', metavar='INT', type=int)
    group_assembly.add_argument('--max-ka', default=51, help='Maximum auto-estimated assembly k-mer size', metavar='INT', type=int)
    group_assembly.add_argument('-e', '--error-threshold', default=2, help='Error threshold', metavar='INT', type=int)
    group_assembly.add_argument('-sb', '--soft-boundary', default='auto', help='Soft boundary (default = auto)', metavar='{INT,auto,unlimited}', type=str)
    group_assembly.add_argument('-i', '--search-depth', default=4096, help='Search depth', metavar='INT', type=int)
    group_assembly.add_argument('--min-coverage', default=0, help='Minimum read depth required for contig generation', metavar='INT', type=int)
    group_assembly.add_argument('--assembly-mode', choices=('reference', 'uce'), default='reference', help='Assembly mode: reference keeps the default reference-guided behavior; uce preserves read-supported flanks around conserved cores')
    group_assembly.add_argument('--uce-side-candidates', default=8, help='One-sided branch candidates to combine in UCE mode (default = 8)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-max-contig-length', default=5000, help='Maximum UCE contig length kept before scoring; use 0 to disable (default = 5000)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-min-read-density', default=0.003, help='Minimum uniquely placed read_count/length for long UCE contigs before scoring (default = 0.003)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-density-check-min-length', default=1000, help='Minimum contig length where the UCE read-density guardrail applies (default = 1000)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-max-depth-cv', default=0, help='Optional maximum k-mer depth coefficient of variation for UCE contigs; 0 disables (default = 0)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-max-depth-ratio', default=0, help='Optional maximum max/median k-mer depth ratio for UCE contigs; 0 disables (default = 0)', metavar='FLOAT', type=float)
    group_assembly.add_argument('--uce-rescue-reads', action='store_true', default=False, help='UCE mode only: after the first assembly, recruit raw reads once using preliminary contigs plus original references, then re-filter and re-assemble')
    group_assembly.add_argument('--uce-rescue-min-contig-length', default=60, help='Minimum preliminary contig length used as a UCE raw-read rescue reference (default = 60)', metavar='INT', type=int)
    group_assembly.add_argument('--uce-rescue-min-density-ratio', default=0.5, help='Minimum rescue/before read-density ratio kept after UCE raw-read rescue (default = 0.5)', metavar='FLOAT', type=float)

    group_consensus = parser.add_argument_group('argument for consensus generation')
    group_consensus.add_argument('-c', '--consensus-threshold', default='0.75', help='Consensus threshold (default = 0.75)', metavar='FLOAT', type=float)

    group_trim = parser.add_argument_group('arguments for sequence trimming')
    group_trim.add_argument('-ts', '--trim-source', choices=('assembly', 'consensus'), default=None, help='Whether to trim the primary assembly or the consensus sequence (default = output of last step, assembly if no other command given)')
    group_trim.add_argument('-tm', '--trim-mode', choices=('all', 'longest', 'terminal', 'isoform'), default='terminal', help='Trim mode (default = terminal)', type=str)
    group_trim.add_argument('-tr', '--trim-retention', default=0, help='Retention length threshold (default = 0.0)', metavar='FLOAT', type=float)

    group_combine = parser.add_argument_group('arguments for sequence alignment and clustering')
    group_combine.add_argument('-cs', '--combine-source', choices=('assembly', 'consensus', 'trimmed'), default=None, help='Whether to combine the primary assembly, the consensus sequences or the trimmed sequences (default = output of last step, assembly if no other command given)')
    group_combine.add_argument('-cd', '--clean-difference', default=1, help='Maximum acceptable pairwise difference in an alignment (default = 1.0)', metavar='FLOAT', type=float)
    group_combine.add_argument('-cn', '--clean-sequences', default=0, help='Number of sequences required in an alignment (default = 0)', metavar='INT', type=int)
    group_combine.add_argument('--msa-program', choices=('clustalo', 'mafft', 'muscle'), default='mafft', help='Program for multiple sequence alignment', type=str)
    group_combine.add_argument('--msa-threads', default=1, help='Threads used by each multiple-sequence-alignment job (default = 1)', metavar='INT', type=int)
    group_combine.add_argument('--alignment-filter', choices=('trimal', 'alifilter', 'none'), default=None, help='Program for filtering aligned loci before tree reconstruction (default = trimal)', type=str)
    group_combine.add_argument('--filter-processes', default=None, help='Maximum number of concurrent alignment filtering jobs (default = -p)', metavar='INT', type=int)
    group_combine.add_argument('--alifilter-model', default=None, help='AliFilter model specification or model.json path when --alignment-filter alifilter is used', metavar='MODEL', type=str)
    group_combine.add_argument('--strict-combine-errors', action='store_true', default=False, help='Stop combine if any locus fails during alignment, cleanup, or alignment filtering')
    group_combine.add_argument('--no-alignment', action='store_true', default=False, help='Do not perform multiple sequence alignment')
    group_combine.add_argument('--no-trimal', action='store_true', default=False, help='Do not run alignment filtering (deprecated; use --alignment-filter none)')

    group_tree = parser.add_argument_group('argument for tree inference')
    group_tree.add_argument('-m', '--tree-method', choices=('coalescent', 'concatenation'), default='coalescent', help='Multi-gene tree reconstruction method (default = coalescent)')
    group_tree.add_argument('-b', '--bootstrap', default=1000, help='Number of bootstrap replicates', metavar='INT', type=int)
    group_tree.add_argument('--phylo-program', choices=('raxmlng', 'iqtree', 'fasttree', 'veryfasttree'), default='fasttree', help='Program for phylogenetic tree reconstruction', type=str)

    group_stats = parser.add_argument_group('arguments for UCE statistics')
    group_stats.add_argument('--stats-no-heatmap', action='store_true', default=False, help='Do not create UCE statistics heatmaps')
    group_stats.add_argument('--stats-count-input-reads', action='store_true', default=False, help='Count input FASTQ reads for InputReads and PctFiltered statistics; can be slow for large datasets')

    parser.add_argument('--min-depth', help=argparse.SUPPRESS)
    parser.add_argument('--max-depth', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.reference_cache_dir and not args.reuse_reference_cache:
        parser.error('--reference-cache-dir requires --reuse-reference-cache')

    args.uce_side_candidates = max(args.uce_side_candidates, 3)
    args.uce_max_contig_length = max(args.uce_max_contig_length, 0)
    args.uce_density_check_min_length = max(args.uce_density_check_min_length, 1)
    args.uce_rescue_min_contig_length = max(args.uce_rescue_min_contig_length, args.kf)

    if args.uce_min_read_density < 0:
        parser.error('--uce-min-read-density must be greater than or equal to 0')

    if args.uce_max_depth_cv < 0:
        parser.error('--uce-max-depth-cv must be greater than or equal to 0')

    if args.uce_max_depth_ratio < 0:
        parser.error('--uce-max-depth-ratio must be greater than or equal to 0')

    if args.uce_rescue_min_density_ratio <= 0:
        parser.error('--uce-rescue-min-density-ratio must be greater than 0')

    if args.no_trimal and args.alignment_filter not in (None, 'none'):
        parser.error('--no-trimal cannot be combined with --alignment-filter trimal or --alignment-filter alifilter')

    if get_alifilter_model(args) and get_alignment_filter(args) != 'alifilter':
        parser.error('--alifilter-model requires --alignment-filter alifilter')

    if args.p < 1:
        parser.error('-p must be at least 1')

    if args.msa_threads < 1:
        parser.error('--msa-threads must be at least 1')

    if args.msa_threads > args.p:
        parser.error('--msa-threads cannot be greater than -p')

    if args.filter_processes is not None and args.filter_processes < 1:
        parser.error('--filter-processes must be at least 1')

    if not args.command:
        if args.assembly_mode == 'uce':
            args.command = ('filter', 'refilter', 'assemble', 'combine', 'tree')
        else:
            args.command = ('filter', 'refilter', 'assemble', 'trim', 'combine', 'tree')

    if args.min_depth is not None:
        print('  Option --min-depth has been removed. Please use --depth-low-water-mark, --error-threshold or --min-coverage instead.')
        print(DEPTH_DEPRECATION_EXPLAINER)
        sys.exit(2)

    if args.max_depth is not None:
        print('  Option --max-depth has been removed. Please use --depth-limit instead.')
        print(DEPTH_DEPRECATION_EXPLAINER)
        sys.exit(2)

    samples = prepare_workdir(args)

    if samples:
        print(f'Running tasks: {", ".join(args.command)}')
        print()
        sys.exit(execute_tasks(args, samples))
    else:
        print('Sample list is empty or invalid, exiting')
        sys.exit(2)
