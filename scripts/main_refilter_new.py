from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from main_refilter_ext import collect_runs_stats, filter_read, parse_record
import argparse
import collections
import contextlib
import itertools
import math
import multiprocessing
import os
import shutil
import struct

FILE_EXTENSION = {
    'fasta': '.fasta',
    'fastq': '.fq'
}

FILE_TYPES = {
    '.fa': 'fasta',
    '.fas': 'fasta',
    '.fasta': 'fasta',
    '.fq': 'fastq',
    '.fastq': 'fastq',
    '.gm2': 'fastq'
}

FORMAT_FUNCTIONS = {
    'fasta': lambda t: f'>{t[0]}\n{t[1]}\n',
    'fastq': lambda t: f'@{t[0]}\n{t[1]}\n+\n{t[2]}\n'
}

READ_ITERATORS = {
   'fasta': SimpleFastaParser,
   'fastq': FastqGeneralIterator
}

BASE_TO_INT = {
    'A': '0', 'C': '1', 'G': '2', 'T': '3', 'U': '3',
    'a': '0', 'c': '1', 'g': '2', 't': '3', 'u': '3',
}

def encode_sequence(seq):
    """Encode a read only when every base is unambiguous.

    Removing ambiguous bases would join its flanks and create artificial
    k-mers, so ambiguous reads are conservatively excluded from matching.
    """
    encoded = []
    for base in seq:
        value = BASE_TO_INT.get(base)
        if value is None:
            return None
        encoded.append(value)
    return ''.join(encoded)

def encode_reference_runs(seq):
    """Yield independent A/C/G/T/U runs without k-mers crossing ambiguity."""
    run = []
    for base in seq:
        value = BASE_TO_INT.get(base)
        if value is None:
            if run:
                yield ''.join(run)
                run = []
            continue
        run.append(value)
    if run:
        yield ''.join(run)

def print_log(log_path, *args, **kwargs):
    if log_path:
        with open(log_path, 'a') as out:
            print(*args, file=out, **kwargs)

    print(*args, **kwargs)

def get_read_dict(se_dir, pe_dir, use_compressed):
    read_dict = {}
    walk_directory = lambda path: ((os.path.dirname(ent.path), *os.path.splitext(ent.name)) for ent in os.scandir(path) if ent.is_file())

    if se_dir:
        if not os.path.isdir(se_dir):
            raise ValueError('Argument --se-dir does not refer to a directory.')

        for dirname, basename, extname in walk_directory(se_dir):
            if extname not in FILE_TYPES:
                continue

            if not use_compressed and extname == '.gm2':
                continue

            if basename in read_dict:
                if os.path.splitext(read_dict[basename][0])[1] == '.gm2':
                    continue

                if extname != '.gm2':
                    raise ValueError(f'Duplicate read group name {basename}.')

            read_dict[basename] = (f'{dirname}/{basename}{extname}', )

    if pe_dir:
        if not os.path.isdir(pe_dir):
            raise ValueError('Argument --pe-dir does not refer to a directory.')

        for dirname, basename, extname in walk_directory(pe_dir):
            if extname not in FILE_TYPES or not basename.endswith('_1'):
                continue

            gene_name = basename[:-2]

            if not use_compressed and extname == '.gm2':
                continue

            if gene_name in read_dict:
                if os.path.splitext(read_dict[gene_name][0])[1] == '.gm2':
                    continue

                if extname != '.gm2':
                    raise ValueError(f'Duplicate read group name {gene_name}.')

            forward_path = os.path.join(dirname, f'{gene_name}_1{extname}')
            reverse_path = os.path.join(dirname, f'{gene_name}_2{extname}')

            if os.path.isfile(reverse_path):
                read_dict[gene_name] = (forward_path, reverse_path)
            else:
                read_dict[gene_name] = (forward_path, )

    return read_dict

def get_ref_dict(ref_dir):
    ref_dict = {}

    with os.scandir(ref_dir) as entries:
        for ent in entries:
            basename, extname = os.path.splitext(ent.name)

            if extname not in FILE_TYPES:
                continue

            if basename in ref_dict:
                raise ValueError(f'Duplicate reference sequence name {basename}.')

            ref_dict[basename] = ent.path

    return ref_dict

def load_reference(ref_path, kmer_size):
    with open(ref_path, 'r') as f:
        ref_set = {seq for _, seq in SimpleFastaParser(f) if len(seq) >= kmer_size}

    if not ref_set:
        return ref_set, 0

    length_list = list(map(len, ref_set))
    effective_len = int(max(length_list) * (math.log10(len(length_list)) + 1))
    return ref_set, effective_len

def gm2_iterator(f, suffix=''):
    read_id = 0
    seq_buf = bytearray(1024)
    phr_buf = bytearray(1024)

    while True:
        rec_hdr = f.read(6)

        if len(rec_hdr) < 6:
            break

        rl1, rl2, sl1, sl2 = struct.unpack('!BHBH', rec_hdr)
        rec_len = (rl1 << 16) | rl2
        has_phr = (sl1 & 0x80) != 0
        seq_len = ((sl1 & 0x7f) << 16) | sl2

        if rec_len == 0:
            continue

        if len(seq_buf) < seq_len:
            seq_buf = bytearray(seq_len)
            phr_buf = bytearray(seq_len)

        record = f.read(rec_len)

        if seq_len == 0:
            continue

        parse_record(record, has_phr, seq_buf, phr_buf, seq_len)
        read_id += 1

        yield (f'read_{read_id}{suffix}', seq_buf[:seq_len].decode("ascii"), phr_buf[:seq_len].decode("ascii") if has_phr else '')

def linked_iterator(iterator, link_size):
    while True:
        linked_reads = []

        for _ in range(link_size):
            try:
                linked_reads.append(next(iterator))
            except StopIteration:
                break

        if not linked_reads:
            break

        if len(linked_reads) != link_size:
            raise ValueError('Interleaved paired-read file has an odd number of records.')

        yield tuple(linked_reads)

def linked_read_iterators(read_iters):
    for linked_reads in itertools.zip_longest(*read_iters, fillvalue=None):
        if any(read is None for read in linked_reads):
            raise ValueError('Paired input files contain different numbers of records.')
        yield linked_reads

def build_kmer_dict(ref_set, kmer_size):
    # Values: 0=unused; 1=forward; 2=reverse; 3=both
    kmer_dict = collections.defaultdict(lambda: 0)
    mask_bin  = (1 << (kmer_size << 1)) - 1

    for seq in ref_set:
        for seq_str in encode_reference_runs(seq):
            if len(seq_str) < kmer_size:
                continue

            seq_str_r = ''.join(str(3 - int(base)) for base in reversed(seq_str))
            seq_int   = int(seq_str, 4)
            seq_int_r = int(seq_str_r, 4)

            for _ in range(0, len(seq_str) - kmer_size + 1):
                kmer_dict[seq_int   & mask_bin] |= 1
                kmer_dict[seq_int_r & mask_bin] |= 2
                seq_int   >>= 2
                seq_int_r >>= 2

    return kmer_dict

def copy_reads(name, out_dir, read_info, file_type):
    gm2_format  = os.path.splitext(read_info[0])[1] == '.gm2'
    output_ext  = FILE_EXTENSION[file_type]
    output_path = os.path.join(out_dir, name + output_ext)
    format_func = FORMAT_FUNCTIONS[file_type]

    with contextlib.ExitStack() as stack:
        if gm2_format:
            read_iters = [gm2_iterator(stack.enter_context(open(path, 'rb')), f'/{i}') for i, path in enumerate(read_info, start=1)]
        else:
            read_iters = [READ_ITERATORS[file_type](stack.enter_context(open(path, 'r'))) for path in read_info]

        output_file = stack.enter_context(open(output_path, 'w'))
        output_file.writelines(format_func(tp) for linked_reads in linked_read_iterators(read_iters) for tp in linked_reads)

def run_length_filter(name, out_dir, ref_set, ref_length, read_info, file_type, kmer_size, keep_temporaries, keep_linked_mates):
    RUN_LEN_CONST = 0.5772156649 / math.log(2) - 1.5
    THR_P95_2T = 1.96
    THR_1e5_1T = 3.74
    TOLERANCE = 1e-5

    gm2_format  = os.path.splitext(read_info[0])[1] == '.gm2'
    output_ext  = FILE_EXTENSION[file_type]
    output_path = os.path.join(out_dir, 'large_files', name + output_ext)
    format_func = FORMAT_FUNCTIONS[file_type]
    open_flags  = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    if os.name == 'nt' and not keep_temporaries:
        open_flags |= os.O_SHORT_LIVED

    kmer_dict = build_kmer_dict(ref_set, kmer_size)

    with contextlib.ExitStack() as stack:
        if gm2_format:
            read_iters = [gm2_iterator(stack.enter_context(open(path, 'rb')), f'/{i}') for i, path in enumerate(read_info, start=1)]
        else:
            read_iters = [READ_ITERATORS[file_type](stack.enter_context(open(path, 'r'))) for path in read_info]

        output_file = stack.enter_context(os.fdopen(os.open(output_path, open_flags), 'w'))

        for linked_reads in linked_read_iterators(read_iters):
            orient = [0] * len(linked_reads)
            group_iter = (
                collect_runs_stats(encoded, kmer_dict, kmer_size) if (encoded := encode_sequence(tp[1])) is not None else [0] * 13
                for tp in linked_reads
            )

            for i, (_, fwd_l, rev_l, _,
                    _, fwd_r, rev_r, _,
                    _, fwd_n, rev_n, amb_n, tot_n) in enumerate(group_iter):

                # Forward hits    Reverse hits    Ambiguous hits    Verdict
                # -----------------------------------------------------------
                # <= 1            <= 1            <= 1              reject
                # <= 1            <= 1            > 1               ambiguous
                # <= 1            > 1             *                 reverse
                # > 1             <= 1            *                 forward
                # > 1             > 1             *                 continue
                if fwd_n <= 1:
                    orient[i] = (rev_n <= 1) * (1 - (amb_n <= 1) * 3) + 2
                    continue
                elif rev_n <= 1:
                    orient[i] = 1
                    continue

                # Note that we count the runs for all four states
                # (mismatch, forward, reverse, ambiguous)
                # but we calculate the expected number of runs with two states
                # This is equivalent to spliting a run into two whenever
                # a mismatch is encountered
                # In principle, E(R) = 2*n1*n2/(n1+n2)+mutation_rate*(n1+n2)+1
                # However, we choose to ignore the mutation term and tolerate
                # a bounded multiplier on E(R), implying much higher stringency
                npr = 2 * fwd_n * rev_n
                nht = fwd_n + rev_n

                # E(R) in standard runs test
                # 2*n1*n2/(n1+n2)+1
                erc = npr / nht + 1

                # Var(R) in standard runs test
                # 2*n1*n2*(2*n1*n2-n1-n2)/(n1+n2)^2*(n1+n2-1)
                vrn = npr * (npr - nht) / (nht * nht * (nht - 1))

                # The direction of runs does not autocorrelate enough
                # i.e. the runs are either random or somehow periodic
                if (fwd_r + rev_r - erc) / math.sqrt(vrn) > -THR_1e5_1T:

                    # Next, we assume some mismatches caused the difference of
                    # the numbers of runs
                    # Let pf be the error rate in the forward matching region
                    # pf = (fwd_r - true_fwd_r) / (fwd_n / (1 - pf))
                    if fwd_r > rev_r:
                        ntt = fwd_n + fwd_r - rev_r
                        rex = fwd_n / ntt
                    else:
                        ntt = fwd_n + rev_r - fwd_r
                        rex = fwd_n / ntt

                    # If we cannot infer that the direction with more runs has
                    # a significant proportion of mismatches against reference
                    # then call a chimera
                    if math.isclose(rex, 1.0, abs_tol=TOLERANCE) or (1 - rex) / math.sqrt(rex * (1 - rex) / ntt) < THR_P95_2T:
                        orient[i] = 0
                        continue

                erl = max(math.log2(tot_n) + RUN_LEN_CONST, 0) + 4
                orient[i] = (fwd_l > erl) + (rev_l > erl) * 2

                # The orientation is unambiguous
                if orient[i] != 3:
                    continue

                # If the longest forward and reverse runs are both much longer
                # than the expected run-length, assume all forward and reverse
                # matches are contiguous and calculate an approximate mutation
                # rate (https://math.stackexchange.com/a/5027331)
                # If the mutation rates are too similar, we consider matches in
                # both directions 'similarly good', thereby calling a chimera
                # As a rule of thumb, we reject reads if the forward region and
                # the reverse region share the same distribution, because their
                # chimeric appearance impedes reference-guided assembly
                lpf = math.exp(math.log(1 / (1 - fwd_l + fwd_n)) / fwd_l)
                lpr = math.exp(math.log(1 / (1 - rev_l + rev_n)) / rev_l)
                fpz = math.isclose(lpf, 0.0, abs_tol=TOLERANCE)
                rpz = math.isclose(lpr, 0.0, abs_tol=TOLERANCE)

                if fpz:
                    orient[i] = 2 - 2 * rpz
                elif rpz:
                    orient[i] = 1
                elif math.isclose(lpf, 1.0, abs_tol=TOLERANCE) and math.isclose(lpr, 1.0, abs_tol=TOLERANCE):
                    orient[i] = 0
                elif abs(lpf - lpr) / math.sqrt(lpf ** 2 * (1 - lpf) / fwd_n + lpr ** 2 * (1 - lpr) / rev_n) < THR_P95_2T:
                    orient[i] = 0

            # Discordant paired reads
            if len(orient) == 2 and 1 <= orient[0] <= 2 and orient[0] == orient[1]:
                continue

            if keep_linked_mates and len(linked_reads) == 2 and any(orient):
                output_file.writelines(format_func(tp) for tp in linked_reads)
            else:
                output_file.writelines(format_func(tp) for i, tp in enumerate(linked_reads) if orient[i])

    return output_path

def kmer_filter(name, out_dir, log_path, ref_set, ref_length, temp_path, file_type, kmer_size, min_depth, max_depth, max_size, keep_temporaries, keep_linked_mates):
    output_ext  = FILE_EXTENSION[file_type]
    output_path = os.path.join(out_dir, name + output_ext)
    format_func = FORMAT_FUNCTIONS[file_type]
    read_iter   = READ_ITERATORS[file_type]
    def read_matches(tp, kmer_dict):
        encoded = encode_sequence(tp[1])
        return encoded is not None and filter_read(encoded, kmer_dict, kmer_size)

    with open(temp_path, 'r') as f:
        if keep_linked_mates:
            total_length = sum(sum(len(tp[1]) for tp in linked_reads)
                               for linked_reads in linked_iterator(read_iter(f), 2))
        else:
            total_length = sum(len(tp[1]) for tp in read_iter(f))

    coverage  = total_length / ref_length
    too_deep  = coverage > max_depth
    too_large = total_length // 1e6 > max_size

    if not too_deep and not too_large:
        return shutil.copyfile(temp_path, output_path)

    min_depth = min(min_depth, max_depth / 4)
    initial_kmer_size = kmer_size

    # NOTE: the largest possible k-mer size is 63 + 6 = 69
    while kmer_size < 64 and (too_deep or too_large):
        last_kmer_size = kmer_size
        last_length = total_length

        if coverage > 8 * max_depth or total_length // 1e6 > 6 * max_size:
            kmer_size += 6
        else:
            kmer_size += 2

        print_log(log_path, f'K-mer size for {name}: {kmer_size}')

        kmer_dict = build_kmer_dict(ref_set, kmer_size)

        with open(temp_path, 'r') as f:
            if keep_linked_mates:
                total_length = sum(sum(len(tp[1]) for tp in linked_reads)
                                   for linked_reads in linked_iterator(read_iter(f), 2)
                                   if any(read_matches(tp, kmer_dict) for tp in linked_reads))
            else:
                total_length = sum(len(tp[1])
                                   for tp in read_iter(f)
                                   if read_matches(tp, kmer_dict))

        coverage  = total_length / ref_length
        too_deep  = coverage > max_depth
        too_large = total_length // 1e6 > max_size

        if coverage < min_depth:
            kmer_size    = last_kmer_size
            total_length = last_length
            coverage     = total_length / ref_length
            too_large    = total_length // 1e6 > max_size
            break

    if kmer_size == initial_kmer_size and not too_large:
        return shutil.copyfile(temp_path, output_path)

    kmer_dict = build_kmer_dict(ref_set, kmer_size)
    interval = max(int(total_length / 1e6 / max_size), 2)
    i = 0

    with open(temp_path, 'r') as f, open(output_path, 'w') as fo:
        if keep_linked_mates:
            for linked_reads in linked_iterator(read_iter(f), 2):
                if any(read_matches(tp, kmer_dict) for tp in linked_reads):
                    i += 1

                    if too_large and i % interval != 0:
                        continue

                    fo.writelines(format_func(tp) for tp in linked_reads)
        else:
            for tp in read_iter(f):
                if read_matches(tp, kmer_dict):
                    i += 1

                    if too_large and i % interval != 0:
                        continue

                    fo.write(format_func(tp))

def filter_gene(task):
    file_ext = os.path.splitext(task.read_path[0])[1]

    if file_ext not in FILE_TYPES:
        print_log(task.log_path, f"File '{task.read_path[0]}' has invalid file type.")
        return

    file_type = FILE_TYPES[file_ext]

    if task.copy_only:
        print_log(task.log_path, f'Writing reads for gene {task.name}.')
        copy_reads(task.name, task.out_dir, task.read_path, file_type)
        return

    print_log(task.log_path, f'Filtering gene {task.name}.')

    ref_set, effective_len = load_reference(task.ref_path, task.kmer_size)

    if not effective_len:
        print_log(task.log_path, f'Gene {task.name} has no valid reference.')
        return

    # On the weird choice of k-mer size
    # Assume the sequencing error rate to be 0.95
    # 1 - 0.95^13 = 0.4867 < 0.5
    # On average, one of two clusters of biological k-mer matches is error-free
    tmp_path = run_length_filter(task.name, task.out_dir, ref_set, effective_len,
                                 task.read_path, file_type, max(task.kmer_size // 2, task.kmer_size - 13) | 1,
                                 task.keep_temporaries, task.keep_linked_mates and len(task.read_path) == 2)

    kmer_filter(task.name, task.out_dir, task.log_path,
                ref_set, effective_len, tmp_path, file_type,
                task.kmer_size, task.min_depth, task.max_depth,
                task.max_size, task.keep_temporaries, task.keep_linked_mates and len(task.read_path) == 2)

    if not task.keep_temporaries:
        os.unlink(tmp_path)

Task = collections.namedtuple('Task', ('name', 'out_dir', 'ref_path', 'read_path',
                                       'log_path', 'min_depth', 'max_depth', 'max_size',
                                       'copy_only', 'keep_temporaries', 'keep_linked_mates',
                                       'kmer_size'))

def run(args):
    tasks = []

    for name, ref_path in ref_dict.items():
        if name not in read_dict:
            print_log(args.log_file, f"No reads for gene {name}.")
            continue

        tasks.append(Task(name, out_dir, ref_path, read_dict[name], args.log_file,
                          args.min_depth, args.max_depth, args.max_size,
                          args.copy_only, args.keep_temporaries, args.keep_linked_mates,
                          args.kmer_size))

    if not tasks:
        print_log(args.log_file, 'No genes with matching reads and references.')
        return

    processes = min(args.processes, len(tasks))

    if processes > 1:
        chunksize = max(1, len(tasks) // (processes * 4))

        with multiprocessing.Pool(processes) as pool:
            for _ in pool.imap_unordered(filter_gene, tasks, chunksize=chunksize):
                pass
    else:
        for task in tasks:
            filter_gene(task)

    if not args.keep_temporaries:
        try:
            os.rmdir(os.path.join(out_dir, 'large_files'))
        except OSError:
            pass

if __name__ == '__main__':
    if os.name == 'nt':
        multiprocessing.freeze_support()

    parser = argparse.ArgumentParser(description='An improved NGS filtering gadget based on k-mers.')

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-qs', '--se-dir', help='Directory with single-read sequencing data')
    input_group.add_argument('-qd', '--pe-dir', help='Directory with paired-end sequencing data')

    parser.add_argument('-r', '--ref-dir', required=True, help='Directory with reference sequences')
    parser.add_argument('-o', '--out-dir', required=True, help='Output directory')

    parser.add_argument('--log-file', default=None, help='Log file')
    parser.add_argument('--min-depth', default=50, help='Min allowed coverage', type=int)
    parser.add_argument('--max-depth', default=768, help='Max allowed coverage', type=int)
    parser.add_argument('--max-size', default=6, help='Max allowed size in million bases', type=int)
    parser.add_argument('--copy-only', action='store_true', help='Interleave paired reads without filtering')
    parser.add_argument('--keep-temporaries', action='store_true', help='Keep temporary files')
    parser.add_argument('--keep-linked-mates', action='store_true', help='For paired-end reads, keep both mates when either mate passes filtering')
    parser.add_argument('--use-gm2-format', action='store_true', help='Read reads from compressed binary format')
    parser.add_argument('-kf', '--kmer-size', default=31, help='K-mer size', type=int)

    parser.add_argument('-p', '--processes', default=1, help='Number of parallel processes', type=int)

    args = parser.parse_args()

    if args.kmer_size < 1:
        parser.error('--kmer-size must be positive')

    if args.min_depth < 0:
        parser.error('--min-depth must be zero or positive')

    if args.max_depth <= 0:
        parser.error('--max-depth must be positive')

    if args.max_size <= 0:
        parser.error('--max-size must be positive')

    if args.processes < 1:
        parser.error('--processes must be positive')

    try:
        read_dict = get_read_dict(args.se_dir, args.pe_dir, args.use_gm2_format)
        ref_dict = get_ref_dict(args.ref_dir)
        out_dir = args.out_dir
        os.makedirs(os.path.join(out_dir, 'large_files'), exist_ok=True)
    except (OSError, ValueError) as e:
        parser.error(str(e))
    else:
        run(args)
