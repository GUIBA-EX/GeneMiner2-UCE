from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter, deque
from itertools import chain
from operator import itemgetter
import argparse
import csv
import gc
import hashlib
import math
import multiprocessing
import os
import pickle
import sys
import time
from typing import NamedTuple

D_BASE_DICT = {'AG':'R','CT':'Y', 'GT':'K', 'GC':'S','AC':'M', 'AT':'W','GA':'R','TC':'Y','TG':'K', 'CG':'S','CA':'M', 'TA':'W',}
ACGT_DICT = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
ACGT_REV   = str.maketrans('ACGT', 'TGCA')
FWD_TRANS  = str.maketrans("ACGTU", "01233")
REV_TRANS  = str.maketrans("ACGTU", "32100")
BIN_DICT   = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}

ref_path_dict = {}  # 序列路径字典
ref_count_dict = {} # 参考序列条数字典
kmer_dict = {}  # kmer字典
ref_reads_count_dict = {}  # reads计数的字典

RAW_SEQ = 0
RAW_WEIGHT = 1
RAW_READ_COUNT = 2
RAW_SUPPORTED_SPAN = 3
RAW_FLANK_BALANCE = 4
RAW_READ_DENSITY = 5
RAW_SUPPORT_FRACTION = 6
RAW_KMER_MEDIAN_DEPTH = 7
RAW_KMER_DEPTH_CV = 8
RAW_KMER_MAX_DEPTH_RATIO = 9
RAW_UNIQUE_READ_COUNT = 10
RAW_MULTI_MAPPING_READ_COUNT = 11
RAW_SUPPORTED_BASES = 12
RAW_SUPPORT_BREADTH = 13
RAW_MAX_SUPPORT_GAP = 14
RAW_ACCEPTED = 15
RAW_REJECTION_REASON = 16

CONTIG_SEQ = 0
CONTIG_SEED_COUNT = 1
CONTIG_POS = 2
CONTIG_WEIGHT = 3
CONTIG_READ_COUNT = 4
CONTIG_SUPPORTED_SPAN = 5
CONTIG_FLANK_BALANCE = 6
CONTIG_READ_DENSITY = 7
CONTIG_SUPPORT_FRACTION = 8
CONTIG_KMER_MEDIAN_DEPTH = 9
CONTIG_KMER_DEPTH_CV = 10
CONTIG_KMER_MAX_DEPTH_RATIO = 11
CONTIG_UNIQUE_READ_COUNT = 12
CONTIG_MULTI_MAPPING_READ_COUNT = 13
CONTIG_SUPPORTED_BASES = 14
CONTIG_SUPPORT_BREADTH = 15
CONTIG_MAX_SUPPORT_GAP = 16
CONTIG_ACCEPTED = 17
CONTIG_REJECTION_REASON = 18


class ReadSupport(NamedTuple):
    total_read_count: int
    unique_read_count: int
    multi_mapping_read_count: int
    supported_extent: int
    supported_bases: int
    breadth: float
    max_gap: int
    flank_balance: float
    left_coord: int
    right_coord: int

def Write_Print(log_path, *log_str, sep = " "):
    """
    记录日志并打印
    """
    line = sep.join(map(str,log_str)).strip()
    with open(log_path, 'a') as out:
        out.write(line + '\n')
    print(line)

def Seq_To_Int(dna_str, trans=FWD_TRANS, rtrans=REV_TRANS):
    """
    将基因转换为整数
    """
    # Do not remove ambiguous bases: doing so joins the two flanks and creates
    # k-mers that never occurred in the input sequence.
    if any(base not in 'ACGTU' for base in dna_str):
        return (), 0

    dna_fw_str = dna_str.translate(trans)
    dna_rc_str = dna_str.translate(rtrans)[::-1]

    if not dna_fw_str:
        return (), 0

    return (int(dna_fw_str, 4), int(dna_rc_str, 4)), len(dna_fw_str)

def Valid_Sequence_Runs(dna_str):
    """Yield contiguous unambiguous sequence runs without joining across ambiguity."""
    for _, run in Valid_Sequence_Runs_With_Positions(dna_str):
        yield run

def Valid_Sequence_Runs_With_Positions(dna_str):
    """Yield the start coordinate and sequence of every unambiguous run."""
    run = []
    run_start = 0
    for pos, base in enumerate(dna_str):
        if base in 'ACGTU':
            if not run:
                run_start = pos
            run.append(base)
        elif run:
            yield run_start, ''.join(run)
            run = []
    if run:
        yield run_start, ''.join(run)

def Int_To_Seq(seq_bin, seq_length, seq_dict=BIN_DICT):
    """
    将整数转换为基因
    :param seq_bin: 序列的整数形式
    :param seq_length: 序列的长度，用来补0
    :return: 返回整数序列的基因形式
    """
    seq_bin_str = bin(seq_bin)[2:].rjust(seq_length << 1, '0')
    return ''.join(seq_dict[seq_bin_str[j << 1:(j << 1) + 2]] for j in range(seq_length))

def Reverse_Int(dna_int, dna_length):
    """
    将基因的整数型转换为反向互补的基因的整数型
    :param dna_int: 序列的整数形式
    :param dna_length: 序列的长度，用来补0
    :return: 返回整数序列的基因形式
    """
    bin_str = bin(dna_int ^ ((1 << (dna_length << 1)) - 1))[2:].rjust(dna_length << 1, '0')
    new_list = [bin_str[i:i + 2] for i in range(0, dna_length << 1, 2)]
    return int(''.join(reversed(new_list)), 2)

def Reverse_Complement_ACGT(seq, table=ACGT_REV):
    """
    简化版反向互补
    """ 
    return seq.translate(table)[::-1]

def Judge_Type(path):
    """
    返回不同文件的类型
    :param infile: 文件路径
    :return: 返回文件类型
    """
    suffix_dict = {'.gz': 0, '.fq': 1, '.fastq': 1,
                   'fa': 2, '.fas': 2, '.fasta': 2}
    return suffix_dict.get(os.path.splitext(path)[-1].lower(), 3)

def Make_Kmer_Dict(_kmer_dict, file_path, kmer_size):
    """
    制作整数型Kmer字典
    :param _kmer_dict: 用来保存kmer的字典变量
    :param file_path: 参考序列的文件名
    :param kmer_size: kmer的长度
    """
    MASK_BIN = (1 << (kmer_size << 1)) - 1  # kmer的掩码
    DEPTH_BIN = 1 << 10  # kmer深度的递增

    with open(file_path, 'r') as f:
        # 1-10位为位置千分比，11-30为深度，31-35为预留符号位，36及以后为文件
        file_id = 1 << 35  # 设置当前文件的符号位

        for _, seq in SimpleFastaParser(f):
            seq = ''.join(filter(str.isalpha, seq)).upper()
            total_kmer_count = len(seq) - kmer_size + 1
            if total_kmer_count <= 0:
                continue

            for run_start, run in Valid_Sequence_Runs_With_Positions(seq):
                intseqs, ref_len = Seq_To_Int(run)
                ref_kmer_count = ref_len - kmer_size + 1
                if ref_kmer_count <= 0:
                    continue

                for x, y in enumerate(intseqs):
                    # 初始化符号位和文件位，反向互补序列的31位符号位为1
                    SIGN_BIN = (1 << 30) + (1 << 10) + file_id if x else (1 << 10) + file_id

                    for j in range(0, ref_kmer_count):
                        temp_int = SIGN_BIN  # 初始化文件位和深度
                        kmer_int = y >> (j << 1) & MASK_BIN  # 获取kmer的整数形式

                        if kmer_int in _kmer_dict:
                            temp_int = _kmer_dict[kmer_int]
                            temp_int += DEPTH_BIN  # 深度加1，深度大于2**20会溢出,不太可能这么深的kmer
                            temp_int |= file_id  # 赋值文件位
                        else:
                            if x:
                                global_j = run_start + j
                            else:
                                global_j = len(seq) - run_start - ref_len + j
                            temp_int += int((global_j + 1) / total_kmer_count * 1000)

                        _kmer_dict[kmer_int] = temp_int

def Reference_File_Cache_Key(file_path, kmer_size):
    stat = os.stat(file_path)
    digest = hashlib.sha256()
    digest.update(os.path.abspath(file_path).encode())
    digest.update(b'\0')
    digest.update(str(stat.st_size).encode())
    digest.update(b'\0')
    digest.update(str(stat.st_mtime_ns).encode())
    return f'{os.path.basename(file_path)}.k{kmer_size}.{digest.hexdigest()[:16]}.pkl'

def Load_Or_Make_Kmer_Dict(_kmer_dict, file_path, kmer_size, cache_dir=None):
    if not cache_dir:
        Make_Kmer_Dict(_kmer_dict, file_path, kmer_size)
        return False

    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, Reference_File_Cache_Key(file_path, kmer_size))

    if os.path.isfile(cache_path):
        try:
            with open(cache_path, 'rb') as handle:
                _kmer_dict.update(pickle.load(handle))
            return True
        except (OSError, EOFError, pickle.PickleError, AttributeError, ValueError):
            _kmer_dict.clear()

    Make_Kmer_Dict(_kmer_dict, file_path, kmer_size)

    temp_path = f'{cache_path}.{os.getpid()}.tmp'
    with open(temp_path, 'wb') as handle:
        pickle.dump(_kmer_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temp_path, cache_path)
    return False

def Get_Ref_Info(ref_path, _ref_path_dict, _ref_count_dict):
    """
    制作保存参考序列路径和长度的字典
    :param ref_path: 参考序列路径
    :param _ref_path_dict: 参考序列路径字典
    :param _ref_count_dict: 参考序列条数字典

    """
    if os.path.isdir(ref_path):
        ref_path_list = [os.path.join(ref_path, ent.name) for ent in os.scandir(ref_path) if ent.is_file()]
    elif os.path.isfile(ref_path):
        ref_path_list = [ref_path]
    else:
        ref_path_list = []

    for file in ref_path_list:
        if Judge_Type(file) != 2:
            continue

        file_name = os.path.splitext(os.path.basename(file))[0]

        with open(file) as f:
            ref_seq_count = sum(1 for _ in SimpleFastaParser(f))

        _ref_count_dict[file_name] = ref_seq_count
        _ref_path_dict[file_name] = file

def Make_Assemble_Dict(file_list, kmer_size, _kmer_dict, _ref_dict, Filted_File_Ext = '.fq'):
    """
    构建拼接用的字典
    :param file_list: 文件列表
    :param kmer_size: kmer的长度
    :param _kmer_dict: 待生成的字典value的格式为[深度，位置（1000以内的整数）]
    :param _ref_dict: 参考序列的字典
    :return: 返回kmer的总数量
    """
    MASK_BIN = (1 << (kmer_size << 1)) - 1 # kmer的掩码
    fasta_file = Filted_File_Ext == '.fasta'
    for file in file_list:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            if fasta_file:
                temp_str = []
                while line and line[0] != '>':
                    temp_str.append(line)
                    line = infile.readline()
                read_seq = ''.join(filter(str.isalpha, ''.join(temp_str).upper()))
            else:
                read_seq = ''.join(filter(str.isalpha, line)).upper()
                infile.readline()
                infile.readline()
                infile.readline()
            kmer_set = set()
            for run in Valid_Sequence_Runs(read_seq):
                intseqs, read_len = Seq_To_Int(run)
                kmer_set.update(x >> (j << 1) & MASK_BIN
                                for x in intseqs
                                for j in range(0, read_len - kmer_size + 1))
            for kmer in kmer_set:
                if kmer in _kmer_dict:
                    _kmer_dict[kmer][0] += 1
                elif kmer in _ref_dict: # kmer的位置
                    temp_int = int(_ref_dict[kmer])
                    temp_depth = (temp_int >> 10) & ((1<<20) -1) #在参考序列中的深度
                    temp_pos = temp_int & 1023
                    is_reverse = bool(temp_int & 1073741824) # 判断是否为反向互补的序列
                    if is_reverse:
                        # 标记为反向的的kmer
                        temp_pos = 1000 - temp_pos
                    _kmer_dict[kmer] = [1, temp_pos, is_reverse, temp_depth]
                else:
                    _kmer_dict[kmer] = [1, 1023, 1, 0]
        infile.close()

def Get_Middle_Fragment(text, slice_len):
    """
    截取reads中间的高质量片段
    """
    start = (len(text) - slice_len) >> 1
    end = start + slice_len
    return text[start:end]

def Make_Reads_Dict(file_list, _reads_dict, Filted_File_Ext = '.fq'):
    """
    截取reads中间的片段，构建高质量的reads字典
    :param file_list: 文件列表
    :param _reads_dict: 待生成的字典value的格式为seq
    :return: 返回切片的长度
    """
    def read_sequences(path):
        with open(path, 'r', encoding='utf-8', errors='ignore') as infile:
            if Filted_File_Ext == '.fasta':
                for _, seq in SimpleFastaParser(infile):
                    yield ''.join(filter(str.isalpha, seq)).upper()
            else:
                while True:
                    header = infile.readline()
                    if not header:
                        break
                    sequence = infile.readline()
                    plus = infile.readline()
                    quality = infile.readline()
                    if not sequence or not plus or not quality:
                        raise ValueError(f'Truncated FASTQ record in {path}')
                    yield ''.join(filter(str.isalpha, sequence)).upper()

    min_length = None
    for file in file_list:
        for seq in read_sequences(file):
            if seq and (min_length is None or len(seq) < min_length):
                min_length = len(seq)

    if min_length is None:
        return 0

    # Use a common slice length that every read can support.  The former
    # first-read heuristic made shorter reads permanently invisible.
    slice_len = max(1, int(min_length * 0.9))

    for file in file_list:
        for read_seq in read_sequences(file):
            if len(read_seq) < slice_len:
                continue
            intseqs = [read_seq, Reverse_Complement_ACGT(read_seq)]
            # A palindromic middle fragment is identical in both orientations;
            # it still represents one physical read and must only be counted once.
            for slice_read in {Get_Middle_Fragment(seq, slice_len) for seq in intseqs}:
                _reads_dict[slice_read] = _reads_dict.get(slice_read, 0) + 1
    return slice_len

def Median(x):
    """
    使用中位数分割列表
    :return: 左半边列表，右半边列表，中位数
    """ 
    x = sorted(x)
    length = len(x)
    mid, rem = divmod(length, 2)    # divmod函数返回商和余数
    if rem:
        return x[:mid], x[mid+1:], x[mid]
    else:
        return x[:mid], x[mid:], (x[mid-1]+x[mid])/2

def Quartile(x):
    """
    获取列表的四分位数
    :return: 左四分位数，中位数，右四分位数，最大值加1
    """ 
    lHalf, rHalf, q2 = Median(x)
    return Median(lHalf)[2], q2, Median(rHalf)[2], max(x) + 1

def Median_Value(x):
    if not x:
        return 0

    return Median(x)[2]

def Get_Weight(_pos, new_pos, weight = 4):
    """
    距离和权重的关系模型，默认权重值为16，最高权重值为256，最低为0
    :param _pos: 当前kmer在参考序列的位置
    :param new_pos: 新的kmer位置
    :param weight: 默认权重
    :return: 返回计算后的权重
    """ 
    return int.bit_length((1024 - abs(_pos - new_pos)) >> 2) if (_pos and new_pos) else weight

def Get_Forward_Contig_v6(_dict, seed, kmer_size, iteration = 1024):
    """
    带权重的DBG贪婪拼接
    :param _pos: 当前kmer在参考序列的位置
    :param seed: 新的kmer位置
    :param kmer_size: kmer的大小
    :param iteration: 最大循环数量
    :param weight: 默认权重
    :return: best_seq, kmer_set, best_kmc, best_pos, best_snp
    """ 
    temp_list, kmer_set, stack_list, pos_list = [seed], set([seed]), [], []
    temp_dict = Counter(temp_list)
    cur_kmc, cur_seq, contigs = deque(), deque(), []
    _pos, node_distance, best_kmc_sum = 0, 0, 0
    MASK = (1 << ((kmer_size << 1) - 2)) - 1
    while iteration:
        next_kmer = (temp_list[-1] & MASK) << 2
        node = [(i, _dict[i][1], _dict[i][0] + _dict[i][3])
                for i in range(next_kmer, next_kmer + 4)
                if i in _dict and not temp_dict[i]]
        node.sort(key=itemgetter(2), reverse=True)
        if not node: 
            iteration -= 1
            cur_kmc_sum = sum(cur_kmc)
            contigs.append((cur_kmc.copy(), cur_seq.copy()))
            if cur_kmc_sum > best_kmc_sum:
                best_kmc_sum = cur_kmc_sum
            for _ in range(node_distance):
                temp_dict[temp_list.pop()] -= 1
                cur_kmc.pop()
                cur_seq.pop()
            if not stack_list:
                break
            node, node_distance, _pos = stack_list.pop()
        if len(node) >= 2:
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0: _pos = node[0][1]
        temp_list.append(node[0][0])
        temp_dict[node[0][0]] += 1
        kmer_set.add(node[0][0])
        pos_list.append(node[0][1])
        cur_kmc.append(node[0][2])
        cur_seq.append(node[0][0] & 3)
        node_distance += 1
    return contigs, kmer_set, pos_list, int(best_kmc_sum)


def Get_Forward_Uce_Backbone(_dict, seed, kmer_size, iteration=1024, lookahead=24):
    """Build one UCE extension without branch-stack backtracking.

    At a bubble, inspect a bounded greedy continuation for every outgoing edge,
    commit the first edge of the longest trace, and discard the alternatives.
    Visited k-mers are never revisited, so cycles terminate instead of causing
    repeated exploration.
    """
    mask = (1 << ((kmer_size << 1) - 2)) - 1
    lookahead = max(1, lookahead)
    path = [seed]
    path_set = {seed}
    discarded_kmers = set()
    kmer_set = {seed}
    pos_list = []
    weights = deque()
    bases = deque()

    def outgoing(current, blocked):
        next_kmer = (current & mask) << 2
        nodes = [
            (candidate, _dict[candidate][1],
             _dict[candidate][0] + _dict[candidate][3])
            for candidate in range(next_kmer, next_kmer + 4)
            if candidate in _dict and candidate not in blocked
            and candidate not in discarded_kmers
        ]
        nodes.sort(key=lambda node: (node[2], node[1] > 0), reverse=True)
        return nodes

    def trace(first_node):
        trace_nodes = []
        trace_seen = set(path_set)
        node = first_node
        for _ in range(lookahead):
            if node[0] in trace_seen:
                break
            trace_nodes.append(node)
            trace_seen.add(node[0])
            following = outgoing(node[0], trace_seen)
            if not following:
                break
            # Look-ahead is deliberately linear: no nested fork or backtrack.
            node = following[0]
        return trace_nodes

    for _ in range(max(0, iteration)):
        nodes = outgoing(path[-1], path_set)
        if not nodes:
            break
        if len(nodes) == 1:
            chosen = nodes[0]
        else:
            traces = [trace(node) for node in nodes]
            winning_trace = max(
                traces,
                key=lambda candidate: (
                    len(candidate),
                    sum(node[2] for node in candidate),
                    candidate[0][2],
                ),
            )
            chosen = winning_trace[0]
            discarded_kmers.update(
                node[0] for node in nodes if node[0] != chosen[0]
            )

        kmer, ref_pos, weight = chosen
        path.append(kmer)
        path_set.add(kmer)
        kmer_set.add(kmer)
        pos_list.append(ref_pos)
        weights.append(weight)
        bases.append(kmer & 3)
    return [(weights, bases)], kmer_set, pos_list, int(sum(weights))


def find_position(dq, n):
    for i in range(len(dq) - 1, -1, -1):
        if dq[i] >= n:
            return i
    return -1

def Locate_Read_Slices(seq, slice_len, reads_dict):
    """Return matching read slices and their unique placement, if any.

    Each dictionary key represents a collection of physical reads with the same
    middle slice.  A key contributes to read support once per contig, regardless
    of how many times the sequence occurs in a repeat.  Repeated placements are
    recorded as ``None`` because they cannot define a supported coordinate.
    """
    if slice_len <= 0 or len(seq) < slice_len:
        return {}

    matches = {}
    for pos in range(len(seq) - slice_len + 1):
        slice_str = seq[pos:pos + slice_len]
        if slice_str not in reads_dict:
            continue
        if slice_str in matches:
            matches[slice_str] = None
        else:
            matches[slice_str] = pos
    return matches

def Process_Contigs(contigs, max_weight, slice_len, reads_dict, soft_boundary = 0, assembly_mode = 'reference'):
    """
    通过将contigs与reads进行map，来检测contig的可靠性
    :param contigs: 拼接过程获取的contigs
    :param max_weight: 最大的权重，只考虑大于最大权重一半的contigs
    :param slice_len: reads的高质量切片的长度
    :param reads_dict: reads的高质量切片的词典
    :return: 按照map上的reads的数量倒序排序过后的contigs
    """ 
    if assembly_mode == 'reference':
        # 基于soft_boundary和四分位点切割序列两端
        for i, contig in enumerate(contigs):
            if len(contig[0]) > 2:
                cut_value = Quartile(contig[0])[0]
                cut_pos = find_position(contig[0], cut_value)
                if cut_pos != -1 and cut_pos + soft_boundary + 1 < len(contig[0]):
                    while len(contig[0]) > cut_pos + soft_boundary + 1:
                        contig[0].pop()
                        contig[1].pop()

    min_weight = max_weight >> (2 if assembly_mode == 'uce' else 1)
    processed_contigs = sorted([[''.join(ACGT_DICT[k] for k in x[1]), sum(x[0]), 0] for x in contigs if sum(x[0]) > min_weight], key=itemgetter(1), reverse=True)
    for x in processed_contigs:
        matches = Locate_Read_Slices(x[0], slice_len, reads_dict)
        x[2] = sum(reads_dict[slice_str] for slice_str in matches)
    if assembly_mode == 'uce':
        processed_contigs.sort(key=lambda x: (len(x[0]), x[2], x[1]), reverse=True)
    else:
        processed_contigs.sort(key=itemgetter(2), reverse=True)
    return processed_contigs


def Calculate_Read_Support(seq, slice_len, reads_dict):
    """
    Measure direct read-slice support for a candidate contig.

    ``supported_extent`` preserves the historical leftmost-to-rightmost
    statistic.  ``supported_bases`` is the union of uniquely placed slice
    intervals and is therefore the quantity used for UCE scoring and quality
    control.  Repetitively placed slices contribute to the total recruited
    read count, but not to positional support.
    """
    contig_len = len(seq)
    matches = Locate_Read_Slices(seq, slice_len, reads_dict)
    total_read_count = sum(reads_dict[slice_str] for slice_str in matches)
    unique_read_count = 0
    multi_mapping_read_count = 0
    intervals = []

    for slice_str, pos in matches.items():
        count = reads_dict[slice_str]
        if pos is None:
            multi_mapping_read_count += count
            continue

        unique_read_count += count
        intervals.append((pos, min(pos + slice_len, contig_len)))

    if not intervals:
        return ReadSupport(
            total_read_count, 0, multi_mapping_read_count, 0, 0, 0.0,
            contig_len, 0.0, contig_len, 0)

    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    left_coord = merged[0][0]
    right_coord = merged[-1][1]
    supported_extent = right_coord - left_coord
    supported_bases = sum(end - start for start, end in merged)
    breadth = supported_bases / contig_len if contig_len else 0.0

    gaps = [left_coord, contig_len - right_coord]
    gaps.extend(merged[i + 1][0] - merged[i][1] for i in range(len(merged) - 1))
    max_gap = max(gaps, default=0)

    if supported_extent:
        left_extension = left_coord
        right_extension = contig_len - right_coord
        if left_extension == 0 and right_extension == 0:
            flank_balance = 1.0
        else:
            flank_balance = min(left_extension, right_extension) / max(left_extension, right_extension)
    else:
        flank_balance = 0.0

    return ReadSupport(
        total_read_count, unique_read_count, multi_mapping_read_count,
        supported_extent, supported_bases, breadth, max_gap, flank_balance,
        left_coord, right_coord)


def Calculate_Kmer_Depth_Stats(seq, kmer_size, assemble_dict):
    seq_str = seq.translate(FWD_TRANS)

    if len(seq_str) < kmer_size:
        return 0, 0, 0

    mask_bin = (1 << (kmer_size << 1)) - 1
    seq_int = int(seq_str, 4)
    counts = []

    for _ in range(0, len(seq_str) - kmer_size + 1):
        kmer = seq_int & mask_bin
        value = assemble_dict.get(kmer)
        counts.append(value[0] if value else 0)
        seq_int >>= 2

    median_depth = Median_Value(counts)

    if not counts or median_depth <= 0:
        return median_depth, 0, 0

    mean_depth = sum(counts) / len(counts)
    if mean_depth <= 0:
        depth_cv = 0
    else:
        variance = sum((x - mean_depth) ** 2 for x in counts) / len(counts)
        depth_cv = math.sqrt(variance) / mean_depth

    max_depth_ratio = max(counts) / median_depth
    return median_depth, depth_cv, max_depth_ratio


def Score_Contig(contig, assembly_mode):
    """
    Reference mode keeps the historical conservative ranking.
    UCE mode prefers candidates that extend farther while retaining read support.
    """
    seq_len = len(contig[CONTIG_SEQ])
    weight = contig[CONTIG_WEIGHT]
    read_count = contig[CONTIG_READ_COUNT]
    supported_span = contig[CONTIG_SUPPORTED_SPAN]
    flank_balance = contig[CONTIG_FLANK_BALANCE]

    if assembly_mode == 'uce':
        unique_read_count = contig[CONTIG_UNIQUE_READ_COUNT]
        supported_bases = contig[CONTIG_SUPPORTED_BASES]
        support_breadth = contig[CONTIG_SUPPORT_BREADTH]
        max_support_gap = contig[CONTIG_MAX_SUPPORT_GAP]
        unique_read_density = unique_read_count / seq_len if seq_len else 0
        depth_cv = contig[CONTIG_KMER_DEPTH_CV]
        max_depth_ratio = contig[CONTIG_KMER_MAX_DEPTH_RATIO]
        density_factor = min(unique_read_density / 0.01, 1.0)
        continuity_factor = 1 / (1 + depth_cv)
        repeat_factor = min(10 / max(max_depth_ratio, 1), 1.0)
        effective_supported_bases = supported_bases * density_factor * continuity_factor * repeat_factor
        gap_fraction = max_support_gap / seq_len if seq_len else 1.0
        return (effective_supported_bases, support_breadth, unique_read_density,
                -gap_fraction, seq_len, unique_read_count, read_count,
                flank_balance, weight)

    return (read_count, weight)


def Build_Uce_Guardrails(args):
    return {
        'max_contig_length': args.uce_max_contig_length,
        'min_read_density': args.uce_min_read_density,
        'density_check_min_length': args.uce_density_check_min_length,
        'max_depth_cv': args.uce_max_depth_cv,
        'max_depth_ratio': args.uce_max_depth_ratio,
    }


def Get_Uce_Guardrail_Reasons(contig_len, read_density, depth_cv, max_depth_ratio, guardrails):
    if not guardrails:
        return []

    reasons = []

    if guardrails['max_contig_length'] > 0 and contig_len > guardrails['max_contig_length']:
        reasons.append('contig_too_long')

    if contig_len >= guardrails['density_check_min_length'] and read_density < guardrails['min_read_density']:
        reasons.append('low_unique_read_density')

    if guardrails['max_depth_cv'] > 0 and depth_cv > guardrails['max_depth_cv']:
        reasons.append('high_depth_cv')

    if guardrails['max_depth_ratio'] > 0 and max_depth_ratio > guardrails['max_depth_ratio']:
        reasons.append('repeat_depth_peak')

    return reasons


def Pass_Uce_Guardrails(contig_len, read_density, depth_cv, max_depth_ratio, guardrails):
    return not Get_Uce_Guardrail_Reasons(
        contig_len, read_density, depth_cv, max_depth_ratio, guardrails)


def Evaluate_Uce_Candidate(contig_len, unique_read_count, supported_bases,
                           unique_read_density, depth_cv, max_depth_ratio,
                           guardrails):
    """Return an auditable UCE acceptance decision and rejection reasons."""
    reasons = []

    if unique_read_count <= 0:
        reasons.append('no_unique_read_support')
    if supported_bases <= 0:
        reasons.append('no_positional_support')

    reasons.extend(Get_Uce_Guardrail_Reasons(
        contig_len, unique_read_density, depth_cv, max_depth_ratio, guardrails))

    return not reasons, reasons


def Build_Contig_Record(raw_contig, seed_count, contig_pos):
    return [
        raw_contig[RAW_SEQ],
        seed_count,
        contig_pos,
        raw_contig[RAW_WEIGHT],
        raw_contig[RAW_READ_COUNT],
        raw_contig[RAW_SUPPORTED_SPAN],
        raw_contig[RAW_FLANK_BALANCE],
        raw_contig[RAW_READ_DENSITY],
        raw_contig[RAW_SUPPORT_FRACTION],
        raw_contig[RAW_KMER_MEDIAN_DEPTH],
        raw_contig[RAW_KMER_DEPTH_CV],
        raw_contig[RAW_KMER_MAX_DEPTH_RATIO],
        raw_contig[RAW_UNIQUE_READ_COUNT],
        raw_contig[RAW_MULTI_MAPPING_READ_COUNT],
        raw_contig[RAW_SUPPORTED_BASES],
        raw_contig[RAW_SUPPORT_BREADTH],
        raw_contig[RAW_MAX_SUPPORT_GAP],
        raw_contig[RAW_ACCEPTED],
        raw_contig[RAW_REJECTION_REASON],
    ]


def Format_Contig_Header(contig, assembly_mode, prefix='contig'):
    header = (
        f'>{prefix}_{len(contig[CONTIG_SEQ])}_{contig[CONTIG_SEED_COUNT]}'
        f'_{contig[CONTIG_POS]}_{contig[CONTIG_WEIGHT]}'
        f'_{contig[CONTIG_READ_COUNT]}_span_{contig[CONTIG_SUPPORTED_SPAN]}'
    )
    if assembly_mode == 'uce':
        header += f'_supported_{contig[CONTIG_SUPPORTED_BASES]}'
    return header + f'_balance_{contig[CONTIG_FLANK_BALANCE]:.3f}'


def Get_Contig_v6(_reads_dict, slice_len, _dict, seed, kmer_size, cov_min, iteration = 1024, soft_boundary = 0, assembly_mode = 'reference', uce_side_candidates = 8, uce_guardrails = None, uce_path_strategy='search', uce_backbone_lookahead=24):
    """
    获取最优的contig
    :param _reads_dict: reads的高质量切片的词典
    :param slice_len: reads的高质量切片的长度
    :param _dict: 用于拼接的kmer字典
    :param seed: 拼接种子
    :param kmer_size: kmer的长度
    :param iteration: 构建contig时允许的最大路径分支数
    :param weight: 没有ref时的默认权重
    :return: contigs的集合，用到所有的kmer的集合，contig的大概位置
    """ 
    if assembly_mode == 'uce' and uce_path_strategy == 'backbone':
        contigs_1, kmer_set_1, pos_list_1, weight_1 = Get_Forward_Uce_Backbone(
            _dict, seed, kmer_size, iteration, uce_backbone_lookahead)
        contigs_2, kmer_set_2, pos_list_2, weight_2 = Get_Forward_Uce_Backbone(
            _dict, Reverse_Int(seed, kmer_size), kmer_size, iteration,
            uce_backbone_lookahead)
    else:
        contigs_1, kmer_set_1, pos_list_1, weight_1 = Get_Forward_Contig_v6(
            _dict, seed, kmer_size, iteration)
        contigs_2, kmer_set_2, pos_list_2, weight_2 = Get_Forward_Contig_v6(
            _dict, Reverse_Int(seed, kmer_size), kmer_size, iteration)
    # 清理位置列表
    pos_list = [x for x in chain(pos_list_1, pos_list_2) if x > 0 and x < 1000]
    # 获取位置中位数
    contig_pos = int(Quartile(pos_list)[1] if len(pos_list)>1 else -1)
    # 获取最可能的两侧的contig
    contigs_1_16 = Process_Contigs(contigs_1, weight_1, slice_len, _reads_dict, soft_boundary, assembly_mode)
    contigs_2_16 = Process_Contigs(contigs_2, weight_2, slice_len, _reads_dict, soft_boundary, assembly_mode)
    processed_contigs = []
    if not contigs_1_16: contigs_1_16.append(['',0,0])
    if not contigs_2_16: contigs_2_16.append(['',0,0])
    candidate_limit = (
        1 if assembly_mode == 'uce' and uce_path_strategy == 'backbone'
        else uce_side_candidates if assembly_mode == 'uce' else 3
    )
    for l in contigs_2_16[:candidate_limit]:
        for r in contigs_1_16[:candidate_limit]:
            c = Reverse_Complement_ACGT(l[0]) + Int_To_Seq(seed, kmer_size) + r[0]
            c_weight = l[1] + r[1]
            contig_len = len(c)
            support = Calculate_Read_Support(c, slice_len, _reads_dict)
            read_density = support.total_read_count / contig_len if contig_len else 0
            support_fraction = support.supported_extent / contig_len if contig_len else 0
            median_depth, depth_cv, max_depth_ratio = (0, 0, 0)
            if assembly_mode == 'uce':
                median_depth, depth_cv, max_depth_ratio = Calculate_Kmer_Depth_Stats(c, kmer_size, _dict)
            positional_count = (support.unique_read_count
                                if assembly_mode == 'uce' else support.total_read_count)
            positional_span = (support.supported_bases
                               if assembly_mode == 'uce' else support.supported_extent)
            cov_dep = positional_count * slice_len / 0.9
            if cov_min > 0:
                if positional_span == 0 or cov_dep / positional_span < cov_min:
                    continue
                if cov_dep / contig_len < cov_min:
                    c = c[support.left_coord:support.right_coord]
                    contig_len = len(c)
                    support = Calculate_Read_Support(c, slice_len, _reads_dict)
                    read_density = support.total_read_count / contig_len if contig_len else 0
                    support_fraction = support.supported_extent / contig_len if contig_len else 0
                    if assembly_mode == 'uce':
                        median_depth, depth_cv, max_depth_ratio = Calculate_Kmer_Depth_Stats(c, kmer_size, _dict)
            accepted = True
            rejection_reasons = []
            if assembly_mode == 'uce':
                unique_read_density = support.unique_read_count / contig_len if contig_len else 0
                accepted, rejection_reasons = Evaluate_Uce_Candidate(
                    contig_len, support.unique_read_count, support.supported_bases,
                    unique_read_density, depth_cv, max_depth_ratio, uce_guardrails)
            # 序列，序列的拼接权重，切片数
            processed_contigs.append([
                c, c_weight, support.total_read_count, support.supported_extent,
                support.flank_balance, read_density, support_fraction,
                median_depth, depth_cv, max_depth_ratio,
                support.unique_read_count, support.multi_mapping_read_count,
                support.supported_bases, support.breadth, support.max_gap,
                accepted, ';'.join(rejection_reasons),
            ])
    return processed_contigs, kmer_set_1 | kmer_set_2, contig_pos

def Calculate_Kmer_Size(ref_path, reads, slice_len, k_min, k_max, error_limit):
    if k_min % 2 == 0:
        k_min += 1

    if slice_len <= k_min:
        return k_min

    mask_bin  = (1 << (k_min << 1)) - 1
    kmer_dict = Counter()
    trans = FWD_TRANS
    rtrans = REV_TRANS

    for seq in reads:
        for run in Valid_Sequence_Runs(seq):
            seq_str   = run.translate(trans)
            seq_str_r = run.translate(rtrans)[::-1]

            if len(seq_str) < k_min:
                continue

            seq_int   = int(seq_str, 4)
            seq_int_r = int(seq_str_r, 4)

            for _ in range(0, len(seq_str) - k_min + 1):
                kmer_dict[seq_int   & mask_bin] += 1
                kmer_dict[seq_int_r & mask_bin] += 1
                seq_int   >>= 2
                seq_int_r >>= 2

    kmer_dict = {k for k, v in kmer_dict.items() if v > error_limit}
    run_length_stats = [0] * (k_max - k_min + 1)
    run_maximum = k_max - k_min + 1

    with open(ref_path, 'r') as f:
        for _, seq in SimpleFastaParser(f):
            seq       = ''.join(filter(str.isalpha, seq)).upper()
            for run in Valid_Sequence_Runs(seq):
                seq_str = run.translate(trans)

                if len(seq_str) < k_min:
                    continue

                seq_int = int(seq_str, 4)

                run_length_list = [0]

                for _ in range(0, len(seq_str) - k_min + 1):
                    if (seq_int & mask_bin) in kmer_dict:
                        run_length_list[-1] += 1
                        if run_length_list[-1] >= run_maximum:
                            run_length_list.append(run_maximum // 2)
                    elif run_length_list[-1] != 0:
                        run_length_list.append(0)
                    seq_int >>= 2

                for k, v in Counter(run_length_list).items():
                    if k == 0:
                        continue

                    kp = k - 1
                    odd = kp % 2
                    kp = kp - odd

                    run_length_stats[kp] += v

                    for i in range(2, kp + 1, 2):
                        run_length_stats[kp - i] += v

    for k, n in reversed(tuple(enumerate(run_length_stats, k_min))):
        if n > 0:
            upper_bound = k
            break
    else:
        return k_min

    lower_bound = (upper_bound + 1) // 2

    k_list = [(k, n * k / (slice_len - k + 1) if slice_len > k else 0)
              for k, n in enumerate(run_length_stats, k_min)
              if lower_bound < k <= upper_bound]
    cutoff = max(map(itemgetter(1), k_list)) / 2

    for k, cov in reversed(k_list):
        if cov > cutoff:
            return k
    else:
        return k_min

def Write_Dict(_dict, file_name):
    """
    写入16进制kmer字典
    :param _dict: kmer字典
    :param file_name: 文件名，
    :return: 成功返回1
    """
    with open(file_name,'w') as f:
        for key, value in _dict.items():
            if type(value) == list:
                f.writelines([str(key), ",", ",".join(map(str,value)), ',', '\n'])
            else:
                f.writelines([str(key), ",", str(value), ",", '\n'])

def Read_Dict(file_name):
    if not os.path.isfile(file_name):
        return {}

    result = {}
    with open(file_name, newline='') as f:
        for row in csv.reader(f):
            if len(row) >= 3 and row[0]:
                result[row[0]] = row[1:3]
    return result

def Read_Uce_Summary(file_name):
    if not os.path.isfile(file_name):
        return {}

    with open(file_name, newline='') as f:
        return {row['locus']: row for row in csv.DictReader(f) if row.get('locus')}

def Write_Uce_Summary(rows, file_name):
    fieldnames = [
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

    with open(file_name, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, '') for name in fieldnames})

def process_key_value(args, key, ref_path, ref_count, iteration, soft_boundary, loop_count, total_count, completed_keys=None):
    contig_best_path = os.path.join(args.o, "results", key + ".fasta")
    contig_all_path = os.path.join(args.o, "contigs_all", key + ".fasta")
    contig_low_path = os.path.join(args.o, "contigs_all_low", key + ".fasta")
    current_ka = args.ka
    limit = args.limit_count

    completed_keys = set() if completed_keys is None else completed_keys
    if key in completed_keys and os.path.isfile(contig_best_path) and os.path.getsize(contig_best_path) > 0:
        return False, key, {"status": "skipped"}

    # 检查是哪种扩展名
    file_extensions = ['.fasta', '.fq']
    Filted_File_Ext = '.fq'
    filtered_file_path = None
    for ext in file_extensions:
        file_path = os.path.join(args.o, 'filtered', key + ext)
        if os.path.exists(file_path):
            filtered_file_path = file_path
            Filted_File_Ext = ext
            break

    # 清理文件
    if not os.path.isfile(os.path.join(args.o, 'filtered', key + Filted_File_Ext)):
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)
        if os.path.isfile(contig_low_path): os.remove(contig_low_path)
        return False, key, {"status": "no filtered file", "value": 0}

    # 获取种子列表
    ref_dict, filtered_dict, reads_dict = {}, {}, {}

    # 获取最大切片长度，建立reads切片字典
    slice_len = Make_Reads_Dict([filtered_file_path], reads_dict)

    if not reads_dict:
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)
        if os.path.isfile(contig_low_path): os.remove(contig_low_path)
        Write_Print(os.path.join(args.o,  "log.txt"), "No reads were obtained for gene", key)
        return False, key, {"status": "no reads", "value": 0}

    # 自动调整soft_boundary
    if soft_boundary == -1:
        soft_boundary = slice_len // 2

    # 如果不指定ka, 估算最大ka，执行动态高精度拼接
    if not current_ka:
        current_ka = Calculate_Kmer_Size(ref_path, reads_dict.keys(), slice_len, args.k_min, args.k_max, limit)

    Write_Print(os.path.join(args.o,  "log.txt"), "Use k=", current_ka, " for assembling gene ", key ,".", sep='')
    Write_Print(os.path.join(args.o,  "log.txt"), 'Assembling', key, loop_count, '/', total_count)

    # 制作参考序列的kmer字典
    Load_Or_Make_Kmer_Dict(ref_dict, ref_path, current_ka, args.assembler_reference_cache_dir)
    # 制作用于拼接的kmer字典
    Make_Assemble_Dict([filtered_file_path], current_ka, filtered_dict, ref_dict)
    # 缩减filtered_dict，保留大于limit和有深度信息的
    if limit > 0:
        filtered_dict = {k: v for k, v in filtered_dict.items() if v[0] > limit or v[3] > 0}

    if len(filtered_dict) < 3:
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)
        if os.path.isfile(contig_low_path): os.remove(contig_low_path)
        Write_Print(os.path.join(args.o,  "log.txt"), 'Could not get enough reads from filter.')
        return False, key, {"status": "insufficient genomic kmers", "value": 0}

    # 纠正深度上限, 获取参考序列的深度修正权重
    # filtered_dict[0排除了上限的过滤深度，1位置，2方向，3修正参考序列深度]
    read_quar = Quartile([v[0] for v in filtered_dict.values()])
    depth_upper = int((read_quar[2] - read_quar[0]) * 1.5 + read_quar[2])
    for k, v in filtered_dict.items():
        if v[3] != 0:
            v[3] = (v[0] > limit) * int(ref_count / (abs(v[3] - ref_count) + 1) * depth_upper) + 1
        v[0] = min(v[0], depth_upper)

    # 在每个参考序列中出现且只出现一次的kmer优先作为种子
    # 长度位置在1~1000之间，与参考序列方向一致v[2] == 0
    seed_list = [(k, v[0], v[1], v[3]) for k, v in filtered_dict.items() if v[1] > 1 and v[1] < 1000 and not v[2]]
    seed_list.sort(key=lambda x: (x[3], x[1]), reverse=True)

    # 必须有seed_list, 否则意味着跟参考序列差别过大
    if not seed_list:
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)
        if os.path.isfile(contig_low_path): os.remove(contig_low_path)
        Write_Print(os.path.join(args.o,  "log.txt"), 'Could not get enough seeds.')
        return False, key, {"status": "no seed", "value": 0}

    # 获取seed集合，用来加速集合操作
    seed_list_len = len(seed_list)
    seed_set = set(i[0] for i in seed_list)

    # 获取contigs
    contigs_all = []
    contigs_all_low = []
    contigs_best = []
    uce_guardrails = Build_Uce_Guardrails(args) if args.assembly_mode == 'uce' else None

    while len(seed_list) > seed_list_len * 0.5: # 已经耗费了大于一半的seed就没必要再做了 
        # org_contigs: 0序列 1序列的拼接权重 2切片数
        org_contigs, kmer_set, contig_pos = Get_Contig_v6(
            reads_dict, slice_len, filtered_dict, seed_list[0][0], current_ka, args.cov_min,
            iteration=iteration, soft_boundary=soft_boundary, assembly_mode=args.assembly_mode,
            uce_side_candidates=args.uce_side_candidates, uce_guardrails=uce_guardrails,
            uce_path_strategy=args.uce_path_strategy,
            uce_backbone_lookahead=args.uce_backbone_lookahead)
        seed_list = [item for item in seed_list if (item[0] not in kmer_set) and (Reverse_Int(item[0], current_ka) not in kmer_set)]
        for contig in org_contigs:
            contig_record = Build_Contig_Record(contig, len(seed_set & kmer_set), contig_pos)
            if args.assembly_mode == 'uce':
                if contig[RAW_ACCEPTED]:
                    contigs_all.append(contig_record)
                else:
                    contigs_all_low.append(contig_record)
            elif contig[RAW_READ_COUNT] * slice_len > len(contig[RAW_SEQ]): # 起码要有reads高质量切片能够覆盖contig，否则就是错误的拼接
                contigs_all.append(contig_record)
            else:
                contigs_all_low.append(contig_record)

    low_qual = not contigs_all
    if low_qual and args.assembly_mode != 'uce':
        contigs_all = contigs_all_low

    # Reference mode preserves the historical fallback.  UCE low-quality
    # candidates remain diagnostic-only and never enter the primary results.
    selection_pool = contigs_all if contigs_all else contigs_all_low
    if selection_pool:
        selection_pool.sort(key=lambda x: Score_Contig(x, args.assembly_mode), reverse=True)
        contigs_best.append(selection_pool[0])
    else:
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)
        if os.path.isfile(contig_low_path): os.remove(contig_low_path)
        Write_Print(os.path.join(args.o, "log.txt"), "Insufficient reads coverage, unable to build contigs.")
        return False, key, {"status": "no contigs", "value": 0}

    write_primary = args.assembly_mode != 'uce' or not low_qual
    if write_primary:
        with open(contig_best_path, 'w') as out:
            for x in contigs_best:
                out.write(Format_Contig_Header(x, args.assembly_mode) + '\n')
                out.write(x[CONTIG_SEQ] + '\n')
        with open(contig_all_path, 'w') as out:
            for x in contigs_all:
                out.write(Format_Contig_Header(x, args.assembly_mode) + '\n')
                out.write(x[CONTIG_SEQ] + '\n')
    else:
        if os.path.isfile(contig_best_path): os.remove(contig_best_path)
        if os.path.isfile(contig_all_path): os.remove(contig_all_path)

    if args.assembly_mode == 'uce' and contigs_all_low:
        contigs_all_low.sort(key=lambda x: Score_Contig(x, args.assembly_mode), reverse=True)
        with open(contig_low_path, 'w') as out:
            for x in contigs_all_low:
                out.write(Format_Contig_Header(x, args.assembly_mode, 'low_support_contig') + '\n')
                out.write(x[CONTIG_SEQ] + '\n')
    elif os.path.isfile(contig_low_path):
        os.remove(contig_low_path)

    best_contig = contigs_best[0]
    ref_dict, filtered_dict = None, None
    gc.collect()
    accepted_locus = args.assembly_mode != 'uce' or not low_qual
    return accepted_locus, key, {
        "status": "low quality" if low_qual else "success",
        "value": best_contig[CONTIG_READ_COUNT],
        "accepted": int(not low_qual),
        "rejection_reason": best_contig[CONTIG_REJECTION_REASON] if low_qual else "",
        "selected_contig_length": len(best_contig[CONTIG_SEQ]),
        "read_supported_span": best_contig[CONTIG_SUPPORTED_SPAN],
        "slice_supported_bases": best_contig[CONTIG_SUPPORTED_BASES],
        "slice_support_breadth": round(best_contig[CONTIG_SUPPORT_BREADTH], 6),
        "max_slice_support_gap": best_contig[CONTIG_MAX_SUPPORT_GAP],
        "read_count": best_contig[CONTIG_READ_COUNT],
        "unique_read_count": best_contig[CONTIG_UNIQUE_READ_COUNT],
        "multi_mapping_read_count": best_contig[CONTIG_MULTI_MAPPING_READ_COUNT],
        "read_density": round(best_contig[CONTIG_READ_DENSITY], 6),
        "unique_read_density": round(best_contig[CONTIG_UNIQUE_READ_COUNT] / len(best_contig[CONTIG_SEQ]), 6),
        "support_fraction": round(best_contig[CONTIG_SUPPORT_FRACTION], 3),
        "flank_balance": round(best_contig[CONTIG_FLANK_BALANCE], 3),
        "kmer_median_depth": round(best_contig[CONTIG_KMER_MEDIAN_DEPTH], 3),
        "kmer_depth_cv": round(best_contig[CONTIG_KMER_DEPTH_CV], 3),
        "kmer_max_depth_ratio": round(best_contig[CONTIG_KMER_MAX_DEPTH_RATIO], 3),
        "candidate_count": len(contigs_all) + len(contigs_all_low),
        "low_quality": int(low_qual),
    }

if __name__ == '__main__':
    if sys.platform.startswith('win'):
        multiprocessing.freeze_support()

    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Assembler by YY 20230314''')
    pars.add_argument('-r', metavar='<str>', type=str, help='''input ref file or dir''', required=True)
    pars.add_argument('-o', metavar='<str>', type=str, help='''out dir''', required=True)
    pars.add_argument('-ka', metavar='<int>', type=int, help='''kmer of assemble''',  default=39)
    pars.add_argument('-k_max', metavar='<int>', type=int, help='''max kmer of assemble''',  default=39)
    pars.add_argument('-k_min', metavar='<int>', type=int, help='''max kmer of assemble''',  default=21)
    pars.add_argument('-limit_count', metavar='<int>', type=int, help='''limit of kmer count''', required=False, default=2)
    pars.add_argument('-iteration', metavar='<int>', type=int, help='''iteration''', required=False, default=8192)
    pars.add_argument('-cov_min', metavar='<int>', type=int, help='''min coverage''', required=False, default=0)
    pars.add_argument('-sb', '--soft_boundary', metavar='<int>', type=int, help='''soft boundary，default = [0], -1时为切片长度的一半''', required=False, default=0)
    pars.add_argument('-p', '--processes', metavar='<int>', type=int, help='Number of processes for multiprocessing', default= 1)#max(multiprocessing.cpu_count()-1,2))
    pars.add_argument('--assembly-mode', choices=('reference', 'uce'), default='reference', help='Assembly mode')
    pars.add_argument('--uce-side-candidates', dest='uce_side_candidates', metavar='<int>', type=int, default=8, help='''number of one-sided branch candidates to combine in UCE mode''')
    pars.add_argument('--uce-path-strategy', choices=('search', 'backbone'), default='backbone', help='UCE path strategy: backbone commits one bounded-lookahead path without backtracking; search preserves legacy branch enumeration')
    pars.add_argument('--uce-backbone-lookahead', dest='uce_backbone_lookahead', metavar='<int>', type=int, default=24, help='bounded greedy look-ahead per UCE backbone branch (default = 24)')
    pars.add_argument('--uce-max-contig-length', dest='uce_max_contig_length', metavar='<int>', type=int, default=0, help='''maximum UCE contig length kept before scoring; 0 disables''')
    pars.add_argument('--uce-min-read-density', dest='uce_min_read_density', metavar='<float>', type=float, default=0.003, help='''minimum uniquely placed read_count/length for long UCE contigs''')
    pars.add_argument('--uce-density-check-min-length', dest='uce_density_check_min_length', metavar='<int>', type=int, default=1000, help='''minimum UCE contig length where read-density guardrail applies''')
    pars.add_argument('--uce-max-depth-cv', dest='uce_max_depth_cv', metavar='<float>', type=float, default=0, help='''optional maximum kmer depth coefficient of variation for UCE contigs; 0 disables''')
    pars.add_argument('--uce-max-depth-ratio', dest='uce_max_depth_ratio', metavar='<float>', type=float, default=0, help='''optional maximum kmer max/median depth ratio for UCE contigs; 0 disables''')
    pars.add_argument('--assembler-reference-cache-dir', dest='assembler_reference_cache_dir', metavar='<str>', default=None, help='''optional directory for cached assembler reference k-mer dictionaries''')
    args = pars.parse_args()
    args.uce_side_candidates = max(args.uce_side_candidates, 3)
    args.uce_backbone_lookahead = max(args.uce_backbone_lookahead, 1)
    args.uce_max_contig_length = max(args.uce_max_contig_length, 0)
    args.uce_density_check_min_length = max(args.uce_density_check_min_length, 1)

    if args.uce_min_read_density < 0:
        pars.error('--uce-min-read-density must be greater than or equal to 0')
    if args.uce_max_depth_cv < 0:
        pars.error('--uce-max-depth-cv must be greater than or equal to 0')
    if args.uce_max_depth_ratio < 0:
        pars.error('--uce-max-depth-ratio must be greater than or equal to 0')

    try:
        # 初始化文件夹
        if not os.path.isdir(os.path.join(args.o, 'results')):
            os.mkdir(os.path.join(args.o, 'results'))
        if not os.path.isdir(os.path.join(args.o, 'contigs_all')):
            os.mkdir(os.path.join(args.o, 'contigs_all'))
        if not os.path.isdir(os.path.join(args.o, 'contigs_all_low')):
            os.mkdir(os.path.join(args.o, 'contigs_all_low'))
        # 载入参考序列信息
        Get_Ref_Info(args.r, ref_path_dict, ref_count_dict)
        t0 = time.time()
    except Exception as e:
        Write_Print(os.path.join(args.o,  "log.txt"), "error:" , e)
        sys.exit(1)

    try:
        Write_Print(os.path.join(args.o,  "log.txt"), '======================== Assemble =========================')

        result_path = os.path.join(args.o, "result_dict.txt")
        uce_summary_path = os.path.join(args.o, "uce_assembly_summary.csv")
        valid_keys = set(ref_path_dict)
        result_dict = {key: value for key, value in Read_Dict(result_path).items() if key in valid_keys}
        uce_summary_by_locus = ({key: value
                                 for key, value in Read_Uce_Summary(uce_summary_path).items()
                                 if key in valid_keys}
                                if args.assembly_mode == 'uce' else {})
        completed_keys = set(result_dict)
        if args.assembly_mode == 'uce':
            completed_keys &= set(uce_summary_by_locus)

        results = []
        if args.processes > 1:
            pool = multiprocessing.Pool(args.processes)
            for loop_count, (key, ref_path) in enumerate(ref_path_dict.items(), start=1):
                results.append(pool.apply_async(process_key_value, (args, key, ref_path, ref_count_dict[key], args.iteration, args.soft_boundary, loop_count, len(ref_path_dict), completed_keys)))
            pool.close()
            pool.join()
        else:
            for loop_count, (key, ref_path) in enumerate(ref_path_dict.items(), start=1):
                results.append(process_key_value(args, key, ref_path, ref_count_dict[key], args.iteration, args.soft_boundary, loop_count, len(ref_path_dict), completed_keys))

        for result in results:
            success, key_update, result_dict_entry = result if type(result) == tuple else result.get()
            if result_dict_entry.get("status") != "skipped":
                result_dict[key_update] = [result_dict_entry["status"], result_dict_entry["value"]]
                if args.assembly_mode == 'uce':
                    uce_summary_by_locus[key_update] = {'locus': key_update, **result_dict_entry}

        Write_Dict(result_dict, result_path)
        if args.assembly_mode == 'uce':
            Write_Uce_Summary([uce_summary_by_locus[key] for key in sorted(uce_summary_by_locus)], uce_summary_path)
        t1 = time.time()
        Write_Print(os.path.join(args.o,  "log.txt"), '\nTime cost:', t1 - t0, '\n') # 拼接所用的时间
    except Exception as e:
        Write_Print(os.path.join(args.o,  "log.txt"), "error:" , e)
        sys.exit(1)
