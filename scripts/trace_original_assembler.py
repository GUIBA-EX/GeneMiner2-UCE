#!/usr/bin/env python3
"""测试工具：导出固定上游 assembler 的图和 seed 顺序。

不修改 ``main_assembler_original.py``，避免破坏它的上游逐字节来源保证。
本程序只复用其函数，输出与 Rust ``--trace-dir`` 可逐行比较的 TSV。
"""

import argparse
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

import main_assembler_original as original


def read_fastq_sequences(path, limit):
    """按原版 FASTQ 读取方式取前若干 read。"""
    sequences = []
    with open(path, encoding="utf-8", errors="ignore") as handle:
        while len(sequences) < limit:
            if not handle.readline():
                break
            sequence = handle.readline()
            handle.readline()
            handle.readline()
            sequences.append(''.join(filter(str.isalpha, sequence)).upper())
    return sequences


def write_read_trace(reads_path, kmer_size, output, read_limit):
    """输出正反向 k-mer 集合的数值排序，先排除插入顺序的影响。"""
    with open(output / "reads.tsv", "w") as handle:
        handle.write("read_index\torientation\tsequence\tkmers_sorted\n")
        for index, sequence in enumerate(read_fastq_sequences(reads_path, read_limit)):
            intseqs, read_len = original.Seq_To_Int(sequence)
            for orientation, encoded in zip(("forward", "reverse"), intseqs):
                mask = (1 << (kmer_size << 1)) - 1
                kmers = sorted((encoded >> (offset << 1)) & mask for offset in range(read_len - kmer_size))
                oriented = sequence if orientation == "forward" else original.Reverse_Complement_ACGT(sequence)
                handle.write(f"{index}\t{orientation}\t{oriented}\t{','.join(map(str, kmers))}\n")


def write_trace(reference, reads, kmer_size, limit, output, read_limit):
    """照原版顺序建图，落下 read、图节点和 seed 证据。"""
    output.mkdir(parents=True, exist_ok=True)
    write_read_trace(reads, kmer_size, output, read_limit)
    ref_dict, graph = {}, {}
    original.Make_Kmer_Dict(ref_dict, str(reference), kmer_size)
    original.Make_Assemble_Dict([str(reads)], kmer_size, graph, ref_dict)
    if limit > 0:
        graph = {key: value for key, value in graph.items() if value[0] > limit or value[3] > 0}
    with open(reference) as handle:
        ref_count = sum(1 for _ in SimpleFastaParser(handle))
    upper = 0
    if graph:
        q1, _, q3, _ = original.Quartile([value[0] for value in graph.values()])
        upper = int((q3 - q1) * 1.5 + q3)
        for value in graph.values():
            if value[3] != 0:
                value[3] = (value[0] > limit) * int(ref_count / (abs(value[3] - ref_count) + 1) * upper) + 1
            value[0] = min(value[0], upper)
    seeds = [(key, value[0], value[1], value[3], rank) for rank, (key, value) in enumerate(graph.items()) if value[1] > 1 and value[1] < 1000 and not value[2]]
    seeds.sort(key=lambda value: (value[3], value[1]), reverse=True)
    with open(output / "graph.tsv", "w") as handle:
        handle.write("rank\tkmer\tdepth\tposition\treverse\tref_depth\n")
        for rank, (key, value) in enumerate(graph.items()):
            handle.write(f"{rank}\t{key}\t{value[0]}\t{value[1]}\t{int(value[2])}\t{value[3]}\n")
    with open(output / "seeds.tsv", "w") as handle:
        handle.write("seed_rank\tkmer\tdepth\tposition\tref_depth\tgraph_rank\n")
        for rank, (key, depth, position, ref_depth, graph_rank) in enumerate(seeds):
            handle.write(f"{rank}\t{key}\t{depth}\t{position}\t{ref_depth}\t{graph_rank}\n")
    return len(graph), len(seeds), upper


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-r", "--reference", required=True, type=Path)
    parser.add_argument("-q", "--reads", required=True, type=Path)
    parser.add_argument("-ka", required=True, type=int)
    parser.add_argument("-limit_count", default=2, type=int)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--read-limit", default=20, type=int)
    args = parser.parse_args()
    graph_count, seed_count, upper = write_trace(args.reference, args.reads, args.ka, args.limit_count, args.output, args.read_limit)
    print(f"graph_nodes={graph_count} seeds={seed_count} depth_upper={upper}")


if __name__ == "__main__":
    main()
