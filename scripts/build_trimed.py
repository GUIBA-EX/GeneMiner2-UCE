from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import argparse
import os
import statistics
import subprocess

class SequenceMatch:
    __slots__ = ('pident', 'length', 'qstart', 'qend', 'reverse')

    def __init__(self, pident, qstart, qend, reverse=False):
        self.pident = pident
        self.qstart = min(qstart, qend)
        self.qend   = max(qstart, qend)
        self.length = self.qend - self.qstart + 1
        self.reverse = reverse

    @classmethod
    def from_line(cls, line):
        qstart, qend = int(line[6]), int(line[7])
        sstart, send = int(line[8]), int(line[9])
        reverse = (qstart > qend) != (sstart > send)
        return cls(float(line[2]), qstart, qend, reverse)

def execute_blastn(query_file, blast_db, executable_path=r"..\analysis\blastn.exe"):
    # Disable NCBI usage reporting to accelerate batch tasks.
    env = os.environ.copy()
    env['BLAST_USAGE_REPORT'] = '0'
    env['DO_NOT_TRACK'] = '1'

    # We do not use E-values because we want to match any possible homologous region.
    # Raw alignment score is probably the simpliest measurement of similarity
    # (confounded by alignment length, but we want long matches anyway).
    # At 99% hit rate, 85% identity and 200bp fragments, the maximal word size is 20.
    proc = subprocess.Popen([executable_path, "-query", query_file, "-db", blast_db,
                             "-outfmt", "6", "-word_size", "20", "-min_raw_gapped_score", "20"],
                            env=env, stdout=subprocess.PIPE, errors='replace', text=True)
    yield from proc.stdout
    returncode = proc.wait()

    if returncode:
        raise subprocess.CalledProcessError(returncode, proc.args)

def execute_magicblast(query_file, blast_db, executable_path=r"..\analysis\magicblast.exe"):
    # Disable NCBI usage reporting to accelerate batch tasks.
    env = os.environ.copy()
    env['BLAST_USAGE_REPORT'] = '0'
    env['DO_NOT_TRACK'] = '1'

    # At 99% hit rate, 85% identity and 50bp exons, the maximal word size is 13.
    proc = subprocess.Popen([executable_path, "-query", query_file, "-db", blast_db,
                             "-outfmt", "tabular", "-word_size", "13", "-score", "20",
                             "-limit_lookup", "F", "-penalty", "-2"],
                            env=env, stdout=subprocess.PIPE, errors='replace', text=True)
    yield from proc.stdout
    returncode = proc.wait()

    if returncode:
        raise subprocess.CalledProcessError(returncode, proc.args)

def process_file(query_file, ref_file, blast_output, output_file, percentage, criterion='all'):
    if os.path.isfile(output_file):
        os.remove(output_file)

    with open(query_file, 'r') as f:
        try:
            header, sequence = next(SimpleFastaParser(f))
        except StopIteration:
            return

    with open(ref_file, 'r') as f:
        try:
            median_length = statistics.median(len(seq) for _, seq in SimpleFastaParser(f))
        except statistics.StatisticsError:
            return

    if criterion == 'all':
        matches = merge_matches(blast_output)
    else:
        matches = read_matches(blast_output)

    if not matches:
        return

    # A single output sequence cannot combine HSPs from opposing orientations.
    # Retain the orientation with the greatest total aligned query length.
    orientation = max((False, True), key=lambda reverse: sum(m.length for m in matches if m.reverse == reverse))
    matches = [m for m in matches if m.reverse == orientation]

    if criterion == 'longest':
        matches = [max(matches, key=lambda match: match.length)]
    elif criterion == 'terminal':
        matches = [SequenceMatch(100, min(m.qstart for m in matches), max(m.qend for m in matches), orientation)]

    matches.sort(key=lambda match: match.qstart)
    combined_sequence = ''.join(sequence[m.qstart - 1:m.qend] for m in matches)

    if orientation:
        combined_sequence = str(Seq(combined_sequence).reverse_complement())

    if len(combined_sequence) / median_length * 100 <= percentage:
        return

    with open(output_file, 'w') as f:
        f.write(f'>{header}\n{combined_sequence}\n')

def read_matches(blast_out, sorting_key=lambda match: match.length):
    return sorted((SequenceMatch.from_line(parts)
                   for parts
                   in (line.split("\t")
                       for line
                       in blast_out
                       if not line.startswith('#'))
                   if len(parts) >= 12),
                   key=sorting_key,
                   reverse=True)

def merge_matches(blast_out):
    matches = read_matches(blast_out, lambda match: match.qstart)
    matches.reverse()

    merged_matches = []

    for reverse in (False, True):
        oriented_matches = [match for match in matches if match.reverse == reverse]
        if not oriented_matches:
            continue

        current_match = oriented_matches[0]

        for match in oriented_matches[1:]:
            if match.qstart <= current_match.qend:
                current_match.qend = max(current_match.qend, match.qend)
                current_match.length = current_match.qend - current_match.qstart + 1
            else:
                merged_matches.append(current_match)
                current_match = match

        merged_matches.append(current_match)

    return merged_matches

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim contigs using BLAST")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-r", "--ref", required=True, help="Reference file")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument("-b", "--blast-db", required=True, help="BLAST database path")
    parser.add_argument("-m", "--mode", type=int, required=False, default=2, help="0: all fragments; 1: longest match; 2: trim terminal; 3: best spliced match")
    parser.add_argument("-p", "--pec", type=int, required=False, default=50, help="Retention percentage")

    args = parser.parse_args()

    if args.mode == 1:
        blast_iter = execute_blastn
        criterion = 'longest'
    elif args.mode == 2:
        blast_iter = execute_blastn
        criterion = 'terminal'
    elif args.mode == 3:
        blast_iter = execute_magicblast
        criterion = 'longest'
    else: # mode == 0
        blast_iter = execute_blastn
        criterion = 'all'

    process_file(args.input, args.ref, blast_iter(args.input, args.blast_db), args.output, args.pec, criterion)
