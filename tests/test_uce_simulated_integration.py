import csv
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from main_assembler import Calculate_Read_Support, Evaluate_Uce_Candidate


def random_dna(length, seed=1731):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(length))


def simulated_reads(sequence, read_length=150, step=10):
    reads = []
    for start in range(0, len(sequence) - read_length + 1, step):
        read = sequence[start:start + read_length]
        reads.append(read)

        # Add one deterministic low-frequency sequencing error while retaining
        # an exact read at every start, so the true path remains identifiable.
        error_pos = (start * 7 + 23) % read_length
        replacement = next(base for base in 'ACGT' if base != read[error_pos])
        reads.append(read[:error_pos] + replacement + read[error_pos + 1:])
    return reads


def write_fastq(path, reads):
    with path.open('w') as handle:
        for index, read in enumerate(reads):
            handle.write(f'@sim_{index}\n{read}\n+\n{"I" * len(read)}\n')


class UceSimulatedIntegrationTests(unittest.TestCase):
    def test_uce_assembler_recovers_known_locus(self):
        truth = random_dna(600)
        reads = simulated_reads(truth)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            refs = root / 'refs'
            sample = root / 'sample'
            filtered = sample / 'filtered'
            refs.mkdir()
            filtered.mkdir(parents=True)
            (refs / 'uce_sim.fasta').write_text(f'>uce_sim\n{truth}\n')
            write_fastq(filtered / 'uce_sim.fq', reads)

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / 'scripts' / 'main_assembler.py'),
                    '-r', str(refs),
                    '-o', str(sample),
                    '-ka', '31',
                    '-limit_count', '1',
                    '-iteration', '4096',
                    '-p', '1',
                    '--assembly-mode', 'uce',
                    '--uce-max-contig-length', '2000',
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )

            with (sample / 'uce_assembly_summary.csv').open(newline='') as handle:
                row = next(csv.DictReader(handle))

            self.assertEqual(row['status'], 'success')
            self.assertEqual(row['accepted'], '1')
            self.assertGreater(float(row['slice_support_breadth']), 0.9)
            self.assertGreater(int(row['unique_read_count']), 0)
            self.assertEqual(int(row['multi_mapping_read_count']), 0)

            with (sample / 'results' / 'uce_sim.fasta').open() as handle:
                _, recovered = next(SimpleFastaParser(handle))

            self.assertEqual(len(recovered), len(truth))
            self.assertIn(recovered, {truth, str(Seq(truth).reverse_complement())})

    def test_sparse_and_repetitive_support_do_not_look_continuous(self):
        candidate = 'A' * 100 + 'C' * 800 + 'G' * 100
        support = Calculate_Read_Support(
            candidate,
            100,
            {'A' * 100: 2, 'G' * 100: 2, 'C' * 100: 20},
        )

        self.assertEqual(support.supported_extent, 1000)
        self.assertEqual(support.supported_bases, 200)
        self.assertAlmostEqual(support.breadth, 0.2)
        self.assertEqual(support.max_gap, 800)
        self.assertEqual(support.multi_mapping_read_count, 20)

        accepted, reasons = Evaluate_Uce_Candidate(
            contig_len=len(candidate),
            unique_read_count=0,
            supported_bases=0,
            unique_read_density=0,
            depth_cv=0,
            max_depth_ratio=1,
            guardrails=None,
        )
        self.assertFalse(accepted)
        self.assertEqual(reasons, ['no_unique_read_support', 'no_positional_support'])


if __name__ == '__main__':
    unittest.main()
