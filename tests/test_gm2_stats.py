import csv
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import gm2_stats


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_summary(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        'locus',
        'status',
        'selected_contig_length',
        'read_supported_span',
        'read_count',
        'flank_balance',
        'candidate_count',
        'low_quality',
    ]

    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


class GeneMinerStatsTests(unittest.TestCase):
    def test_build_stats_from_uce_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / 'ref'
            out = root / 'out'

            write_text(ref / 'uce-1.fasta', '>a\nAAAA\n>b\nAAAAAA\n')
            write_text(ref / 'uce-2.fasta', '>a\nAAAAAAAA\n')
            write_summary(out / '1_A' / 'uce_assembly_summary.csv', [
                {
                    'locus': 'uce-1',
                    'status': 'success',
                    'selected_contig_length': '5',
                    'read_supported_span': '4',
                    'read_count': '10',
                    'flank_balance': '0.5',
                    'candidate_count': '2',
                    'low_quality': '0',
                },
                {
                    'locus': 'uce-2',
                    'status': 'low quality',
                    'selected_contig_length': '4',
                    'read_supported_span': '2',
                    'read_count': '1',
                    'flank_balance': '0.2',
                    'candidate_count': '1',
                    'low_quality': '1',
                },
            ])
            write_text(out / '1_A' / 'ref_reads_count_dict.txt', 'uce-1,12\nuce-2,3\n')

            refs, loci, sample_names, matrices, sample_stats, locus_stats = gm2_stats.build_stats(
                str(out),
                str(ref),
                {'1_A': ('reads_1.fq.gz', 'reads_2.fq.gz')},
            )

            self.assertEqual(refs, {'uce-1': 5, 'uce-2': 8})
            self.assertEqual(loci, ['uce-1', 'uce-2'])
            self.assertEqual(sample_names, ['1_A'])
            self.assertEqual(matrices['lengths']['1_A']['uce-1'], 5)
            self.assertEqual(matrices['read_counts']['1_A']['uce-2'], 1)
            self.assertEqual(matrices['filtered_counts']['1_A']['uce-1'], 12)
            self.assertEqual(sample_stats[0]['LociSuccess'], 1)
            self.assertEqual(sample_stats[0]['LociLowQuality'], 1)
            self.assertEqual(sample_stats[0]['LociAt75pct'], 1)
            self.assertEqual(locus_stats[0]['SuccessSamples'], 1)

    def test_run_writes_core_tables_without_heatmaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / 'ref'
            out = root / 'out'

            write_text(ref / 'uce-1.fasta', '>a\nAAAA\n')
            write_summary(out / '1_A' / 'uce_assembly_summary.csv', [
                {
                    'locus': 'uce-1',
                    'status': 'success',
                    'selected_contig_length': '4',
                    'read_supported_span': '4',
                    'read_count': '8',
                    'flank_balance': '1.0',
                    'candidate_count': '1',
                    'low_quality': '0',
                },
            ])
            write_text(out / '1_A' / 'ref_reads_count_dict.txt', 'uce-1,8\n')

            args = SimpleNamespace(
                o=str(out),
                r=str(ref),
                stats_count_input_reads=False,
                stats_no_heatmap=True,
            )
            gm2_stats.run(args, {'1_A': ('reads_1.fq.gz', 'reads_2.fq.gz')})

            self.assertTrue((out / 'uce_seq_lengths.tsv').is_file())
            self.assertTrue((out / 'uce_stats.tsv').is_file())
            self.assertTrue((out / 'uce_read_counts.tsv').is_file())
            self.assertTrue((out / 'uce_filtered_read_counts.tsv').is_file())
            self.assertTrue((out / 'uce_locus_stats.tsv').is_file())


if __name__ == '__main__':
    unittest.main()
