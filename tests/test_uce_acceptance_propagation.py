import csv
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from unix_command import combine_genes, uce_summary_row_is_accepted, write_uce_contigs_for_phyluce


def write_summary(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=['locus', 'status', 'accepted', 'low_quality'],
        )
        writer.writeheader()
        writer.writerows(rows)


class UceAcceptancePropagationTests(unittest.TestCase):
    def test_acceptance_reader_supports_new_and_legacy_summaries(self):
        self.assertTrue(uce_summary_row_is_accepted(
            {'status': 'success', 'accepted': '1', 'low_quality': '0'}))
        self.assertFalse(uce_summary_row_is_accepted(
            {'status': 'success', 'accepted': '0', 'low_quality': '0'}))
        self.assertTrue(uce_summary_row_is_accepted(
            {'status': 'success', 'low_quality': '0'}))
        self.assertFalse(uce_summary_row_is_accepted(
            {'status': 'low quality', 'low_quality': '1'}))

    def test_combine_and_phyluce_export_skip_rejected_stale_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / 'ref'
            out = root / 'out'
            results = out / 'sample_1' / 'results'
            ref.mkdir()
            results.mkdir(parents=True)

            (ref / 'accepted.fasta').write_text('>ref\nAAAA\n')
            (ref / 'rejected.fasta').write_text('>ref\nCCCC\n')
            (results / 'accepted.fasta').write_text('>contig_4_1_1_1_2\nAAAA\n')
            (results / 'rejected.fasta').write_text('>contig_4_1_1_1_2\nCCCC\n')
            write_summary(out / 'sample_1' / 'uce_assembly_summary.csv', [
                {'locus': 'accepted', 'status': 'success', 'accepted': '1', 'low_quality': '0'},
                {'locus': 'rejected', 'status': 'low quality', 'accepted': '0', 'low_quality': '1'},
            ])

            args = SimpleNamespace(
                o=str(out),
                r=str(ref),
                no_alignment=True,
                combine_source='assembly',
                assembly_mode='uce',
                p=1,
            )
            samples = {'sample_1': ('reads_1.fq', 'reads_2.fq')}

            combine_genes(args, samples)
            write_uce_contigs_for_phyluce(args, samples)

            self.assertTrue((out / 'combined_results' / 'accepted.fasta').is_file())
            self.assertFalse((out / 'combined_results' / 'rejected.fasta').exists())
            phyluce_text = (out / 'uce_contigs' / 'sample_1.contigs.fasta').read_text()
            self.assertIn('AAAA', phyluce_text)
            self.assertNotIn('CCCC', phyluce_text)


if __name__ == '__main__':
    unittest.main()
