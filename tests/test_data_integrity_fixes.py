from pathlib import Path
import io
import sys
import tempfile
import types
import unittest
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import build_trimed
import main_assembler
import main_refilter_new
import unix_command


def load_build_consensus():
    try:
        import build_consensus
        return build_consensus
    except ModuleNotFoundError:
        pyplot = types.ModuleType('matplotlib.pyplot')
        pyplot.plot = lambda *args, **kwargs: None
        pyplot.savefig = lambda *args, **kwargs: None
        matplotlib = types.ModuleType('matplotlib')
        matplotlib.pyplot = pyplot
        numpy = types.ModuleType('numpy')
        numpy.linspace = lambda *args, **kwargs: []
        stats = types.ModuleType('scipy.stats')
        stats.gaussian_kde = lambda values: None
        scipy = types.ModuleType('scipy')
        scipy.stats = stats
        sys.modules.update({
            'matplotlib': matplotlib,
            'matplotlib.pyplot': pyplot,
            'numpy': numpy,
            'scipy': scipy,
            'scipy.stats': stats,
        })
        import build_consensus
        return build_consensus


class DataIntegrityFixTests(unittest.TestCase):
    def test_reference_kmers_do_not_cross_ambiguous_bases(self):
        kmer_dict = main_refilter_new.build_kmer_dict({'AAAANCCCC'}, 4)
        self.assertEqual(len(kmer_dict), 4)

    def test_assembler_rejects_ambiguous_sequences_and_keeps_exact_k_read(self):
        self.assertEqual(main_assembler.Seq_To_Int('AAAANCCCC'), ((), 0))

        with tempfile.TemporaryDirectory() as tmp:
            reference = Path(tmp) / 'reference.fasta'
            reference.write_text('>ref\nAAAANCCCC\n')
            reference_kmers = {}
            main_assembler.Make_Kmer_Dict(reference_kmers, str(reference), 4)
            self.assertEqual(len(reference_kmers), 4)

            reads = Path(tmp) / 'reads.fq'
            reads.write_text('@r\nACGTA\n+\nFFFFF\n@amb\nAAAANCCCC\n+\nFFFFFFFFF\n')
            observed = {}
            main_assembler.Make_Assemble_Dict([str(reads)], 5, observed, {})
            self.assertGreaterEqual(len(observed), 2)

    def test_assembler_counts_repeated_and_palindromic_read_slices_once(self):
        unique_support = main_assembler.Calculate_Read_Support('CCCAAAGGG', 3, {'AAA': 2})
        self.assertEqual(unique_support, (2, 3, 1.0, 3, 6))

        read_count, supported_span, _, left, right = main_assembler.Calculate_Read_Support(
            'AAAAAAAAAA', 3, {'AAA': 1}
        )
        self.assertEqual(read_count, 1)
        self.assertEqual((supported_span, left, right), (0, 10, 0))

        with tempfile.TemporaryDirectory() as tmp:
            reads = Path(tmp) / 'reads.fq'
            reads.write_text('@palindrome\nATAT\n+\nFFFF\n')
            observed = {}
            self.assertEqual(main_assembler.Make_Reads_Dict([str(reads)], observed), 3)
            self.assertEqual(observed, {'ATA': 1})

    def test_consensus_preserves_insertions_and_all_references(self):
        consensus = load_build_consensus()
        sam = (
            '@HD\tVN:1.6\n'
            '@SQ\tSN:ref1\tLN:4\n'
            '@SQ\tSN:ref2\tLN:4\n'
            'r1\t0\tref1\t1\t60\t2M1I2M\t*\t0\t0\tACTGT\tFFFFF\n'
        )
        sequences, coverages, insertions, refname, _ = consensus.process_sam_header(io.StringIO(sam))
        consensus.parse_sam_file(io.StringIO(sam), sequences, coverages, insertions, refname, 150)
        self.assertEqual(insertions['ref1'], [(1, 'T')])

        sequences = consensus.reformat_sequences(sequences, coverages, insertions)
        with tempfile.TemporaryDirectory() as tmp:
            consensus.save_fastas(sequences, '-', coverages, tmp, 'indel', 1, insertions, 0, [0.75])
            self.assertIn('ACTGT', (Path(tmp) / 'indel.fasta').read_text())

        sequences = {
            'ref1': [[[1, ['A']]]],
            'ref2': [[[1, ['C']]]],
        }
        coverages = {'ref1': [1], 'ref2': [1]}
        insertions = {'ref1': {}, 'ref2': {}}
        with tempfile.TemporaryDirectory() as tmp:
            consensus.save_fastas(sequences, '-', coverages, tmp, 'out', 1, insertions, 0, [0.75])
            output = (Path(tmp) / 'out.fasta').read_text()
        self.assertIn('reference:ref1', output)
        self.assertIn('reference:ref2', output)

    def test_consensus_ignores_non_primary_alignments(self):
        consensus = load_build_consensus()
        sam = (
            '@SQ\tSN:ref\tLN:1\n'
            'primary\t0\tref\t1\t60\t1M\t*\t0\t0\tA\tF\n'
            'secondary\t256\tref\t1\t0\t1M\t*\t0\t0\t*\t*\n'
            'supplementary\t2048\tref\t1\t60\t1M\t*\t0\t0\tC\tF\n'
            'unmapped\t4\tref\t1\t0\t1M\t*\t0\t0\tG\tF\n'
        )
        sequences, coverages, insertions, refname, _ = consensus.process_sam_header(io.StringIO(sam))
        _, reads_mapped = consensus.parse_sam_file(
            io.StringIO(sam), sequences, coverages, insertions, refname, 150
        )
        self.assertEqual(reads_mapped, 1)
        self.assertEqual(sequences['ref'][0], {'-': 0, 'A': 1, 'C': 0, 'G': 0, 'N': 0, 'T': 0})

    def test_consensus_mapping_disables_secondary_alignments(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / '1_A' / 'results').mkdir(parents=True)
            (root / '1_A' / 'filtered').mkdir()
            (root / '1_A' / 'results' / 'gene.fasta').write_text('>gene\nAAAA\n')
            (root / '1_A' / 'filtered' / 'gene.fq').write_text('@r\nAAAA\n+\nFFFF\n')
            commands = []

            def fake_run(command, **kwargs):
                commands.append(command)
                if command[0] == 'minimap2':
                    Path(command[command.index('-o') + 1]).write_text('@SQ\tSN:gene\tLN:4\n')
                return SimpleNamespace(returncode=0)

            args = SimpleNamespace(o=str(root), r='references', consensus_threshold=0.75, p=1)
            with mock.patch.object(unix_command, 'find_executable', side_effect=lambda name, internal=False: name), \
                    mock.patch.object(unix_command, 'get_ref_genes', return_value=[('gene', '.fasta')]), \
                    mock.patch.object(unix_command.subprocess, 'run', side_effect=fake_run):
                unix_command.generate_consensus(args, {'1_A': ('reads.fq', 'reads.fq')})

            self.assertIn('--secondary=no', commands[0])

    def test_reverse_blast_coordinates_are_normalized(self):
        match = build_trimed.SequenceMatch.from_line(['q', 's', '99', '9', '0', '0', '10', '2', '1', '9', '0', '10'])
        self.assertEqual((match.qstart, match.qend, match.length, match.reverse), (2, 10, 9, True))

        with tempfile.TemporaryDirectory() as tmp:
            query = Path(tmp) / 'query.fasta'
            reference = Path(tmp) / 'reference.fasta'
            output = Path(tmp) / 'output.fasta'
            query.write_text('>query\nAACCGGTT\n')
            reference.write_text('>reference\nAACCGGTT\n')
            blast_line = 'q\ts\t100\t4\t0\t0\t8\t5\t1\t4\t0\t100\n'
            build_trimed.process_file(str(query), str(reference), iter([blast_line]), str(output), 0, 'longest')
            self.assertIn('\nAACC\n', output.read_text())

    def test_sam_padding_does_not_advance_reference(self):
        consensus = load_build_consensus()
        self.assertEqual(consensus.parse_cigar('1M1P1M', 'AC', 0), ('AC', [], 0))
        with self.assertRaisesRegex(ValueError, 'Malformed CIGAR'):
            consensus.parse_cigar('1Mbad', 'A', 0)

    def test_paired_read_length_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'different numbers'):
            list(main_refilter_new.linked_read_iterators([iter(['r1', 'r2']), iter(['r1'])]))

        with self.assertRaisesRegex(ValueError, 'odd number'):
            list(main_refilter_new.linked_iterator(iter(['r1']), 2))

    def test_rescue_revert_restores_matching_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            sample = Path(tmp) / 'sample'
            backup = Path(tmp) / 'backup'
            for root in (sample, backup):
                for subdir in ('results', 'contigs_all', 'contigs_all_low', 'filtered', 'filtered_pe'):
                    (root / subdir).mkdir(parents=True, exist_ok=True)

            (backup / 'results' / 'locus.fasta').write_text('>first\nAAAA\n')
            (sample / 'results' / 'locus.fasta').write_text('>rescue\nCCCC\n')
            (backup / 'filtered' / 'locus.fasta').write_text('>first_reads\nAAAA\n')
            (sample / 'filtered' / 'locus.fasta').write_text('>rescue_reads\nCCCC\n')
            (backup / 'filtered_pe' / 'locus_1.gm2').write_text('first')
            (sample / 'filtered_pe' / 'locus_1.gm2').write_text('rescue')
            (backup / 'filtered_pe' / 'locus_extra_1.gm2').write_text('other-first')
            (sample / 'filtered_pe' / 'locus_extra_1.gm2').write_text('other-rescue')
            (backup / 'ref_reads_count_dict.txt').write_text('locus,100,\n')
            (sample / 'ref_reads_count_dict.txt').write_text('locus,10,\nother,5,\n')

            before = {'locus': {'locus': 'locus', 'status': 'success', 'selected_contig_length': '100', 'read_count': '100',
                                'read_density': '1', 'support_fraction': '0.8', 'kmer_median_depth': '4',
                                'kmer_depth_cv': '0.2', 'kmer_max_depth_ratio': '2'}}
            after = {'locus': {'locus': 'locus', 'status': 'success', 'selected_contig_length': '100', 'read_count': '10'}}
            reverted = unix_command.revert_low_density_rescue_loci(str(sample), str(backup), before, after, 0.5)

            self.assertEqual(reverted, {'locus'})
            self.assertIn('first_reads', (sample / 'filtered' / 'locus.fasta').read_text())
            self.assertEqual((sample / 'filtered_pe' / 'locus_1.gm2').read_text(), 'first')
            self.assertEqual((sample / 'filtered_pe' / 'locus_extra_1.gm2').read_text(), 'other-rescue')
            self.assertIn('locus,100,', (sample / 'ref_reads_count_dict.txt').read_text())
            summary = (sample / 'uce_assembly_summary.csv').read_text()
            self.assertIn('read_density', summary.splitlines()[0])
            self.assertIn(',0.8,', summary)

    def test_locus_read_file_matching_is_exact(self):
        self.assertTrue(unix_command.locus_file_name_matches('gene.fq', 'gene'))
        self.assertFalse(unix_command.locus_file_name_matches('gene_1.fq', 'gene'))
        self.assertTrue(unix_command.locus_file_name_matches('gene_1.gm2', 'gene', paired=True))
        self.assertFalse(unix_command.locus_file_name_matches('gene_1_1.gm2', 'gene', paired=True))
        self.assertTrue(unix_command.locus_file_name_matches('gene_1_1.gm2', 'gene_1', paired=True))

    def test_assembler_resume_skips_only_with_existing_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for subdir in ('results', 'contigs_all', 'contigs_all_low', 'filtered'):
                (root / subdir).mkdir()
            (root / 'results' / 'gene.fasta').write_text('>contig\nAAAA\n')
            args = SimpleNamespace(o=str(root), ka=3, limit_count=0)

            skipped = main_assembler.process_key_value(args, 'gene', 'unused', 1, 1, 0, 1, 1, {'gene'})
            self.assertEqual(skipped[2]['status'], 'skipped')

            recomputed = main_assembler.process_key_value(args, 'gene', 'unused', 1, 1, 0, 1, 1, set())
            self.assertNotEqual(recomputed[2]['status'], 'skipped')


if __name__ == '__main__':
    unittest.main()
