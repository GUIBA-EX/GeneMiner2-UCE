from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import (
    build_fused_ucefilter_command,
    build_uce_rescue_filter_commands,
    get_uce_rescue_parallelism,
    get_uce_sample_parallelism,
    preserve_alignment_shadow,
    run_ordered_sample_stages,
    ucefilter_candidate_subdir,
    use_ucefilter,
)


class UceRescueParallelismTests(unittest.TestCase):
    def test_ucefilter_is_default_but_mito_and_legacy_are_excluded(self):
        self.assertTrue(use_ucefilter(SimpleNamespace(assembly_mode='uce')))
        self.assertEqual(ucefilter_candidate_subdir(SimpleNamespace(assembly_mode='uce')), 'filtered')
        self.assertFalse(use_ucefilter(SimpleNamespace(assembly_mode='original')))
        self.assertFalse(use_ucefilter(SimpleNamespace(assembly_mode='uce', is_mito_workflow=True)))
        self.assertFalse(use_ucefilter(SimpleNamespace(assembly_mode='uce', legacy_uce_filter=True)))

    def test_uce_pipeline_budgets_two_decode_workers_per_sample(self):
        self.assertEqual(get_uce_sample_parallelism(32, 45), (10, 1))
        self.assertEqual(get_uce_sample_parallelism(16, 45), (5, 1))
        self.assertEqual(get_uce_sample_parallelism(32, 2), (2, 1))
        self.assertEqual(get_uce_sample_parallelism(2, 45), (1, 1))
        self.assertEqual(get_uce_sample_parallelism(1, 45), (1, 1))

    def test_sample_pipeline_runs_stages_in_order_once(self):
        calls = []
        def record(stage):
            return lambda sample: calls.append((sample, stage))
        run_ordered_sample_stages(
            'sample_a', True, True, True, True,
            record('filter'), record('refilter'), record('assemble'), record('rescue'),
        )
        self.assertEqual(calls, [
            ('sample_a', 'filter'), ('sample_a', 'refilter'),
            ('sample_a', 'assemble'), ('sample_a', 'rescue'),
        ])

    def test_uses_four_by_four_when_threads_allow(self):
        self.assertEqual(get_uce_rescue_parallelism(32, 45), (4, 4))
        self.assertEqual(get_uce_rescue_parallelism(16, 45), (4, 4))

    def test_scales_down_for_small_thread_counts(self):
        self.assertEqual(get_uce_rescue_parallelism(8, 45), (2, 4))
        self.assertEqual(get_uce_rescue_parallelism(4, 45), (1, 4))
        self.assertEqual(get_uce_rescue_parallelism(2, 45), (1, 2))

    def test_scales_down_for_small_sample_counts(self):
        self.assertEqual(get_uce_rescue_parallelism(32, 2), (2, 4))

    def test_fused_ucefilter_uses_split_recruit_and_verify_references(self):
        args = SimpleNamespace(
            kf=21, step_size=4, max_reads=0,
            depth_low_water_mark=50, depth_limit=768, file_size_limit=6,
        )
        command = build_fused_ucefilter_command(
            'uce_filter', '/tmp/assembly_refs', '/tmp/terminal_baits',
            '/tmp/sample', '/tmp/R1.fq.gz', '/tmp/R2.fq.gz', args,
        )
        self.assertEqual(command[command.index('-r') + 1], '/tmp/assembly_refs')
        self.assertEqual(
            command[command.index('--recruit-references') + 1],
            '/tmp/terminal_baits',
        )
        self.assertEqual(command[command.index('--selection') + 1], 'auto')
        self.assertEqual(command[command.index('--reference-role') + 1], 'bait')
        self.assertEqual(command[command.index('--threads') + 1], '1')
        self.assertEqual(command[command.index('--memory-limit-mib') + 1], '256')

    def test_fused_ucefilter_marks_rescue_assembly_references_as_contigs(self):
        args = SimpleNamespace(
            kf=21, step_size=4, max_reads=0,
            depth_low_water_mark=50, depth_limit=768, file_size_limit=6,
        )
        command = build_fused_ucefilter_command(
            'uce_filter', '/tmp/assembly_refs', '/tmp/terminal_baits',
            '/tmp/sample', '/tmp/R1.fq.gz', '/tmp/R2.fq.gz', args,
            reference_role='contig',
        )
        self.assertEqual(command[command.index('--reference-role') + 1], 'contig')

    def test_legacy_uce_filter_retains_gm2_mode(self):
        args = SimpleNamespace(kf=21, step_size=4, max_reads=0, assembly_mode='uce', legacy_uce_filter=True)
        _, reads_cmd = build_uce_rescue_filter_commands(
            'MainFilterNew', '/tmp/rescue_refs', '/tmp/sample', '/tmp/R1.fq.gz',
            '/tmp/R2.fq.gz', args, '/tmp/sample/filter.dict',
        )
        self.assertEqual(reads_cmd[reads_cmd.index('-m') + 1], '5')
        self.assertEqual(reads_cmd[reads_cmd.index('-subdir') + 1], 'filtered_pe')

    def test_alignment_shadow_is_explicit_and_bounded(self):
        args = SimpleNamespace(
            kf=31, step_size=4, max_reads=0,
            depth_low_water_mark=50, depth_limit=768, file_size_limit=6,
            uce_alignment_shadow=True, uce_shadow_per_locus=16,
            uce_shadow_band=24, uce_shadow_terminal_window=120,
        )
        command = build_fused_ucefilter_command(
            'uce_filter', '/tmp/refs', '/tmp/refs', '/tmp/sample',
            '/tmp/R1.fq.gz', '/tmp/R2.fq.gz', args,
        )
        self.assertIn('--alignment-shadow', command)
        self.assertEqual(command[command.index('--shadow-per-locus') + 1], '16')
        self.assertEqual(command[command.index('--shadow-band') + 1], '24')
        self.assertEqual(command[command.index('--terminal-window') + 1], '120')

    def test_alignment_shadow_is_preserved_before_rescue_overwrite(self):
        with tempfile.TemporaryDirectory() as root:
            sample = Path(root) / 'sample'
            rescue = sample / 'uce_rescue_round_1'
            sample.mkdir()
            (sample / 'alignment_shadow.tsv').write_text('raw\n')
            (sample / 'alignment_shadow_summary.tsv').write_text('summary\n')
            preserve_alignment_shadow(str(sample), str(sample), '_initial')
            preserve_alignment_shadow(str(sample), str(rescue))
            self.assertEqual((sample / 'alignment_shadow_initial.tsv').read_text(), 'raw\n')
            self.assertEqual((rescue / 'alignment_shadow.tsv').read_text(), 'raw\n')
            self.assertEqual(
                (rescue / 'alignment_shadow_summary.tsv').read_text(), 'summary\n'
            )

    def test_rescue_filter_commands_preserve_max_reads(self):
        args = SimpleNamespace(kf=21, step_size=4, max_reads=1000, assembly_mode='uce')
        _, reads_cmd = build_uce_rescue_filter_commands(
            "MainFilterNew",
            "/tmp/rescue_refs",
            "/tmp/sample",
            "/tmp/R1.fq.gz",
            "/tmp/R2.fq.gz",
            args,
            "/tmp/sample/uce_rescue_kmer_dict_k21.dict",
        )

        self.assertEqual(reads_cmd[reads_cmd.index("-m_reads") + 1], "1000")


if __name__ == "__main__":
    unittest.main()
