from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import unix_command


class MitoWorkflowTests(unittest.TestCase):
    def test_circular_consistency_accepts_rotation_and_reverse_complement_only(self):
        sequence = "ACGTTGCA" * 1500
        rotated = sequence[137:] + sequence[:137]
        self.assertTrue(unix_command.circularly_consistent(sequence, rotated))
        self.assertTrue(unix_command.circularly_consistent(sequence, "TGCAACGT" * 1500))
        complement = str.maketrans("ACGT", "TGCA")
        self.assertTrue(unix_command.circularly_consistent(sequence, sequence.translate(complement)[::-1]))

    def test_circular_consistency_rejects_one_base_difference(self):
        sequence = "ACGTTGCA" * 1500
        near_identical = "T" + sequence[1:]
        self.assertFalse(unix_command.circularly_consistent(sequence, near_identical))

    def test_circular_consistency_rejects_materially_different_sequence(self):
        sequence = "ACGTTGCA" * 1500
        different = "TGCATGCA" * 1500
        self.assertFalse(unix_command.circularly_consistent(sequence, different))

    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/mito_workflow")
    def test_prepare_passes_tiles_to_rust_helper(self, _find, run):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank = root / "mito.gb"
            genbank.write_text("LOCUS test\nORIGIN\n        1 acgt\n//\n")
            args = SimpleNamespace(
                mito_genbank=str(genbank),
                mito_flank=300,
                mito_tile_length=1200,
                mito_tile_step=600,
                o=str(root / "output"),
                r="",
            )
            unix_command.prepare_mito_reference(args)
            command = run.call_args.args[0]
            self.assertEqual(command[0:2], ["/gm2/mito_workflow", "prepare-reference"])
            self.assertEqual(command[command.index("--tile-length") + 1], "1200")
            self.assertEqual(command[command.index("--tile-step") + 1], "600")
            self.assertTrue(args.r.endswith(".gm2_mito_reference"))
            self.assertTrue(run.call_args.kwargs["check"])

    def test_mito_rescue_reference_keeps_all_baits_and_nonredundant_contigs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ref = root / "reference"
            sample = root / "sample"
            ref.mkdir()
            (sample / "contigs_all").mkdir(parents=True)
            (ref / "mitochondrion.fasta").write_text(">bait_a\nAAAA\n>bait_b\nCCCC\n")
            (sample / "contigs_all" / "mitochondrion.fasta").write_text(
                ">duplicate\nAAAA\n>seed_one\nGGGGTTTT\n>short\nAC\n>seed_two\nTTTTGGGG\n"
            )
            rescue = root / "rescue"
            self.assertEqual(unix_command.build_mito_rescue_refs(str(ref), str(sample), str(rescue), 4), 2)
            content = (rescue / "mitochondrion.fasta").read_text()
            self.assertIn(">bait_a", content)
            self.assertIn(">bait_b", content)
            self.assertIn(">mito_gm2_seed_2", content)
            self.assertIn(">mito_gm2_seed_4", content)
            self.assertEqual(content.count("AAAA"), 1)
            manifest = (rescue / "mito_rescue_seeds.tsv").read_text()
            self.assertIn("contig_index\tcontig_id\tsource_length\tsegment_index\tseed_length\tdecision", manifest)
            self.assertIn("1\tduplicate\t4\t1\t4\tduplicate", manifest)
            self.assertIn("3\tshort\t2\t1\t2\tshort", manifest)

    def test_mito_rescue_filter_uses_text_paired_output(self):
        args = SimpleNamespace(kf=25, step_size=1, max_reads=0)
        _, reads = unix_command.build_uce_rescue_filter_commands(
            "MainFilterNew", "refs", "sample", "r1.fq", "r2.fq", args, "dict", is_mito=True
        )
        self.assertEqual(reads[reads.index("-m") + 1], "4")

    def test_adaptive_stages_allow_partial_and_stop_after_two_exact_circles(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calls, policies, cache_flags = [], [], []
            args = SimpleNamespace(
                o=str(root / "output"), r=str(root / "reference"), max_reads=0,
                mito_initial_reads=1, mito_max_reads=2, reference_cache_dir=None, reuse_reference_cache=False,
            )

            def fake_filter(current_args, _samples, *_stages):
                calls.append(current_args.max_reads)
                cache_flags.append(current_args.reuse_reference_cache)

            def fake_finalize(current_args, _samples, require_circular=True):
                policies.append(require_circular)
                out = Path(current_args.o) / "sample" / "mito"
                out.mkdir(parents=True, exist_ok=True)
                sequence = "ACGTTGCA" * 100
                if current_args.max_reads == 2:
                    sequence = sequence[17:] + sequence[:17]
                (out / "mitochondrial_assembly_summary.tsv").write_text("metric\tvalue\nstatus\tcircular\n")
                (out / "mitochondrial_assembly.fasta").write_text(f">mito\n{sequence}\n")

            with mock.patch.object(unix_command, "do_filter_assemble", side_effect=fake_filter), \
                 mock.patch.object(unix_command, "write_uce_outputs"), \
                 mock.patch.object(unix_command, "run_mito_finalize", side_effect=fake_finalize):
                unix_command.run_mito_adaptive(args, {"sample": ("r1.fq", "r2.fq")})

            self.assertEqual(calls, [1, 2])
            self.assertEqual(policies, [False, False])
            self.assertEqual(cache_flags, [True, True])
            self.assertTrue((root / "output" / "sample" / "mito" / "mitochondrial_assembly.fasta").is_file())
            self.assertEqual(args.o, str(root / "output"))
            self.assertEqual(args.max_reads, 0)
            self.assertIsNone(args.reference_cache_dir)
            self.assertFalse(args.reuse_reference_cache)

    def test_adaptive_partial_stage_reaches_the_next_depth_before_failing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calls, policies = [], []
            args = SimpleNamespace(
                o=str(root / "output"), r=str(root / "reference"), max_reads=0,
                mito_initial_reads=1, mito_max_reads=2, reference_cache_dir=None, reuse_reference_cache=False,
            )

            def fake_filter(current_args, _samples, *_stages):
                calls.append(current_args.max_reads)

            def fake_finalize(current_args, _samples, require_circular=True):
                policies.append(require_circular)
                out = Path(current_args.o) / "sample" / "mito"
                out.mkdir(parents=True, exist_ok=True)
                status = "partial_multi_contig" if current_args.max_reads == 1 else "circular"
                (out / "mitochondrial_assembly_summary.tsv").write_text(f"metric\tvalue\nstatus\t{status}\n")
                (out / "mitochondrial_assembly.fasta").write_text(">mito\nACGTTGCAACGTTGCA\n")

            with mock.patch.object(unix_command, "do_filter_assemble", side_effect=fake_filter), \
                 mock.patch.object(unix_command, "write_uce_outputs"), \
                 mock.patch.object(unix_command, "run_mito_finalize", side_effect=fake_finalize):
                with self.assertRaisesRegex(RuntimeError, "stable circular assembly"):
                    unix_command.run_mito_adaptive(args, {"sample": ("r1.fq", "r2.fq")})

            self.assertEqual(calls, [1, 2])
            self.assertEqual(policies, [False, False])
            self.assertTrue((root / "output" / "sample" / "mito" / "mitochondrial_assembly.fasta").is_file())

    def test_adaptive_stable_partial_stops_without_deeper_scan(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calls = []
            args = SimpleNamespace(
                o=str(root / "output"), r=str(root / "reference"), max_reads=0,
                mito_initial_reads=1, mito_max_reads=8, reference_cache_dir=None, reuse_reference_cache=False,
            )

            def fake_filter(current_args, _samples, *_stages):
                calls.append(current_args.max_reads)

            def fake_finalize(current_args, _samples, require_circular=True):
                out = Path(current_args.o) / "sample" / "mito"
                out.mkdir(parents=True, exist_ok=True)
                (out / "mitochondrial_assembly_summary.tsv").write_text(
                    "metric\tvalue\nstatus\tpartial_multi_contig\n"
                )
                (out / "mitochondrial_assembly.fasta").write_text(">contig_1\nAAAA\n>contig_2\nCCCC\n")

            with mock.patch.object(unix_command, "do_filter_assemble", side_effect=fake_filter), \
                 mock.patch.object(unix_command, "write_uce_outputs"), \
                 mock.patch.object(unix_command, "run_mito_finalize", side_effect=fake_finalize):
                with self.assertRaisesRegex(RuntimeError, "stable non-circular assembly"):
                    unix_command.run_mito_adaptive(args, {"sample": ("r1.fq", "r2.fq")})

            self.assertEqual(calls, [1, 2])
            self.assertTrue((root / "output" / "sample" / "mito" / "mitochondrial_assembly.fasta").is_file())

    def test_mito_rescue_seed_qc_is_reference_agnostic_and_conservative(self):
        seen = set()
        distant = "ACGTGATCCGTAACGTTGCA"
        self.assertEqual(unix_command.mito_rescue_seed_reason(distant, 4, seen), "")
        seen.add(min(distant, unix_command.reverse_complement_dna(distant)))
        self.assertEqual(unix_command.mito_rescue_seed_reason("A" * 16, 4, seen), "uninformative_low_complexity")
        self.assertEqual(unix_command.mito_rescue_seed_segments("AACNCCGG"), ["AAC", "CCGG"])
        self.assertEqual(
            unix_command.mito_rescue_seed_reason(unix_command.reverse_complement_dna(distant), 4, seen),
            "duplicate",
        )

    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/mito_workflow")
    def test_finalize_passes_mate_link_policy_to_rust_helper(self, _find, run):
        args = SimpleNamespace(
            mito_genbank="mito.gb",
            mito_min_overlap=41,
            mito_min_overlap_identity=0.98,
            mito_terminal_window=500,
            mito_link_kmer=31,
            mito_min_link_hits=2,
            mito_min_pair_support=3,
            mito_bridge_kmer=31,
            mito_bridge_min_depth=2,
            mito_max_bridge=1000,
            r="output/.gm2_mito_reference",
            o="output",
        )
        unix_command.run_mito_finalize(args, {"1_Sample": ("r1.fq", "r2.fq")})
        self.assertEqual(run.call_count, 1)
        command = run.call_args.args[0]
        self.assertEqual(command[0:2], ["/gm2/mito_workflow", "finalize"])
        self.assertEqual(command[command.index("--minimum-overlap") + 1], "41")
        self.assertEqual(command[command.index("--minimum-identity") + 1], "0.98")
        self.assertEqual(command[command.index("--terminal-window") + 1], "500")
        self.assertEqual(command[command.index("--minimum-pair-support") + 1], "3")
        self.assertTrue(command[command.index("--contigs") + 1].endswith("contigs_all/mitochondrion.fasta"))
        self.assertTrue(command[command.index("--paired-reads") + 1].endswith("filtered/mitochondrion.fq"))
        self.assertEqual(command[command.index("--minimum-junction-support") + 1], "3")
        self.assertEqual(command[command.index("--require-circular") + 1], "true")
        self.assertTrue(run.call_args.kwargs["check"])


if __name__ == "__main__":
    unittest.main()
