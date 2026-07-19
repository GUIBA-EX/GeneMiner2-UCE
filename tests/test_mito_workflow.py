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
    def test_circular_consistency_accepts_rotation_and_near_identical_consensus(self):
        sequence = "ACGTTGCA" * 1500
        rotated = sequence[137:] + sequence[:137]
        self.assertTrue(unix_command.circularly_consistent(sequence, rotated))
        near_identical = rotated[:4000] + "A" + rotated[4001:]
        self.assertTrue(unix_command.circularly_consistent(sequence, near_identical))

    def test_circular_consistency_handles_a_difference_in_the_first_seed(self):
        sequence = "ACGTTGCA" * 1500
        near_identical = "T" + sequence[1:]
        self.assertTrue(unix_command.circularly_consistent(sequence, near_identical))

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
