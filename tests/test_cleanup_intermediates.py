from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
mock = unittest.mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import unix_command


class CleanupIntermediatesTests(unittest.TestCase):
    def test_uce_cleanup_keeps_final_outputs_and_audit_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            sample = output / "sample"
            (sample / "filtered").mkdir(parents=True)
            (sample / "filtered" / "locus.fq").write_text("reads")
            (sample / "filtered_pe").mkdir()
            (sample / "filtered_pe" / "locus_1.fq").write_text("reads")
            (sample / "results").mkdir()
            (sample / "results" / "locus.fasta").write_text(">locus\nACGT\n")
            (sample / "contigs_all").mkdir()
            (sample / "contigs_all" / "locus.fasta").write_text(">locus\nACGT\n")
            round_root = sample / "uce_rescue_round_1"
            (round_root / "assembly_refs").mkdir(parents=True)
            (round_root / "terminal_baits").mkdir()
            (round_root / "filter_k31.dict").write_text("index")
            (sample / "uce_terminal_rescue_diagnostics.tsv").write_text("audit")
            args = SimpleNamespace(o=str(output), is_mito_workflow=False, assembly_mode="uce", uce_rescue_reads=True)

            unix_command.cleanup_intermediates(args, {"sample": ("r1", "r2")})

            self.assertFalse((sample / "filtered").exists())
            self.assertFalse((sample / "filtered_pe").exists())
            self.assertFalse((round_root / "assembly_refs").exists())
            self.assertFalse((round_root / "terminal_baits").exists())
            self.assertFalse((round_root / "filter_k31.dict").exists())
            self.assertTrue((sample / "results" / "locus.fasta").is_file())
            self.assertTrue((sample / "contigs_all" / "locus.fasta").is_file())
            self.assertTrue((sample / "uce_terminal_rescue_diagnostics.tsv").is_file())
            manifest = (output / "cleanup_manifest.tsv").read_text()
            self.assertIn("sample/filtered", manifest)
            self.assertIn("uce_rescue_round_1/assembly_refs", manifest)

    def test_mito_cleanup_removes_only_finalization_input_and_adaptive_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            sample = output / "sample"
            (sample / "filtered").mkdir(parents=True)
            (sample / "filtered" / "mitochondrion.fq").write_text("reads")
            (sample / "filtered" / "keep.fq").write_text("keep")
            (sample / "mito").mkdir()
            (sample / "mito" / "mitochondrial_assembly.fasta").write_text(">mito\nACGT\n")
            (output / ".mito_adaptive" / "10m").mkdir(parents=True)
            (output / ".mito_adaptive" / "10m" / "scratch").write_text("stage")
            (output / ".gm2_mito_reference").mkdir()
            args = SimpleNamespace(o=str(output), is_mito_workflow=True, assembly_mode="uce", uce_rescue_reads=True)

            unix_command.cleanup_intermediates(args, {"sample": ("r1", "r2")})

            self.assertFalse((sample / "filtered" / "mitochondrion.fq").exists())
            self.assertTrue((sample / "filtered" / "keep.fq").is_file())
            self.assertTrue((sample / "mito" / "mitochondrial_assembly.fasta").is_file())
            self.assertFalse((output / ".mito_adaptive").exists())
            self.assertTrue((output / ".gm2_mito_reference").is_dir())

    def test_execute_tasks_cleans_only_after_later_consumers_succeed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference"
            reference.mkdir()
            args = SimpleNamespace(
                command=("filter", "assemble", "consensus"),
                r=str(reference), o=str(root / "output"), assembly_mode="original",
                is_mito_workflow=False, is_gene_workflow=False,
                cleanup_intermediates=True,
            )
            order = []
            with mock.patch.object(unix_command, "do_filter_assemble", side_effect=lambda *args: order.append("assemble")), \
                 mock.patch.object(unix_command, "generate_consensus", side_effect=lambda *args: order.append("consensus")), \
                 mock.patch.object(unix_command, "cleanup_intermediates", side_effect=lambda *args: order.append("cleanup")):
                self.assertEqual(unix_command.execute_tasks(args, {"sample": ("r1", "r2")}), 0)
            self.assertEqual(order, ["assemble", "consensus", "cleanup"])

    def test_execute_tasks_remains_compatible_without_the_cleanup_attribute(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference"
            reference.mkdir()
            args = SimpleNamespace(r=str(reference), command=("stats",))
            with mock.patch.object(unix_command, "run_stats") as stats:
                self.assertEqual(unix_command.execute_tasks(args, {"sample": ("r1", "r2")}), 0)
            stats.assert_called_once()

    def test_execute_tasks_does_not_clean_after_a_failed_workflow(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference"
            reference.mkdir()
            args = SimpleNamespace(
                command=("filter", "assemble"),
                r=str(reference), o=str(root / "output"), assembly_mode="original",
                is_mito_workflow=False, is_gene_workflow=False,
                cleanup_intermediates=True,
            )
            with mock.patch.object(unix_command, "do_filter_assemble", side_effect=RuntimeError("failed")), \
                 mock.patch.object(unix_command, "cleanup_intermediates") as cleanup:
                self.assertEqual(unix_command.execute_tasks(args, {"sample": ("r1", "r2")}), 1)
            cleanup.assert_not_called()

    def test_cleanup_skips_a_symlink_instead_of_following_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            sample = output / "sample"
            sample.mkdir(parents=True)
            external = root / "external"
            external.mkdir()
            (external / "reads.fq").write_text("reads")
            (sample / "filtered").symlink_to(external, target_is_directory=True)
            args = SimpleNamespace(o=str(output), is_mito_workflow=False, assembly_mode="original", uce_rescue_reads=False)

            unix_command.cleanup_intermediates(args, {"sample": ("r1", "r2")})

            self.assertTrue((external / "reads.fq").is_file())
            self.assertTrue((sample / "filtered").is_symlink())
            self.assertIn("skipped_symlink", (output / "cleanup_manifest.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
