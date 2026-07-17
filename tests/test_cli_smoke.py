from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import unix_command


ROOT = Path(__file__).resolve().parents[1]


class CliSmokeTests(unittest.TestCase):
    def test_legacy_script_help(self):
        proc = subprocess.run(
            [sys.executable, "scripts/unix_command.py", "-h"],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("--assembly-mode", proc.stdout)
        self.assertIn("--assembler-implementation", proc.stdout)
        self.assertIn("{auto,rust,original}", proc.stdout)
        self.assertIn("--uce-path-strategy", proc.stdout)
        self.assertIn("--uce-backbone-lookahead", proc.stdout)
        self.assertIn("--uce-rescue-reads", proc.stdout)
        self.assertIn("--uce-rescue-min-density-ratio", proc.stdout)
        self.assertIn("--uce-max-contig-length", proc.stdout)
        self.assertIn("--reuse-reference-cache", proc.stdout)
        self.assertIn("--population-minibwa", proc.stdout)
        self.assertIn("--population-min-call-rate", proc.stdout)
        self.assertIn("--population-reference-strategy", proc.stdout)
        self.assertIn("--population-ld-r2", proc.stdout)
        self.assertIn("--population-admixture-k-max", proc.stdout)
        self.assertIn("population", proc.stdout)

    def test_default_alifilter_model_uses_builtin_model(self):
        for value in (None, "", "default", " DEFAULT "):
            args = SimpleNamespace(alifilter_model=value)
            self.assertIsNone(unix_command.get_alifilter_model(args))

    def test_custom_alifilter_model_is_preserved(self):
        args = SimpleNamespace(alifilter_model="/path/to/model.json")
        self.assertEqual(unix_command.get_alifilter_model(args), "/path/to/model.json")

    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/main_population")
    def test_population_command_invokes_rust_driver_with_minibwa(self, _find, run):
        args = SimpleNamespace(
            o="out", f="samples.tsv", p=4,
            population_reference_strategy="sqcl-longest",
            population_min_mapq=20, population_min_baseq=20,
            population_min_dp=5, population_min_gq=20,
            population_min_qual=20.0, population_min_call_rate=0.8,
            population_min_mac=2, population_start_at="reference", population_stop_after="selection",
            population_reference_fasta=None,
            population_ld_window=50, population_ld_step=5, population_ld_r2=0.2,
            population_admixture_k_min=2, population_admixture_k_max=6,
            population_admixture_cv=10,
            population_minibwa="minibwa", population_samtools="samtools",
            population_bcftools="bcftools", population_plink="plink",
            population_admixture="admixture",
            population_skip_mark_duplicates=False, population_skip_plink=True,
            population_skip_admixture=False,
        )
        unix_command.run_population(args)
        command = run.call_args.args[0]
        self.assertEqual(command[0], "/gm2/main_population")
        self.assertIn("--minibwa", command)
        self.assertEqual(command[command.index("--minibwa") + 1], "minibwa")
        self.assertEqual(command[command.index("--reference-strategy") + 1], "sqcl-longest")
        self.assertEqual(command[command.index("--start-at") + 1], "reference")
        self.assertEqual(command[command.index("--ld-r2") + 1], "0.2")
        self.assertEqual(command[command.index("--admixture") + 1], "admixture")
        self.assertIn("--skip-plink", command)
        self.assertTrue(run.call_args.kwargs["check"])

if __name__ == "__main__":
    unittest.main()
