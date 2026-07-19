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
        self.assertIn("{auto,uce-rust,original,original-rust}", proc.stdout)
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

    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/gm2_stats")
    def test_single_end_stats_path_is_not_counted_twice(self, _find, run):
        args = SimpleNamespace(
            o="out", r="refs", stats_count_input_reads=True, stats_no_heatmap=True
        )
        unix_command.run_stats(args, {"1_A": ("reads.fq.gz", "reads.fq.gz")})
        command = run.call_args.args[0]
        sample = command.index("--sample")
        self.assertEqual(command[sample + 1:sample + 4], ["1_A", "reads.fq.gz", ""])

    def test_failed_samples_are_excluded_from_rescue(self):
        samples = {"1_A": (), "2_B": (), "3_C": ()}
        failures = [("2_B", "assemble", "failed")]
        self.assertEqual(
            unix_command.get_rescue_sample_names(samples, failures),
            ["1_A", "3_C"],
        )

    @mock.patch.object(unix_command, "run_stats", side_effect=OSError("broken file"))
    @mock.patch.object(unix_command.os.path, "isdir", return_value=True)
    def test_execute_tasks_handles_file_errors(self, _isdir, _stats):
        args = SimpleNamespace(r="refs", command=("stats",))
        self.assertEqual(unix_command.execute_tasks(args, {"1_A": ()}), 1)

    @mock.patch.object(unix_command, "write_failed_samples")
    @mock.patch.object(unix_command.os.path, "isdir", return_value=True)
    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/main_assembler-original-rust")
    def test_reference_auto_uses_original_rust_by_default(
        self, find, run, _isdir, _write_failures
    ):
        args = SimpleNamespace(
            o="out", r="refs", p=1, kf=21, step_size=4,
            reuse_reference_cache=False, reference_cache_dir=None,
            uce_rescue_reads=False, assembly_mode="original", soft_boundary="auto",
            assembler_implementation="auto", ka=39, min_ka=21, max_ka=39,
            error_threshold=2, search_depth=4096, min_coverage=0,
            uce_side_candidates=8, uce_max_contig_length=0,
            uce_min_read_density=0.003, uce_density_check_min_length=1000,
            uce_max_depth_cv=0, uce_max_depth_ratio=0,
        )

        with self.assertRaises(RuntimeError):
            unix_command.do_filter_assemble(
                args, {"1_A": ("r1.fq", "r2.fq")}, False, False, True
            )

        find.assert_called_once_with("main_assembler-original-rust", internal=True)
        command = run.call_args.args[0]
        self.assertEqual(command[0], "/gm2/main_assembler-original-rust")
        self.assertNotIn("--assembly-mode", command)
        self.assertNotIn("--uce-side-candidates", command)


    @mock.patch.object(unix_command, "get_reference_kmer_dict_path", return_value="/cache/filter.dict")
    @mock.patch.object(unix_command, "write_failed_samples")
    @mock.patch.object(unix_command.os.path, "isdir", return_value=True)
    @mock.patch.object(unix_command.subprocess, "run")
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/main_assembler-original-rust")
    def test_reference_original_rust_receives_versioned_cache(self, find, run, _isdir, _write_failures, _filter_cache):
        args = SimpleNamespace(
            o="out", r="refs", p=1, kf=21, step_size=4,
            reuse_reference_cache=True, reference_cache_dir="/cache",
            uce_rescue_reads=False, assembly_mode="original", soft_boundary="auto",
            assembler_implementation="original-rust", ka=39, min_ka=21, max_ka=39,
            error_threshold=2, search_depth=4096, min_coverage=0,
            uce_side_candidates=8, uce_max_contig_length=0,
            uce_min_read_density=0.003, uce_density_check_min_length=1000,
            uce_max_depth_cv=0, uce_max_depth_ratio=0,
        )
        with self.assertRaises(RuntimeError):
            unix_command.do_filter_assemble(
                args, {"1_A": ("r1.fq", "r2.fq")}, False, False, True
            )
        find.assert_called_once_with("main_assembler-original-rust", internal=True)
        command = run.call_args.args[0]
        self.assertEqual(command[0], "/gm2/main_assembler-original-rust")
        self.assertEqual(
            command[command.index("--assembler-reference-cache-dir") + 1],
            "/cache/assembler",
        )
        self.assertNotIn("--assembly-mode", command)
        self.assertNotIn("--uce-side-candidates", command)
    @mock.patch.object(unix_command, "write_failed_samples")
    @mock.patch.object(unix_command.os.path, "isdir", return_value=True)
    @mock.patch.object(
        unix_command.subprocess,
        "run",
        side_effect=subprocess.CalledProcessError(1, ["main_assembler-rust"]),
    )
    @mock.patch.object(unix_command, "find_executable", return_value="/gm2/main_assembler-rust")
    def test_uce_rust_failures_do_not_fall_back_to_python(
        self, find, run, _isdir, _write_failures
    ):
        for assembly_mode in ("uce",):
            with self.subTest(assembly_mode=assembly_mode):
                find.reset_mock()
                run.reset_mock()
                args = SimpleNamespace(
                    o="out", r="refs", p=1, kf=21, step_size=4,
                    reuse_reference_cache=False, reference_cache_dir=None,
                    uce_rescue_reads=False, assembly_mode=assembly_mode,
                    soft_boundary="auto", assembler_implementation="auto",
                    ka=21, min_ka=21, max_ka=21,
                    error_threshold=2, search_depth=4096, min_coverage=0,
                    uce_side_candidates=8, uce_max_contig_length=0,
                    uce_min_read_density=0.003,
                    uce_density_check_min_length=1000,
                    uce_max_depth_cv=0, uce_max_depth_ratio=0,
                )

                with self.assertRaises(RuntimeError):
                    unix_command.do_filter_assemble(
                        args, {"1_A": ("r1.fq", "r2.fq")}, False, False, True
                    )

                find.assert_called_once_with("main_assembler-rust", internal=True)
                self.assertEqual(run.call_count, 1)
                self.assertEqual(run.call_args.args[0][0], "/gm2/main_assembler-rust")

if __name__ == "__main__":
    unittest.main()
