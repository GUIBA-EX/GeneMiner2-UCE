from pathlib import Path
from types import SimpleNamespace
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import build_assembler_command


class UceRescueAssemblerReferenceTests(unittest.TestCase):
    def test_assembler_command_uses_supplied_reference_dir(self):
        args = SimpleNamespace(
            ka=0,
            min_ka=17,
            max_ka=31,
            error_threshold=1,
            search_depth=4096,
            min_coverage=0,
            assembly_mode="uce",
            uce_side_candidates=8,
            uce_max_contig_length=5000,
            uce_min_read_density=0.003,
            uce_density_check_min_length=1000,
            uce_max_depth_cv=0,
            uce_max_depth_ratio=0,
            assembler_reference_cache_dir=None,
        )

        cmd = build_assembler_command(
            "main_assembler",
            args,
            "/tmp/out/sample",
            "/tmp/out/sample/uce_rescue_refs",
            "10000",
            4,
        )

        self.assertEqual(cmd[cmd.index("-r") + 1], "/tmp/out/sample/uce_rescue_refs")
        self.assertEqual(cmd[cmd.index("-o") + 1], "/tmp/out/sample")
        self.assertEqual(cmd[cmd.index("-p") + 1], "4")
        self.assertEqual(cmd[cmd.index("--uce-max-contig-length") + 1], "5000")

    def test_assembler_command_forwards_reference_cache_dir(self):
        args = SimpleNamespace(
            ka=0,
            min_ka=17,
            max_ka=31,
            error_threshold=1,
            search_depth=4096,
            min_coverage=0,
            assembly_mode="uce",
            uce_side_candidates=8,
            uce_max_contig_length=5000,
            uce_min_read_density=0.003,
            uce_density_check_min_length=1000,
            uce_max_depth_cv=0,
            uce_max_depth_ratio=0,
            assembler_reference_cache_dir="/tmp/cache/assembler",
            r="/tmp/ref",
        )

        cmd = build_assembler_command(
            "main_assembler",
            args,
            "/tmp/out/sample",
            "/tmp/ref",
            "10000",
            4,
        )

        self.assertEqual(
            cmd[cmd.index("--assembler-reference-cache-dir") + 1],
            "/tmp/cache/assembler",
        )

    def test_assembler_command_does_not_cache_dynamic_rescue_refs(self):
        args = SimpleNamespace(
            ka=0,
            min_ka=17,
            max_ka=31,
            error_threshold=1,
            search_depth=4096,
            min_coverage=0,
            assembly_mode="uce",
            uce_side_candidates=8,
            uce_max_contig_length=5000,
            uce_min_read_density=0.003,
            uce_density_check_min_length=1000,
            uce_max_depth_cv=0,
            uce_max_depth_ratio=0,
            assembler_reference_cache_dir="/tmp/cache/assembler",
            r="/tmp/ref",
        )

        cmd = build_assembler_command(
            "main_assembler",
            args,
            "/tmp/out/sample",
            "/tmp/out/sample/uce_rescue_refs",
            "10000",
            4,
        )

        self.assertNotIn("--assembler-reference-cache-dir", cmd)


if __name__ == "__main__":
    unittest.main()
