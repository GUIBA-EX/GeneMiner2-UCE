from pathlib import Path
from types import SimpleNamespace
import csv
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import build_assembler_command, build_uce_rescue_refs


class UceRescueAssemblerReferenceTests(unittest.TestCase):
    def test_rescue_references_include_only_accepted_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / 'ref'
            sample = root / 'sample'
            results = sample / 'results'
            rescue = sample / 'uce_rescue_refs'
            ref.mkdir()
            results.mkdir(parents=True)

            (ref / 'accepted.fasta').write_text('>ref\nAAAA\n')
            (ref / 'rejected.fasta').write_text('>ref\nCCCC\n')
            (results / 'accepted.fasta').write_text('>contig\nAAAAGGGG\n')
            (results / 'rejected.fasta').write_text('>contig\nCCCCTTTT\n')

            with (sample / 'uce_assembly_summary.csv').open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=['locus', 'status', 'accepted', 'low_quality'])
                writer.writeheader()
                writer.writerow({'locus': 'accepted', 'status': 'success', 'accepted': '1', 'low_quality': '0'})
                writer.writerow({'locus': 'rejected', 'status': 'low quality', 'accepted': '0', 'low_quality': '1'})

            added = build_uce_rescue_refs(str(ref), str(sample), str(rescue), 4)

            self.assertEqual(added, 1)
            self.assertIn('AAAAGGGG', (rescue / 'accepted.fasta').read_text())
            self.assertNotIn('CCCCTTTT', (rescue / 'rejected.fasta').read_text())

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
            uce_path_strategy="search",
            uce_backbone_lookahead=7,
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
        self.assertEqual(cmd[cmd.index("--uce-path-strategy") + 1], "search")
        self.assertEqual(cmd[cmd.index("--uce-backbone-lookahead") + 1], "7")

        original_cmd = build_assembler_command(
            "main_assembler-original", args, "/tmp/out/sample",
            "/tmp/out/sample/uce_rescue_refs", "10000", 4, original=True,
        )
        self.assertNotIn("--uce-path-strategy", original_cmd)
        self.assertNotIn("--uce-backbone-lookahead", original_cmd)

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
