from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import (
    do_filter_assemble,
    get_assembler_reference_cache_dir,
    get_reference_kmer_dict_path,
    reference_cache_key,
)
from main_assembler import Load_Or_Make_Kmer_Dict


class ReferenceCacheTests(unittest.TestCase):
    def test_reference_cache_path_uses_fingerprint_when_enabled(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "ref"
            ref.mkdir()
            (ref / "uce-1.fasta").write_text(">a\nACGTACGT\n")

            args = SimpleNamespace(
                reuse_reference_cache=True,
                reference_cache_dir=None,
                r=str(ref),
                kf=21,
                step_size=4,
            )

            path = get_reference_kmer_dict_path(args, str(root / "out"))

            self.assertIn(".gm2_reference_cache", path)
            self.assertIn(reference_cache_key(str(ref), 21, 4), path)
            self.assertTrue(path.endswith(".dict"))

    def test_reference_cache_path_preserves_legacy_output_when_disabled(self):
        args = SimpleNamespace(
            reuse_reference_cache=False,
            reference_cache_dir=None,
            r="/tmp/ref",
            kf=21,
            step_size=4,
        )

        self.assertEqual(
            get_reference_kmer_dict_path(args, "/tmp/out"),
            "/tmp/out/kmer_dict_k21.dict",
        )

    def test_assembler_cache_dir_is_nested_under_reference_cache(self):
        args = SimpleNamespace(
            reuse_reference_cache=True,
            reference_cache_dir="/tmp/cache",
        )

        self.assertEqual(
            get_assembler_reference_cache_dir(args, "/tmp/out"),
            "/tmp/cache/assembler",
        )

    def test_assembler_reference_cache_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "uce-1.fasta"
            cache = root / "cache"
            ref.write_text(">a\nACGTACGTACGT\n")

            first = {}
            second = {}

            self.assertFalse(Load_Or_Make_Kmer_Dict(first, str(ref), 5, str(cache)))
            self.assertTrue(Load_Or_Make_Kmer_Dict(second, str(ref), 5, str(cache)))
            self.assertEqual(first, second)
            self.assertEqual(len(list(cache.glob("*.pkl"))), 1)

    def test_reference_cache_directory_exists_before_first_filter_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "ref"
            out = root / "out"
            ref.mkdir()
            out.mkdir()
            (ref / "uce-1.fasta").write_text(">a\nACGTACGT\n")

            args = SimpleNamespace(
                o=str(out),
                r=str(ref),
                kf=21,
                step_size=4,
                reuse_reference_cache=True,
                reference_cache_dir=None,
                uce_rescue_reads=False,
                assembly_mode="original",
                p=1,
                soft_boundary="auto",
            )

            def fake_run(command, check):
                cache_path = Path(command[command.index("-lkd") + 1])
                self.assertTrue(cache_path.parent.is_dir())
                cache_path.write_text("cache")

            with patch("unix_command.find_executable", return_value="MainFilterNew"):
                with patch("unix_command.subprocess.run", side_effect=fake_run):
                    do_filter_assemble(args, {}, True, False, False)


if __name__ == "__main__":
    unittest.main()
