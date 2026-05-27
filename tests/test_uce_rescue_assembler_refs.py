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


if __name__ == "__main__":
    unittest.main()
