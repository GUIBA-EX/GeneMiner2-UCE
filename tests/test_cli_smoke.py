from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import unittest

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
        self.assertIn("--uce-rescue-reads", proc.stdout)
        self.assertIn("--uce-rescue-min-density-ratio", proc.stdout)
        self.assertIn("--uce-max-contig-length", proc.stdout)
        self.assertIn("--reuse-reference-cache", proc.stdout)

    def test_default_alifilter_model_uses_builtin_model(self):
        for value in (None, "", "default", " DEFAULT "):
            args = SimpleNamespace(alifilter_model=value)
            self.assertIsNone(unix_command.get_alifilter_model(args))

    def test_custom_alifilter_model_is_preserved(self):
        args = SimpleNamespace(alifilter_model="/path/to/model.json")
        self.assertEqual(unix_command.get_alifilter_model(args), "/path/to/model.json")

if __name__ == "__main__":
    unittest.main()
