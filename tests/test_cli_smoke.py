from pathlib import Path
import subprocess
import sys
import unittest


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

if __name__ == "__main__":
    unittest.main()
