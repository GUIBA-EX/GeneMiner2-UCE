import csv
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "rust" / "main_assembler" / "Cargo.toml"
BINARY = ROOT / "rust" / "main_assembler" / "target" / "debug" / "main_assembler"


def random_dna(length, seed=8128):
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def write_input(root, truth):
    references = root / "references"
    output = root / "output"
    filtered = output / "filtered"
    references.mkdir(parents=True)
    filtered.mkdir(parents=True)
    (references / "uce_test.fasta").write_text(f">uce_test\n{truth}\n")

    with (filtered / "uce_test.fq").open("w") as handle:
        read_index = 0
        for start in range(0, len(truth) - 150 + 1, 15):
            read = truth[start : start + 150]
            for copy in range(3):
                handle.write(
                    f"@read_{read_index}_{copy}\n{read}\n+\n{'I' * len(read)}\n"
                )
            read_index += 1
    return references, output


def selected_sequence(output):
    lines = (output / "results" / "uce_test.fasta").read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">"))


class RustMainAssemblerParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("cargo") is None:
            raise unittest.SkipTest("cargo is not installed")
        subprocess.run(
            ["cargo", "build", "--quiet", "--manifest-path", str(MANIFEST)],
            check=True,
            cwd=ROOT,
        )

    def test_rust_matches_python_backbone_on_synthetic_uce(self):
        truth = random_dna(600)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            py_reference, py_output = write_input(work / "python", truth)
            rs_reference, rs_output = write_input(work / "rust", truth)
            common = [
                "-ka",
                "21",
                "-limit_count",
                "1",
                "-iteration",
                "1024",
                "-p",
                "1",
                "--assembly-mode",
                "uce",
                "--uce-path-strategy",
                "backbone",
                "--uce-backbone-lookahead",
                "24",
                "--uce-max-contig-length",
                "0",
            ]

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "main_assembler.py"),
                    "-r",
                    str(py_reference),
                    "-o",
                    str(py_output),
                    *common,
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )
            subprocess.run(
                [
                    str(BINARY),
                    "-r",
                    str(rs_reference),
                    "-o",
                    str(rs_output),
                    *common,
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )

            self.assertEqual(selected_sequence(rs_output), selected_sequence(py_output))
            with (py_output / "uce_assembly_summary.csv").open(newline="") as handle:
                py_row = next(csv.DictReader(handle))
            with (rs_output / "uce_assembly_summary.csv").open(newline="") as handle:
                rs_row = next(csv.DictReader(handle))
            for field in (
                "status",
                "accepted",
                "selected_contig_length",
                "slice_supported_bases",
                "unique_read_count",
                "candidate_count",
            ):
                self.assertEqual(rs_row[field], py_row[field], field)

    def test_rust_matches_python_reference_mode(self):
        truth = random_dna(500, seed=991)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            py_reference, py_output = write_input(work / "python", truth)
            rs_reference, rs_output = write_input(work / "rust", truth)
            common = [
                "-ka",
                "21",
                "-limit_count",
                "1",
                "-iteration",
                "1024",
                "-p",
                "1",
                "--assembly-mode",
                "original",
            ]

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "main_assembler.py"),
                    "-r",
                    str(py_reference),
                    "-o",
                    str(py_output),
                    *common,
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )
            subprocess.run(
                [
                    str(BINARY),
                    "-r",
                    str(rs_reference),
                    "-o",
                    str(rs_output),
                    *common,
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )

            rs_sequence = selected_sequence(rs_output)
            py_sequence = selected_sequence(py_output)
            self.assertEqual(len(rs_sequence), len(py_sequence))
            self.assertIn(rs_sequence, truth)
            self.assertIn(py_sequence, truth)
            self.assertEqual(
                (rs_output / "result_dict.txt").read_text(),
                (py_output / "result_dict.txt").read_text(),
            )

if __name__ == "__main__":
    unittest.main()
