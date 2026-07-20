import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "rust" / "build_consensus" / "Cargo.toml"
BINARY = ROOT / "rust" / "build_consensus" / "target" / "debug" / "build_consensus"


class RustConsensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("cargo") is None:
            raise unittest.SkipTest("cargo is not installed")
        subprocess.run(["cargo", "build", "--quiet", "--manifest-path", str(MANIFEST)], check=True, cwd=ROOT)

    def test_gzip_sam_preserves_insertions_and_ignores_non_primary_records(self):
        sam = (
            "@HD\tVN:1.6\n"
            "@SQ\tSN:ref1\tLN:4\n"
            "@SQ\tSN:ref2\tLN:4\n"
            "r1\t0\tref1\t1\t60\t2M1I2M\t*\t0\t0\tACTGT\tFFFFF\n"
            "r2\t0\tref2\t1\t60\t4M\t*\t0\t0\tCCCC\tFFFF\n"
            "secondary\t256\tref1\t1\t60\t4M\t*\t0\t0\tTTTT\tFFFF\n"
            "supplementary\t2048\tref1\t1\t60\t4M\t*\t0\t0\tTTTT\tFFFF\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            source = work / "input.sam.gz"
            with gzip.open(source, "wt") as out:
                out.write(sam)
            subprocess.run([str(BINARY), "-i", str(source), "-o", str(work), "-p", "consensus",
                            "-c", "0.75", "-m", "1", "-s", "0"], check=True, cwd=ROOT)
            output = (work / "consensus.fasta").read_text()
            self.assertIn("reference:ref1", output)
            self.assertIn("ACTGT", output)
            self.assertIn("reference:ref2", output)
            self.assertIn("CCCC", output)
            self.assertNotIn("TTTT", output)


    def test_gzip_default_prefix_matches_legacy_splitext(self):
        sam = "@SQ\tSN:ref\tLN:1\nread\t0\tref\t1\t60\t1M\t*\t0\t0\tA\tF\n"
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            source = work / "input.sam.gz"
            with gzip.open(source, "wt") as out:
                out.write(sam)
            subprocess.run([str(BINARY), "-i", str(source), "-o", str(work), "-m", "1", "-s", "0"], check=True, cwd=ROOT)
            output = work / "input.sam.fasta"
            self.assertTrue(output.is_file())
            self.assertTrue(output.read_text().startswith(">input.sam|c25 reference:ref"))

    def test_mutation_density_png_is_written_when_supported(self):
        sam = "@SQ\tSN:ref\tLN:3\n" + "".join(
            f"r{i}\t0\tref\t1\t60\t3M\t*\t0\t0\t{sequence}\tFFF\n"
            for i, sequence in enumerate(("AAA", "AAA", "AAA", "CCC", "CCC"), 1)
        )
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            source = work / "input.sam"
            source.write_text(sam)
            subprocess.run([str(BINARY), "-i", str(source), "-o", str(work), "-p", "density", "-m", "1", "-s", "1"], check=True, cwd=ROOT)
            image = work / "density.png"
            self.assertTrue(image.is_file())
            self.assertGreater(image.stat().st_size, 100)


if __name__ == "__main__":
    unittest.main()
