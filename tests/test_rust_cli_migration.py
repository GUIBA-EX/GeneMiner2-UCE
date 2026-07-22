from pathlib import Path
import csv
import os
import random
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "rust" / "geneminer2_cli" / "Cargo.toml"


class RustCliMigrationTests(unittest.TestCase):
    def test_rust_preview_help_is_available(self):
        proc = subprocess.run(
            ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "--help"],
            cwd=ROOT,
            text=True,
            env={**os.environ, "GENEMINER2_ENGINE": "rust"},
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Rust migration launcher", proc.stdout)

    def test_output_comparator_accepts_identical_directories(self):
        fixture = ROOT / "tests" / "__pycache__"
        proc = subprocess.run(
            ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--bin", "compare_outputs", "--", str(fixture), str(fixture)],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Compatibility check passed", proc.stdout)

    def test_native_uce_pipeline_recovers_a_paired_synthetic_locus(self):
        rng = random.Random(20260722)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'
            references.mkdir()
            (references / 'uce_demo.fasta').write_text(f'>uce_demo\n{truth}\n')
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n')
                    second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{r1}\t{r2}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "--assembly-mode", "uce", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out'), "-kf", "31", "-ka", "31", "--min-ka", "31", "--max-ka", "31"],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with (tmp / 'out' / '1_Demo' / 'uce_assembly_summary.csv').open(newline='') as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row['accepted'], '1')
            self.assertGreaterEqual(int(row['selected_contig_length']), 650)

    def test_native_gene_pipeline_runs_filter_to_cohort(self):
        rng = random.Random(77)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        complement = str.maketrans('ACGT', 'TGCA')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'
            references.mkdir()
            (references / 'gene_demo.fasta').write_text(f'>gene_demo\n{truth}\n')
            (tmp / 'out' / '99_Stale').mkdir(parents=True)
            r1, r2 = tmp / 'r1.fq', tmp / 'r2.fq'
            with r1.open('w') as first, r2.open('w') as second:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(f'@read{index}/1\n{read}\n+\n' + 'I' * 150 + '\n')
                    second.write(f'@read{index}/2\n{read.translate(complement)[::-1]}\n+\n' + 'I' * 150 + '\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{r1}\t{r2}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "gene", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out'), "-kf", "31", "-ka", "31", "--min-ka", "31", "--max-ka", "31"],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=120,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((tmp / 'out' / '1_Demo' / 'results' / 'gene_demo.fasta').is_file())
            self.assertTrue((tmp / 'out' / 'gene' / 'family_summary.tsv').is_file())
            self.assertNotIn('99_Stale', (tmp / 'out' / 'gene' / 'family_count_matrix.tsv').read_text())

    def test_native_gene_pipeline_supports_single_end_reads(self):
        rng = random.Random(78)
        truth = ''.join(rng.choice('ACGT') for _ in range(700))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            references = tmp / 'references'
            references.mkdir()
            (references / 'gene_demo.fasta').write_text(f'>gene_demo' + chr(10) + truth + chr(10))
            r1 = tmp / 'r1.fq'
            with r1.open('w') as first:
                for index, start in enumerate(range(0, 551, 10)):
                    read = truth[start:start + 150]
                    first.write(
                        f'@read{index}' + chr(10) + read + chr(10) + '+' + chr(10)
                        + 'I' * 150 + chr(10)
                    )
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo' + chr(9) + str(r1) + chr(10))
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "gene", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out'), "-kf", "31", "-ka", "31", "--min-ka", "31", "--max-ka", "31"],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=120,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((tmp / 'out' / '1_Demo' / 'results' / 'gene_demo.fasta').is_file())

    def test_native_uce_refilter_alone_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            reads, references = tmp / 'reads.fq', tmp / 'references'
            reads.write_text('@read/1\nACGTACGTACGTACGTACGTACGTACGTACG\n+\n' + 'I' * 32 + '\n')
            references.mkdir()
            (references / 'uce.fasta').write_text('>uce\nACGTACGTACGTACGTACGTACGTACGTACG\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{reads}\t{reads}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "refilter", "--assembly-mode", "uce", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out')],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('UCE refilter is fused', proc.stderr)

    def test_native_uce_rejects_single_end_before_launching_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            reads, references = tmp / 'reads.fq', tmp / 'references'
            reads.write_text('@read\nACGTACGTACGTACGTACGTACGTACGTACG\n+\n' + 'I' * 32 + '\n')
            references.mkdir()
            (references / 'uce.fasta').write_text('>uce\nACGTACGTACGTACGTACGTACGTACGTACG\n')
            sample_list = tmp / 'samples.tsv'
            sample_list.write_text(f'demo\t{reads}\n')
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(MANIFEST), "--", "--assembly-mode", "uce", "-f", str(sample_list), "-r", str(references), "-o", str(tmp / 'out')],
                cwd=ROOT,
                text=True,
                env={**os.environ, "GENEMINER2_ENGINE": "rust", "GM2_COMPONENT_DIR": str(ROOT / 'cli' / 'bin')},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('requires paired-end input', proc.stderr)


if __name__ == "__main__":
    unittest.main()
