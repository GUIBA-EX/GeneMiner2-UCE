from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TOOLS_MANIFEST = ROOT / "rust" / "gm2_tools" / "Cargo.toml"
CLI_MANIFEST = ROOT / "rust" / "geneminer2_cli" / "Cargo.toml"


class RadWorkflowTests(unittest.TestCase):
    @staticmethod
    def write_loci(path):
        path.write_text(
            "radA.trimmed x ACGTACGTACGTACGTACGTNNNGGGGTTTTGGGGTTTTGGGG\n"
            "radB.trimmed x ACGTACGTACGTACGTACGTNNNGGGGTTTTGGGGTTTTGGGG\n"
            "/ / |42| /\n"
        )

    def run_tool(self, *args):
        return subprocess.run(
            ["cargo", "run", "--quiet", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow", "--", *args],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60,
        )

    def test_reference_and_finalize_keep_arms_independent(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            loci, reference, recovery, matrix = tmp / "input.loci", tmp / "reference", tmp / "recovery", tmp / "matrix"
            self.write_loci(loci)
            proc = self.run_tool("reference", "--loci", str(loci), "--out", str(reference))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((reference / "arms" / "42__R1.fasta").is_file())
            self.assertTrue((reference / "arms" / "42__R2.fasta").is_file())
            result = recovery / "1_New" / "results"; result.mkdir(parents=True)
            (result / "42__R1.fasta").write_text(">new\nACGTACGTACGTACGTACGT\n")
            (result / "42__R2.fasta").write_text(">new\nGGGGTTTTGGGGTTTTGGGG\n")
            proc = self.run_tool(
                "finalize", "--reference", str(reference), "--recovery", str(recovery), "--out", str(matrix),
                "--sample", "1_New", "--min-arm-breadth", "0.8",
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = (matrix / "rad_sample_locus.tsv").read_text()
            self.assertIn("rad_missing_wgs_recovered", report)
            strict = (matrix / "paired_arms" / "R1" / "42.fasta").read_text()
            self.assertIn(">1_New", strict)
            self.assertIn("independent observations", (matrix / "README.txt").read_text())

    def test_partial_arm_is_not_in_strict_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            loci, reference, recovery, matrix = tmp / "input.loci", tmp / "reference", tmp / "recovery", tmp / "matrix"
            self.write_loci(loci)
            self.assertEqual(self.run_tool("reference", "--loci", str(loci), "--out", str(reference)).returncode, 0)
            result = recovery / "1_New" / "results"; result.mkdir(parents=True)
            (result / "42__R1.fasta").write_text(">new\nACGTACGTACGTACGTACGT\n")
            proc = self.run_tool("finalize", "--reference", str(reference), "--recovery", str(recovery), "--out", str(matrix), "--sample", "1_New")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("partial_arm_recovery", (matrix / "rad_sample_locus.tsv").read_text())
            self.assertIn(">1_New", (matrix / "recovered_arms" / "R1" / "42.fasta").read_text())
            self.assertNotIn(">1_New", (matrix / "paired_arms" / "R1" / "42.fasta").read_text())

    def test_validate_requires_unique_full_length_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, output = tmp / "reference", tmp / "recovery", tmp / "validated"
            arms = reference / "arms"; arms.mkdir(parents=True)
            seq1, seq2 = "ACGTACGTACGTACGTACGT", "TTTTGGGGTTTTGGGGTTTT"
            (arms / "A__R1.fasta").write_text(f">radA\n{seq1}\n")
            (arms / "A__R2.fasta").write_text(f">radA\n{seq2}\n")
            # A distinct locus is required to exercise the cross-locus competitor.
            (arms / "B__R1.fasta").write_text(">radB\nCCCCAAAACCCCAAAACCCC\n")
            (arms / "B__R2.fasta").write_text(">radB\nAAAACCCCAAAACCCCAAAA\n")
            results = recovery / "new" / "results"; results.mkdir(parents=True)
            (results / "A__R1.fasta").write_text(f">candidate\n{seq1}\n")
            (results / "A__R2.fasta").write_text(f">candidate\n{seq2}\n")
            proc = self.run_tool("validate", "--reference", str(reference), "--recovery", str(recovery), "--out", str(output), "--min-identity", "0.9", "--min-breadth", "0.8", "--min-delta", "0.05")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = (output / "rad_validation.tsv").read_text()
            self.assertIn("new\tA\tR1\tvalidated", report)
            self.assertIn(">new", (output / "strict_arms" / "R1" / "A.fasta").read_text())

    def test_validate_rejects_ambiguous_and_short_arms(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, output = tmp / "reference", tmp / "recovery", tmp / "validated"
            arms = reference / "arms"; arms.mkdir(parents=True)
            shared, mate = "ACGTACGTACGTACGTACGT", "TTTTGGGGTTTTGGGGTTTT"
            # A/B deliberately share R1 so A's R1 has no cross-locus score margin.
            for locus in ("A", "B"):
                (arms / f"{locus}__R1.fasta").write_text(f">rad{locus}\n{shared}\n")
                (arms / f"{locus}__R2.fasta").write_text(f">rad{locus}\n{mate}\n")
            results = recovery / "new" / "results"; results.mkdir(parents=True)
            (results / "A__R1.fasta").write_text(f">candidate\n{shared}\n")
            (results / "A__R2.fasta").write_text(">candidate\nTTTTGGGG\n")
            proc = self.run_tool("validate", "--reference", str(reference), "--recovery", str(recovery), "--out", str(output))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = (output / "rad_validation.tsv").read_text()
            self.assertIn("new\tA\tR1\tambiguous_paralog", report)
            self.assertIn("new\tA\tR2\tinsufficient_coverage", report)
            self.assertNotIn(">new", (output / "strict_arms" / "R1" / "A.fasta").read_text())

    def test_denovo_builds_paired_candidate_arms(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); out = tmp / "reference"
            r1a, r2a, r1b, r2b = [tmp / name for name in ("a_R1.fq", "a_R2.fq", "b_R1.fq", "b_R2.fq")]
            def fq(sequence):
                return "".join(f"@r{i}\nTGCAG{sequence}\n+\n{'I' * (len(sequence) + 5)}\n" for i in range(3))
            r1a.write_text(fq("ACGTACGTACGTACGTACGT")); r2a.write_text(fq("GGGGTTTTGGGGTTTTGGGG"))
            r1b.write_text(fq("ACGTACGTACGTACGTACGT")); r2b.write_text(fq("GGGGTTTTGGGGTTTTGGGG"))
            proc = self.run_tool(
                "denovo", "--out", str(out),
                "--sample", "radA", "--read1", str(r1a), "--read2", str(r2a),
                "--sample", "radB", "--read1", str(r1b), "--read2", str(r2b),
                "--overhang", "TGCAG", "--overhang-r2", "TGCAG", "--kmer", "15", "--min-count", "3", "--min-samples", "2", "--min-length", "20",
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            manifest = (out / "locus_manifest.tsv").read_text()
            self.assertIn("denovo_candidate", manifest)
            arms = sorted((out / "arms").glob("*.fasta"))
            self.assertEqual(len(arms), 2)
            self.assertIn(">radA", arms[0].read_text())
            self.assertIn("pass", (out / "denovo_probe_evidence.tsv").read_text())

    def test_rad_probe_accepts_existing_loci_without_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / "components"; components.mkdir()
            tool = ROOT / "rust" / "gm2_tools" / "target" / "debug" / "rad_workflow"
            build = subprocess.run(["cargo", "build", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertEqual(build.returncode, 0, build.stderr)
            (components / "rad_workflow").symlink_to(tool)
            loci = tmp / "input.loci"; self.write_loci(loci)
            output = tmp / "probe"
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(CLI_MANIFEST), "--", "rad-probe", "--ipyrad-loci", str(loci), "-o", str(output)],
                cwd=ROOT, text=True, env={**os.environ, "GM2_COMPONENT_DIR": str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / "rad_reference" / "arms" / "42__R1.fasta").is_file())
            self.assertIn("existing ipyrad .loci", (output / "rad_reference" / "PROVENANCE.txt").read_text())

    def test_rad_probe_runs_ipyrad_params_and_discovers_loci(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / "components"; components.mkdir()
            tool = ROOT / "rust" / "gm2_tools" / "target" / "debug" / "rad_workflow"
            build = subprocess.run(["cargo", "build", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertEqual(build.returncode, 0, build.stderr)
            (components / "rad_workflow").symlink_to(tool)
            project = tmp / "project"; params = tmp / "params.txt"
            params.write_text(f"demo ## [0] [assembly_name]\n{project} ## [1] [project_dir]\n")
            fake = tmp / "fake_ipyrad"
            fake.write_text("#!/bin/sh\nmkdir -p '" + str(project / "demo_outfiles") + "'\nprintf 'radA x ACGTACGTACGTACGTNNNGGGGTTTTGGGGTTTT\n/ / |42| /\n' > '" + str(project / "demo_outfiles" / "demo.loci") + "'\n")
            fake.chmod(0o755)
            output = tmp / "probe"
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(CLI_MANIFEST), "--", "rad-probe", "--ipyrad-params", str(params), "--ipyrad-executable", str(fake), "-o", str(output)],
                cwd=ROOT, text=True, env={**os.environ, "GM2_COMPONENT_DIR": str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / "rad_reference" / "arms" / "42__R2.fasta").is_file())
            self.assertIn("steps=1234567", (output / "rad_reference" / "PROVENANCE.txt").read_text())

    def test_cli_rad_validate_dispatches(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / "components"; components.mkdir()
            tool = ROOT / "rust" / "gm2_tools" / "target" / "debug" / "rad_workflow"
            build = subprocess.run(["cargo", "build", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertEqual(build.returncode, 0, build.stderr)
            (components / "rad_workflow").symlink_to(tool)
            arms = tmp / "reference" / "arms"; arms.mkdir(parents=True)
            (arms / "A__R1.fasta").write_text(">rad\nACGTACGTACGTACGTACGT\n")
            (arms / "A__R2.fasta").write_text(">rad\nTTTTGGGGTTTTGGGGTTTT\n")
            results = tmp / "recovery" / "new" / "results"; results.mkdir(parents=True)
            (results / "A__R1.fasta").write_text(">new\nACGTACGTACGTACGTACGT\n")
            (results / "A__R2.fasta").write_text(">new\nTTTTGGGGTTTTGGGGTTTT\n")
            output = tmp / "out"
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(CLI_MANIFEST), "--", "rad-validate", "--rad-probe", str(tmp / "reference"), "--rad-recovery", str(tmp / "recovery"), "-o", str(output)],
                cwd=ROOT, text=True, env={**os.environ, "GM2_COMPONENT_DIR": str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((output / "rad_validated" / "rad_validation.tsv").is_file())

    def test_cli_rad_probe_denovo_uses_paired_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / "components"; components.mkdir()
            tool = ROOT / "rust" / "gm2_tools" / "target" / "debug" / "rad_workflow"
            build = subprocess.run(["cargo", "build", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertEqual(build.returncode, 0, build.stderr)
            (components / "rad_workflow").symlink_to(tool)
            reads = []
            for sample in ("a", "b"):
                r1, r2 = tmp / f"{sample}_R1.fq", tmp / f"{sample}_R2.fq"
                r1.write_text("".join(f"@r{i}\nTGCAGACGTACGTACGTACGTACGT\n+\n{'I' * 25}\n" for i in range(3)))
                r2.write_text("".join(f"@r{i}\nGGGGTTTTGGGGTTTTGGGG\n+\n{'I' * 20}\n" for i in range(3)))
                reads.append((sample, r1, r2))
            samples = tmp / "rad.tsv"; samples.write_text("".join(f"{sample}\t{r1}\t{r2}\n" for sample, r1, r2 in reads))
            output = tmp / "probe"
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(CLI_MANIFEST), "--", "rad-probe", "--rad-denovo", "-f", str(samples), "-o", str(output), "--rad-overhang", "TGCAG", "--rad-kmer", "15", "--rad-min-count", "3", "--rad-min-samples", "2", "--rad-min-length", "20"],
                cwd=ROOT, text=True, env={**os.environ, "GM2_COMPONENT_DIR": str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("denovo_candidate_probe", (output / "rad_reference" / "PROVENANCE.txt").read_text())
            self.assertTrue(list((output / "rad_reference" / "arms").glob("*__R1.fasta")))

    def test_cli_runs_reference_recovery_and_finalize(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); components = tmp / "components"; components.mkdir()
            tool = ROOT / "rust" / "gm2_tools" / "target" / "debug" / "rad_workflow"
            build = subprocess.run(["cargo", "build", "--manifest-path", str(TOOLS_MANIFEST), "--bin", "rad_workflow"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertEqual(build.returncode, 0, build.stderr)
            (components / "rad_workflow").symlink_to(tool)
            (components / "MainFilterNew").write_text("#!/bin/sh\nexit 0\n")
            (components / "main_refilter_new").write_text("#!/bin/sh\nwhile [ \"$#\" -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then mkdir -p \"$2\"; shift 2; else shift; fi; done\n")
            (components / "main_assembler-original-rust").write_text("#!/bin/sh\nwhile [ \"$#\" -gt 0 ]; do case \"$1\" in -r) ref=\"$2\"; shift 2;; -o) out=\"$2\"; shift 2;; *) shift;; esac; done\nmkdir -p \"$out/results\"\nfor file in \"$ref\"/*.fasta; do cp \"$file\" \"$out/results/$(basename \"$file\")\"; done\n")
            for name in ("MainFilterNew", "main_refilter_new", "main_assembler-original-rust"):
                (components / name).chmod(0o755)
            loci = tmp / "input.loci"; self.write_loci(loci)
            reads = tmp / "reads.fq"; reads.write_text("@r\nACGT\n+\nIIII\n")
            samples = tmp / "samples.tsv"; samples.write_text(f"new\t{reads}\n")
            output = tmp / "out"
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", str(CLI_MANIFEST), "--", "rad", "--ipyrad-loci", str(loci), "-f", str(samples), "-o", str(output), "-p", "1"],
                cwd=ROOT, text=True, env={**os.environ, "GM2_COMPONENT_DIR": str(components)}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("rad_missing_wgs_recovered", (output / "rad_matrix" / "rad_sample_locus.tsv").read_text())

    def test_validate_keeps_the_exact_selected_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, output = tmp / "reference", tmp / "recovery", tmp / "validated"
            arms = reference / "arms"; arms.mkdir(parents=True)
            r1, r2 = "ACGTACGTACGTACGTACGT", "TTTTGGGGTTTTGGGGTTTT"
            (arms / "A__R1.fasta").write_text(f">radA\n{r1}\n")
            (arms / "A__R2.fasta").write_text(f">radA\n{r2}\n")
            (arms / "B__R1.fasta").write_text(">radB\nCCCCAAAACCCCAAAACCCC\n")
            (arms / "B__R2.fasta").write_text(">radB\nAAAACCCCAAAACCCCAAAA\n")
            results = recovery / "new" / "results"; results.mkdir(parents=True)
            # Duplicate headers occur in some candidate FASTA files.  The second
            # record is the one selected by alignment and must be the one emitted.
            (results / "A__R1.fasta").write_text(f">candidate\nCCCCCCCCCCCCCCCCCCCC\n>candidate\n{r1}\n")
            (results / "A__R2.fasta").write_text(f">candidate\n{r2}\n")
            proc = self.run_tool("validate", "--reference", str(reference), "--recovery", str(recovery), "--out", str(output))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            strict = (output / "strict_arms" / "R1" / "A.fasta").read_text()
            self.assertIn(f">new\n{r1}", strict)
            self.assertNotIn(">new\nCCCCCCCCCCCCCCCCCCCC", strict)
    def test_reference_rejects_normalized_sample_collision_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); loci, reference = tmp / "input.loci", tmp / "reference"
            loci.write_text(
                "a/b.trimmed x ACGTACGTACGTACGTACGTNNNGGGGTTTTGGGGTTTTGGGG\n"
                "a_b.trimmed x ACGTACGTACGTACGTACGTNNNGGGGTTTTGGGGTTTTGGGG\n"
                "/ / |42| /\n"
            )
            proc = self.run_tool("reference", "--loci", str(loci), "--out", str(reference))
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("safe-name normalization", proc.stderr)
            self.assertFalse(reference.exists())

    def test_validate_rejects_malformed_candidate_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, output = tmp / "reference", tmp / "recovery", tmp / "validated"
            arms = reference / "arms"; arms.mkdir(parents=True)
            r1, r2 = "ACGTACGTACGTACGTACGT", "TTTTGGGGTTTTGGGGTTTT"
            (arms / "A__R1.fasta").write_text(f">radA\n{r1}\n")
            (arms / "A__R2.fasta").write_text(f">radA\n{r2}\n")
            results = recovery / "new" / "results"; results.mkdir(parents=True)
            (results / "A__R1.fasta").write_text(r1 + "\n")
            (results / "A__R2.fasta").write_text(f">candidate\n{r2}\n")
            proc = self.run_tool("validate", "--reference", str(reference), "--recovery", str(recovery), "--out", str(output))
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("before a header", proc.stderr)
    def test_validate_rejects_incomplete_reference_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, output = tmp / "reference", tmp / "recovery", tmp / "validated"
            arms = reference / "arms"; arms.mkdir(parents=True)
            (arms / "A__R1.fasta").write_text(">radA\nACGTACGTACGTACGTACGT\n")
            (recovery / "new" / "results").mkdir(parents=True)
            proc = self.run_tool("validate", "--reference", str(reference), "--recovery", str(recovery), "--out", str(output))
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("both R1 and R2", proc.stderr)
            self.assertFalse(output.exists())

    def test_finalize_rejects_recovery_sample_already_in_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); reference, recovery, matrix = tmp / "reference", tmp / "recovery", tmp / "matrix"
            arms = reference / "arms"; arms.mkdir(parents=True)
            r1, r2 = "ACGTACGTACGTACGTACGT", "TTTTGGGGTTTTGGGGTTTT"
            (arms / "A__R1.fasta").write_text(f">radA\n{r1}\n")
            (arms / "A__R2.fasta").write_text(f">radA\n{r2}\n")
            (recovery / "radA" / "results").mkdir(parents=True)
            proc = self.run_tool("finalize", "--reference", str(reference), "--recovery", str(recovery), "--out", str(matrix), "--sample", "radA")
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("already exists in the RAD reference", proc.stderr)
            self.assertFalse(matrix.exists())

if __name__ == "__main__":
    unittest.main()
