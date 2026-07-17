from pathlib import Path
import csv
import random
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
POPULATION_BIN = ROOT / "cli" / "bin" / "main_population"


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def alternative_alleles(reference_base):
    return [base for base in "ACGT" if base != reference_base][:2]


def write_fastq_pair(read1_path, read2_path, loci, sample_index):
    with read1_path.open("w") as read1, read2_path.open("w") as read2:
        read_number = 0
        for locus, reference in loci.items():
            reference_sequence = reference
            sequence = list(reference)
            variant_position = {"uce1": 300, "uce2": 400, "uce3": 500}[locus]
            alternate1, alternate2 = alternative_alleles(reference[variant_position - 1])
            if locus in {"uce1", "uce2"} and sample_index >= 2:
                sequence[variant_position - 1] = alternate1
            elif locus == "uce3" and sample_index in {1, 2}:
                sequence[variant_position - 1] = alternate1
            elif locus == "uce3" and sample_index == 3:
                sequence[variant_position - 1] = alternate2
            sequence = "".join(sequence)
            fragment_starts = (
                [200]
                if locus == "uce1" and sample_index == 1
                else range(0, len(sequence) - 300 + 1, 25)
            )
            for fragment_index, fragment_start in enumerate(fragment_starts):
                # Sample C carries both alleles.  Mixing read pairs from the
                # reference and alternate haplotypes exercises heterozygous
                # genotype calling and supplies ADMIXTURE with a diploid panel.
                fragment_source = (
                    reference_sequence
                    if sample_index == 2 and fragment_index % 2 == 0
                    else sequence
                )
                fragment = fragment_source[fragment_start:fragment_start + 300]
                r1 = fragment[:150]
                r2 = reverse_complement(fragment[-150:])
                read_number += 1
                name = f"{locus}_{sample_index}_{read_number}"
                read1.write(f"@{name}/1\n{r1}\n+\n{'I' * len(r1)}\n")
                read2.write(f"@{name}/2\n{r2}\n+\n{'I' * len(r2)}\n")


@unittest.skipUnless(
    POPULATION_BIN.is_file()
    and all(shutil.which(tool) for tool in ("minibwa", "samtools", "bcftools", "plink")),
    "real population tools or compiled main_population are unavailable",
)
class PopulationRealToolSmokeTest(unittest.TestCase):
    def test_real_minibwa_bcftools_plink_pipeline(self):
        rng = random.Random(70195)
        loci = {
            "uce1": "".join(rng.choice("ACGT") for _ in range(800)),
            "uce2": "".join(rng.choice("ACGT") for _ in range(800)),
            "uce3": "".join(rng.choice("ACGT") for _ in range(800)),
        }
        summary_fields = [
            "locus", "status", "accepted", "selected_contig_length",
            "slice_supported_bases", "slice_support_breadth",
            "max_slice_support_gap", "unique_read_count",
            "unique_read_density", "low_quality",
        ]

        with tempfile.TemporaryDirectory(prefix="gm2-population-real-") as temporary:
            root = Path(temporary)
            output = root / "output"
            sample_rows = []
            for sample_index, sample_name in enumerate(("A", "B", "C", "D")):
                read1 = root / f"{sample_name}_R1.fastq"
                read2 = root / f"{sample_name}_R2.fastq"
                write_fastq_pair(read1, read2, loci, sample_index)
                sample_rows.append((sample_name, read1, read2))
                sample_dir = output / f"{sample_index + 1}_{sample_name}"
                results_dir = sample_dir / "results"
                results_dir.mkdir(parents=True)
                with (sample_dir / "uce_assembly_summary.csv").open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=summary_fields)
                    writer.writeheader()
                    for locus in loci:
                        writer.writerow({
                            "locus": locus,
                            "status": "success",
                            "accepted": "true",
                            "selected_contig_length": 800,
                            "slice_supported_bases": 800 if sample_index == 0 else 790,
                            "slice_support_breadth": 1.0,
                            "max_slice_support_gap": 0,
                            "unique_read_count": 42,
                            "unique_read_density": 0.0525,
                            "low_quality": "false",
                        })
                        (results_dir / f"{locus}.fasta").write_text(f">{sample_name}_{locus}\n{loci[locus]}\n")

            samples_tsv = root / "samples.tsv"
            with samples_tsv.open("w") as handle:
                for sample_name, read1, read2 in sample_rows:
                    handle.write(f"{sample_name}\t{read1}\t{read2}\n")

            completed = subprocess.run(
                [
                    str(POPULATION_BIN),
                    "--output", str(output),
                    "--samples", str(samples_tsv),
                    "--threads", "2",
                    "--min-call-rate", "0.75",
                    "--min-mac", "2",
                    "--skip-admixture",
                ],
                cwd=ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

            final_vcf = output / "population" / "structure" / "one_snp_per_uce.vcf.gz"
            query = subprocess.run(
                ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n", str(final_vcf)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            ).stdout.strip().splitlines()
            self.assertEqual(len(query), 3)
            records = {line.split("\t")[0]: line.split("\t") for line in query}
            self.assertEqual(set(records), {"uce1", "uce2", "uce3"})
            self.assertEqual(records["uce1"][2:], ["0/0", "./.", "0/1", "1/1"])
            self.assertEqual(records["uce2"][2:], ["0/0", "0/0", "0/1", "1/1"])
            filtered_vcf = output / "population" / "variants" / "cohort.filtered.vcf.gz"
            filtered_records = subprocess.run(
                ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(filtered_vcf)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            ).stdout.strip().splitlines()
            uce3_records = [line for line in filtered_records if line.startswith("uce3\t500\t")]
            self.assertEqual(len(uce3_records), 2)
            sample_ids = subprocess.run(
                ["bcftools", "query", "-l", str(final_vcf)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            ).stdout.strip().splitlines()
            self.assertEqual(sample_ids, ["A", "B", "C", "D"])
            self.assertTrue((output / "population" / "structure" / "population.bed").is_file())
            self.assertTrue((output / "population" / "structure" / "population_pca.eigenvec").is_file())
            self.assertTrue((output / "population" / "structure" / "all_snps.vcf.gz").is_file())
            self.assertTrue((output / "population" / "structure" / "all_snps_pca.eigenvec").is_file())
            self.assertTrue((output / "population" / "structure" / "ld_pruned.vcf.gz").is_file())
            self.assertTrue((output / "population" / "structure" / "ld_pruned_pca.eigenvec").is_file())
            self.assertTrue((output / "population" / "structure" / "panel_summary.tsv").is_file())
            self.assertTrue((output / "population" / "reference" / "reference_contribution.tsv").is_file())
            mapping_qc = (output / "population" / "mapping" / "mapping_qc.tsv").read_text()
            self.assertIn("mapping_rate", mapping_qc)
            self.assertIn("coverage_breadth", mapping_qc)


@unittest.skipUnless(
    all(shutil.which(tool) for tool in ("plink", "admixture")),
    "PLINK or ADMIXTURE is unavailable",
)
class AdmxtureRealToolSmokeTest(unittest.TestCase):
    def test_real_admixture_accepts_plink_binary_panel(self):
        """Confirm that the installed ADMIXTURE binary can fit a PLINK panel."""
        with tempfile.TemporaryDirectory(prefix="gm2-admixture-real-") as temporary:
            workdir = Path(temporary)
            prefix = workdir / "population"
            subprocess.run(
                [
                    "plink", "--dummy", "8", "40", "--make-bed",
                    "--out", str(prefix),
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            completed = subprocess.run(
                ["admixture", "--cv=2", "-j2", "population.bed", "2"],
                cwd=workdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue((workdir / "population.2.Q").is_file())
            self.assertTrue((workdir / "population.2.P").is_file())
            self.assertIn("CV error (K=2)", completed.stdout + completed.stderr)


if __name__ == "__main__":
    unittest.main()
