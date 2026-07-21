from pathlib import Path
import csv
import random
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import unix_command


class UceTerminalRescueTests(unittest.TestCase):
    @staticmethod
    def dna(length, seed):
        rng = random.Random(seed)
        return "".join(rng.choice("ACGT") for _ in range(length))

    def test_terminal_reconcile_keeps_supported_side_and_preserves_candidates(self):
        left = self.dna(40, 1)
        core = self.dna(120, 2)
        right = self.dna(40, 3)
        assembled = left + core + right

        with tempfile.TemporaryDirectory() as tmp:
            sample = Path(tmp) / "sample"
            backup = Path(tmp) / "backup"
            for root in (sample, backup):
                (root / "results").mkdir(parents=True)
            (sample / "filtered").mkdir()
            (sample / "contigs_all").mkdir()
            (sample / "contigs_all_low").mkdir()
            (backup / "results" / "locus.fasta").write_text(f">old\n{core}\n")
            (sample / "results" / "locus.fasta").write_text(f">new\n{assembled}\n")
            candidate = ">candidate_1\n" + assembled + "\n>candidate_2\nACGT\n"
            (sample / "contigs_all" / "locus.fasta").write_text(candidate)
            (sample / "contigs_all_low" / "locus.fasta").write_text(candidate)

            spanning = assembled[:100]
            with (sample / "filtered" / "locus.fq").open("w") as out:
                for fragment in ("frag1", "frag2"):
                    out.write(f"@{fragment}/1\n{spanning}\n+\n" + "I" * len(spanning) + "\n")

            after = {
                "selected_contig_length": str(len(assembled)),
                "unique_read_count": "2",
            }
            evidence, status = unix_command.terminal_reconcile_locus(
                str(sample), str(backup), "locus", after
            )

            self.assertEqual(status, "accepted")
            self.assertTrue(evidence["left"]["accepted"])
            self.assertFalse(evidence["right"]["accepted"])
            selected = unix_command.read_first_fasta_sequence(sample / "results" / "locus.fasta")
            self.assertEqual(selected, left + core)
            self.assertEqual((sample / "contigs_all" / "locus.fasta").read_text(), candidate)
            self.assertEqual((sample / "contigs_all_low" / "locus.fasta").read_text(), candidate)

    def test_round_audit_contains_terminal_evidence(self):
        records = []
        evidence = {
            "locus": {
                "left": {"length": 40, "breadth": 0.9, "max_gap": 2, "fragments": 3, "bridges": 2, "accepted": True},
                "right": {"length": 35, "breadth": 0.4, "max_gap": 20, "fragments": 1, "bridges": 0, "accepted": False},
            }
        }
        before = {"locus": {"status": "success", "selected_contig_length": "100", "unique_read_count": "2"}}
        after = {"locus": {"status": "success", "selected_contig_length": "140", "unique_read_count": "4"}}
        unix_command.append_uce_rescue_round_records(
            records, "sample", 2, before, after,
            {"locus": "terminal_left_kept_right_trimmed"},
            terminal_evidence=evidence,
        )
        self.assertEqual(records[0]["left_bridges"], 2)
        self.assertEqual(records[0]["left_accepted"], 1)
        self.assertEqual(records[0]["right_accepted"], 0)

    def test_only_growing_loci_enter_terminal_round(self):
        before = {
            "growing": {"accepted": "1", "selected_contig_length": "100", "unique_read_count": "2"},
            "stable": {"accepted": "1", "selected_contig_length": "100", "unique_read_count": "2"},
        }
        after = {
            "growing": {"accepted": "1", "selected_contig_length": "140", "unique_read_count": "3"},
            "stable": {"accepted": "1", "selected_contig_length": "100", "unique_read_count": "2"},
            "new": {"accepted": "1", "selected_contig_length": "80", "unique_read_count": "2"},
        }
        self.assertEqual(
            unix_command.select_terminal_rescue_loci(before, after),
            {"growing", "new"},
        )


if __name__ == "__main__":
    unittest.main()
