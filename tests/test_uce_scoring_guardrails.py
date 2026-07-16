from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from main_assembler import Evaluate_Uce_Candidate, Format_Contig_Header, Pass_Uce_Guardrails, Score_Contig


def guard_args(**overrides):
    values = {
        "max_contig_length": 5000,
        "min_read_density": 0.003,
        "density_check_min_length": 1000,
        "max_depth_cv": 0,
        "max_depth_ratio": 0,
    }
    values.update(overrides)
    return values


class UceScoringGuardrailTests(unittest.TestCase):
    def test_reference_header_remains_backward_compatible(self):
        contig = [None] * 19
        contig[0] = 'AAAA'
        contig[1] = 1
        contig[2] = 2
        contig[3] = 3
        contig[4] = 4
        contig[5] = 4
        contig[6] = 1.0
        contig[14] = 4

        self.assertEqual(
            Format_Contig_Header(contig, 'reference'),
            '>contig_4_1_2_3_4_span_4_balance_1.000',
        )
        self.assertIn('_supported_4_', Format_Contig_Header(contig, 'uce'))

    def test_long_low_density_contig_is_rejected(self):
        self.assertFalse(
            Pass_Uce_Guardrails(
                contig_len=3000,
                read_density=0.001,
                depth_cv=0.2,
                max_depth_ratio=2,
                guardrails=guard_args(),
            )
        )

    def test_short_low_density_contig_is_not_rejected_by_density_guardrail(self):
        self.assertTrue(
            Pass_Uce_Guardrails(
                contig_len=500,
                read_density=0.001,
                depth_cv=0.2,
                max_depth_ratio=2,
                guardrails=guard_args(),
            )
        )

    def test_uce_score_penalizes_weak_read_density(self):
        weak_long = [
            "A" * 3000, 1, 0, 10, 3, 2900, 0.5, 0.001, 0.967,
            2, 0.1, 1.0, 3, 0, 300, 0.1, 2700, True, "",
        ]
        dense_short = [
            "A" * 800, 1, 0, 10, 8, 700, 0.5, 0.01, 0.875,
            2, 0.1, 1.0, 8, 0, 700, 0.875, 100, True, "",
        ]

        self.assertGreater(Score_Contig(dense_short, "uce"), Score_Contig(weak_long, "uce"))

    def test_uce_acceptance_requires_unique_positional_support(self):
        accepted, reasons = Evaluate_Uce_Candidate(
            contig_len=800,
            unique_read_count=0,
            supported_bases=0,
            unique_read_density=0,
            depth_cv=0.2,
            max_depth_ratio=2,
            guardrails=guard_args(),
        )

        self.assertFalse(accepted)
        self.assertEqual(reasons, ['no_unique_read_support', 'no_positional_support'])

    def test_uce_acceptance_uses_unique_read_density(self):
        accepted, reasons = Evaluate_Uce_Candidate(
            contig_len=3000,
            unique_read_count=3,
            supported_bases=300,
            unique_read_density=0.001,
            depth_cv=0.2,
            max_depth_ratio=2,
            guardrails=guard_args(),
        )

        self.assertFalse(accepted)
        self.assertIn('low_unique_read_density', reasons)


if __name__ == "__main__":
    unittest.main()
