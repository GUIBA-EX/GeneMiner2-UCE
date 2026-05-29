from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from main_assembler import Pass_Uce_Guardrails, Score_Contig


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
        weak_long = ["A" * 3000, 1, 0, 10, 3, 2900, 0.5, 0.001, 0.967, 2, 0.1, 1.0]
        dense_short = ["A" * 800, 1, 0, 10, 8, 700, 0.5, 0.01, 0.875, 2, 0.1, 1.0]

        self.assertGreater(Score_Contig(dense_short, "uce"), Score_Contig(weak_long, "uce"))


if __name__ == "__main__":
    unittest.main()
