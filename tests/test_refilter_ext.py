from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from main_refilter_ext import collect_runs_stats, filter_read


class RefilterExtensionFallbackTests(unittest.TestCase):
    def test_filter_read_matches_known_kmer(self):
        self.assertTrue(filter_read("0123", {int("012", 4): 1}, 3))
        self.assertFalse(filter_read("0123", {int("333", 4): 1}, 3))

    def test_collect_runs_stats_counts_hits(self):
        stats = collect_runs_stats("0123", {int("012", 4): 1}, 3)
        self.assertEqual(stats[12], 2)
        self.assertGreaterEqual(stats[9], 1)


if __name__ == "__main__":
    unittest.main()
