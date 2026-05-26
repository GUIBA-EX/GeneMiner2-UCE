from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import get_uce_rescue_parallelism


class UceRescueParallelismTests(unittest.TestCase):
    def test_uses_four_by_four_when_threads_allow(self):
        self.assertEqual(get_uce_rescue_parallelism(32, 45), (4, 4))
        self.assertEqual(get_uce_rescue_parallelism(16, 45), (4, 4))

    def test_scales_down_for_small_thread_counts(self):
        self.assertEqual(get_uce_rescue_parallelism(8, 45), (2, 4))
        self.assertEqual(get_uce_rescue_parallelism(4, 45), (1, 4))
        self.assertEqual(get_uce_rescue_parallelism(2, 45), (1, 2))

    def test_scales_down_for_small_sample_counts(self):
        self.assertEqual(get_uce_rescue_parallelism(32, 2), (2, 4))


if __name__ == "__main__":
    unittest.main()
