from pathlib import Path
from types import SimpleNamespace
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import build_uce_rescue_filter_commands, get_uce_rescue_parallelism


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

    def test_rescue_filter_commands_use_assigned_threads(self):
        args = SimpleNamespace(kf=21, step_size=4, max_reads=0)
        dict_cmd, reads_cmd = build_uce_rescue_filter_commands(
            "MainFilterNew",
            "/tmp/rescue_refs",
            "/tmp/sample",
            "/tmp/R1.fq.gz",
            "/tmp/R2.fq.gz",
            args,
            "/tmp/sample/uce_rescue_kmer_dict_k21.dict",
            4,
        )

        self.assertEqual(dict_cmd[dict_cmd.index("-m") + 1], "4")
        self.assertEqual(reads_cmd[reads_cmd.index("-m") + 1], "4")
        self.assertNotIn("5", reads_cmd)

    def test_rescue_filter_commands_preserve_max_reads(self):
        args = SimpleNamespace(kf=21, step_size=4, max_reads=1000)
        _, reads_cmd = build_uce_rescue_filter_commands(
            "MainFilterNew",
            "/tmp/rescue_refs",
            "/tmp/sample",
            "/tmp/R1.fq.gz",
            "/tmp/R2.fq.gz",
            args,
            "/tmp/sample/uce_rescue_kmer_dict_k21.dict",
            2,
        )

        self.assertEqual(reads_cmd[reads_cmd.index("-m_reads") + 1], "1000")


if __name__ == "__main__":
    unittest.main()
