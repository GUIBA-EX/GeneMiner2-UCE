from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import density_ratio_or_blank, read_density_or_blank, rescue_density_below_ratio


class UceRescueFallbackTests(unittest.TestCase):
    def test_large_density_drop_triggers_fallback(self):
        before = {'selected_contig_length': '326', 'read_count': '14'}
        after = {'selected_contig_length': '14426', 'read_count': '13'}
        self.assertTrue(rescue_density_below_ratio(before, after, 0.5))

    def test_density_increase_keeps_rescue(self):
        before = {'selected_contig_length': '532', 'read_count': '10'}
        after = {'selected_contig_length': '1476', 'read_count': '38'}
        self.assertFalse(rescue_density_below_ratio(before, after, 0.5))

    def test_moderate_density_drop_keeps_rescue(self):
        before = {'selected_contig_length': '100', 'read_count': '10'}
        after = {'selected_contig_length': '160', 'read_count': '9'}
        self.assertAlmostEqual(density_ratio_or_blank(before, after), 0.5625)
        self.assertFalse(rescue_density_below_ratio(before, after, 0.5))

    def test_blank_density_does_not_trigger_fallback(self):
        self.assertEqual(read_density_or_blank({'selected_contig_length': '0', 'read_count': '5'}), '')
        self.assertFalse(rescue_density_below_ratio({}, {'selected_contig_length': '10', 'read_count': '1'}, 0.5))

    def test_unique_read_density_is_preferred(self):
        before = {
            'selected_contig_length': '100',
            'read_count': '100',
            'unique_read_count': '10',
        }
        after = {
            'selected_contig_length': '100',
            'read_count': '1000',
            'unique_read_count': '1',
        }

        self.assertAlmostEqual(density_ratio_or_blank(before, after), 0.1)
        self.assertTrue(rescue_density_below_ratio(before, after, 0.5))
        self.assertEqual(read_density_or_blank({'unique_read_density': '0.125'}), 0.125)


if __name__ == "__main__":
    unittest.main()
