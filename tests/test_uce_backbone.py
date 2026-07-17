from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from main_assembler import (
    ACGT_DICT,
    Get_Forward_Uce_Backbone,
    Get_Forward_Contig_v6,
    Int_To_Seq,
    Seq_To_Int,
)


def encode(sequence):
    return Seq_To_Int(sequence)[0][0]


def graph(entries):
    return {
        encode(sequence): [depth, 100, False, 0]
        for sequence, depth in entries
    }


def extension_sequence(contigs):
    return ''.join(ACGT_DICT[base] for base in contigs[0][1])


class UceBackboneTests(unittest.TestCase):
    def test_commits_longer_bubble_arm_without_backtracking(self):
        kmer_graph = graph([
            ("AAAA", 5),
            ("AAAC", 2),
            ("AACG", 2),
            ("ACGT", 2),
            ("CGTT", 2),
            ("GTTT", 2),
            ("AAAG", 20),
            ("AAGA", 20),
            ("AGAC", 20),
        ])

        contigs, visited, _, _ = Get_Forward_Uce_Backbone(
            kmer_graph, encode("AAAA"), 4, iteration=100, lookahead=24
        )
        legacy_contigs, _, _, _ = Get_Forward_Contig_v6(
            kmer_graph, encode("AAAA"), 4, iteration=100
        )

        self.assertEqual(len(contigs), 1)
        self.assertGreater(len(legacy_contigs), len(contigs))
        self.assertEqual(extension_sequence(contigs), "CGTTT")
        self.assertIn(encode("AAAC"), visited)
        self.assertNotIn(encode("AAAG"), visited)

    def test_equal_length_bubble_prefers_cumulative_support(self):
        kmer_graph = graph([
            ("CCCC", 5),
            ("CCCA", 2),
            ("CCAA", 2),
            ("CCCG", 20),
            ("CCGT", 20),
        ])

        contigs, visited, _, _ = Get_Forward_Uce_Backbone(
            kmer_graph, encode("CCCC"), 4, iteration=100, lookahead=24
        )

        self.assertEqual(extension_sequence(contigs), "GT")
        self.assertIn(encode("CCCG"), visited)
        self.assertNotIn(encode("CCCA"), visited)

    def test_cycle_stops_without_revisiting_kmers(self):
        kmer_graph = graph([
            ("AAA", 5),
            ("AAC", 5),
            ("ACA", 5),
            ("CAA", 5),
        ])

        contigs, visited, _, _ = Get_Forward_Uce_Backbone(
            kmer_graph, encode("AAA"), 3, iteration=100, lookahead=24
        )

        self.assertEqual(extension_sequence(contigs), "CAA")
        self.assertEqual(
            {Int_To_Seq(kmer, 3) for kmer in visited},
            {"AAA", "AAC", "ACA", "CAA"},
        )


if __name__ == "__main__":
    unittest.main()
