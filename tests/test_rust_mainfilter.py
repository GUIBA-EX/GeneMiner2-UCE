import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "rust" / "main_filter_new" / "Cargo.toml"
BINARY = ROOT / "rust" / "main_filter_new" / "target" / "debug" / "MainFilterNew"


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def decode_gm2(path):
    records = []
    data = path.read_bytes()
    position = 0
    while position < len(data):
        header = data[position : position + 6]
        if len(header) != 6:
            raise AssertionError("truncated GM2 header")
        position += 6
        record_len = (header[0] << 16) | (header[1] << 8) | header[2]
        has_quality = bool(header[3] & 0x80)
        sequence_len = ((header[3] & 0x7F) << 16) | (header[4] << 8) | header[5]
        payload = data[position : position + record_len]
        position += record_len

        sequence = []
        payload_position = 0
        last_chunk = 0
        while len(sequence) < sequence_len:
            encoded = payload[payload_position]
            payload_position += 1
            for _ in range(2 if encoded & 0x80 else 1):
                chunk = (encoded & 0x7F) ^ last_chunk
                sequence.extend(chr((chunk & 31) + 64) for _ in range((chunk >> 5) + 1))
                last_chunk = chunk

        quality = []
        if has_quality:
            while len(quality) < sequence_len:
                encoded = payload[payload_position]
                payload_position += 1
                if encoded & 0x80:
                    value = payload[payload_position]
                    payload_position += 1
                    repeat = ((encoded & 0x7F) | (value & 0x80)) + 1
                    quality.extend(chr(value & 0x7F) for _ in range(repeat))
                else:
                    quality.append(chr(encoded))
        records.append(("".join(sequence[:sequence_len]), "".join(quality[:sequence_len])))
    return records


class RustMainFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("cargo") is None:
            raise unittest.SkipTest("cargo is not installed")
        subprocess.run(
            ["cargo", "build", "--quiet", "--manifest-path", str(MANIFEST)],
            check=True,
            cwd=ROOT,
        )

    def test_paired_simulation_cache_text_and_gm2(self):
        locus_a = "ACGTTGCATGTCAGTACGATCGTACCTGACGTAGCTAGCA"
        locus_b = "TTGGAACCTTCCGGAATTCGATATCGCGATATGCCGTAGC"
        read_a = locus_a[4:28]
        read_b = locus_b[3:27]
        read_a_reverse = reverse_complement(locus_a[8:32])
        read_1 = [read_a, read_b, read_a_reverse, "N" * 24]
        read_2 = ["N" * 24, "C" * 24, "N" * 24, "G" * 24]

        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            references = work / "references"
            references.mkdir()
            (references / "locus_a.fasta").write_text(f">a\n{locus_a}\n")
            (references / "locus_b.fasta").write_text(f">b\n{locus_b}\n")
            q1 = work / "reads_1.fq.gz"
            q2 = work / "reads_2.fq"
            with gzip.open(q1, "wt") as output:
                for index, sequence in enumerate(read_1, 1):
                    output.write(f"@pair{index}/1\n{sequence}\n+\n{'I' * len(sequence)}\n")
            with q2.open("w") as output:
                for index, sequence in enumerate(read_2, 1):
                    output.write(f"@pair{index}/2\n{sequence}\n+\n{'J' * len(sequence)}\n")

            dictionary = work / "reference.dict"
            index_output = work / "index"
            subprocess.run(
                [
                    str(BINARY), "-r", str(references), "-o", str(index_output),
                    "-kf", "16", "-s", "1", "-gr", "-lkd", str(dictionary), "-m", "2",
                ],
                check=True,
                cwd=ROOT,
            )
            self.assertEqual(dictionary.read_bytes()[:4], b"GM2K")

            text_output = work / "text"
            subprocess.run(
                [
                    str(BINARY), "-r", str(references), "-q1", str(q1), "-q2", str(q2),
                    "-o", str(text_output), "-kf", "16", "-s", "1", "-lkd", str(dictionary),
                    "-subdir", "filtered_pe", "-m", "4",
                ],
                check=True,
                cwd=ROOT,
            )
            self.assertEqual(
                (text_output / "ref_reads_count_dict.txt").read_text().splitlines(),
                ["locus_a,4", "locus_b,2"],
            )
            locus_a_r1 = (text_output / "filtered_pe" / "locus_a_1.fq").read_text()
            locus_a_r2 = (text_output / "filtered_pe" / "locus_a_2.fq").read_text()
            self.assertIn("@pair1/1", locus_a_r1)
            self.assertIn("@pair3/1", locus_a_r1)
            self.assertNotIn("@pair2/1", locus_a_r1)
            self.assertIn("@pair1/2", locus_a_r2)
            self.assertIn("@pair3/2", locus_a_r2)

            gm2_output = work / "gm2"
            subprocess.run(
                [
                    str(BINARY), "-r", str(references), "-q1", str(q1), "-q2", str(q2),
                    "-o", str(gm2_output), "-kf", "16", "-s", "1", "-gr",
                    "-lkd", str(dictionary), "-subdir", "filtered_pe", "-m", "5", "-lb",
                ],
                check=True,
                cwd=ROOT,
            )
            self.assertEqual(
                decode_gm2(gm2_output / "filtered_pe" / "locus_a_1.gm2"),
                [(read_a, "I" * 24), (read_a_reverse, "I" * 24)],
            )
            self.assertEqual(
                decode_gm2(gm2_output / "filtered_pe" / "locus_a_2.gm2"),
                [("N" * 24, "J" * 24), ("N" * 24, "J" * 24)],
            )

    def test_remaining_output_modes(self):
        reference = "ACGTTGCATGTCAGTACGATCGTACCTGACGTAGCTAGCA"
        matching = reference[4:28]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref_dir = work / "references"
            ref_dir.mkdir()
            (ref_dir / "locus.fasta").write_text(f">locus\n{reference}\n")
            q1 = work / "reads_1.fq"
            q2 = work / "reads_2.fq"
            q1.write_text(f"@hit/1\n{matching}\n+\n{'I' * 24}\n")
            q2.write_text(f"@hit/2\n{'N' * 24}\n+\n{'J' * 24}\n")

            for mode in (0, 1, 3):
                output = work / f"mode_{mode}"
                subprocess.run(
                    [
                        str(BINARY), "-r", str(ref_dir), "-q1", str(q1), "-q2", str(q2),
                        "-o", str(output), "-kf", "16", "-s", "1", "-gr", "-m", str(mode),
                    ],
                    check=True,
                    cwd=ROOT,
                )
                self.assertEqual(
                    (output / "ref_reads_count_dict.txt").read_text(),
                    "locus,2\n",
                )

            interleaved = (work / "mode_0" / "filtered" / "locus.fq").read_text()
            self.assertIn("@hit/1", interleaved)
            self.assertIn("@hit/2", interleaved)
            self.assertTrue((work / "mode_1" / "filtered" / "all_1.fq").is_file())
            self.assertTrue((work / "mode_1" / "filtered" / "all_2.fq").is_file())
            self.assertEqual(list((work / "mode_3" / "filtered").iterdir()), [])

    def test_long_kmer_cache_roundtrip(self):
        reference = "ACGTTGCATGTCAGTACGATCGTACCTGACGTAGCTAGCATGGCATACGTTAGCCATGCACTGA"
        read = reference[8:56]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref_dir = work / "references"
            ref_dir.mkdir()
            (ref_dir / "long_locus.fasta").write_text(f">long\n{reference}\n")
            reads = work / "reads.fq"
            reads.write_text(f"@long\n{read}\n+\n{'I' * len(read)}\n")
            dictionary = work / "long.dict"

            subprocess.run(
                [
                    str(BINARY), "-r", str(ref_dir), "-o", str(work / "index"),
                    "-kf", "33", "-gr", "-lkd", str(dictionary), "-m", "2",
                ],
                check=True,
                cwd=ROOT,
            )
            output = work / "filtered"
            subprocess.run(
                [
                    str(BINARY), "-r", str(ref_dir), "-q1", str(reads), "-o", str(output),
                    "-kf", "33", "-lkd", str(dictionary), "-m", "3",
                ],
                check=True,
                cwd=ROOT,
            )
            self.assertEqual(
                (output / "ref_reads_count_dict.txt").read_text(),
                "long_locus,1\n",
            )

    def test_one_encoded_record_is_reused_for_multiple_loci(self):
        reference = "ACGTTGCATGTCAGTACGATCGTACCTGACGTAGCTAGCA"
        read = reference[4:36]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref_dir = work / "references"
            ref_dir.mkdir()
            (ref_dir / "locus_a.fasta").write_text(f">a\n{reference}\n")
            (ref_dir / "locus_b.fasta").write_text(f">b\n{reference}\n")
            reads = work / "reads.fq"
            reads.write_text(f"@shared\n{read}\n+\n{'I' * len(read)}\n")
            output = work / "output"

            subprocess.run(
                [
                    str(BINARY), "-r", str(ref_dir), "-q1", str(reads), "-o", str(output),
                    "-kf", "16", "-s", "3", "-gr", "-subdir", "filtered", "-m", "5",
                ],
                check=True,
                cwd=ROOT,
            )

            expected = [(read, "I" * len(read))]
            self.assertEqual(decode_gm2(output / "filtered" / "locus_a_1.gm2"), expected)
            self.assertEqual(decode_gm2(output / "filtered" / "locus_b_1.gm2"), expected)
            self.assertEqual(
                (output / "ref_reads_count_dict.txt").read_text().splitlines(),
                ["locus_a,1", "locus_b,1"],
            )


if __name__ == "__main__":
    unittest.main()
