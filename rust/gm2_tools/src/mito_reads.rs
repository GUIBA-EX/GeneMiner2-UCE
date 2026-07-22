use crate::fastx::{open_input, FastxFormat, FastxReader, FastxRecord};
use std::io::BufRead;
use std::path::Path;

type FastqRecord = FastxRecord;

/// The finalizer only needs bases after paired identifiers have been checked.
/// Dropping FASTQ header and quality strings here substantially reduces peak
/// memory for high-coverage mitochondrial read pools.
#[derive(Clone, Debug)]
pub struct ReadPair {
    pub first: Vec<u8>,
    pub second: Vec<u8>,
}

fn open_fastq(path: &Path) -> Result<FastxReader, String> {
    FastxReader::open(path, FastxFormat::Fastq).map_err(|error| error.to_string())
}

fn read_record(reader: &mut FastxReader) -> Result<Option<FastqRecord>, String> {
    let Some(mut record) = reader.next_record().map_err(|error| error.to_string())? else {
        return Ok(None);
    };
    // Preserve the previous mitochondrial path's case-normalized comparison.
    record.sequence.make_ascii_uppercase();
    Ok(Some(record))
}

fn fragment_id(header: &[u8]) -> Vec<u8> {
    let header = header.strip_prefix(b"@").unwrap_or(header);
    let id = header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
        .to_vec()
}

fn read_fastq_line(reader: &mut dyn BufRead, line: &mut Vec<u8>) -> Result<bool, String> {
    line.clear();
    if reader
        .read_until(b'\n', line)
        .map_err(|error| error.to_string())?
        == 0
    {
        return Ok(false);
    }
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
    Ok(true)
}

/// Visit interleaved FASTQ pairs without retaining records.  Eight reusable
/// buffers avoid per-read header/quality allocation in the common mito path.
pub fn visit_interleaved_pairs(
    path: &Path,
    mut visit: impl FnMut(&[u8], &[u8]),
) -> Result<usize, String> {
    let mut reader = open_input(path).map_err(|error| error.to_string())?;
    let mut left = std::array::from_fn::<_, 4, _>(|_| Vec::with_capacity(512));
    let mut right = std::array::from_fn::<_, 4, _>(|_| Vec::with_capacity(512));
    let mut count = 0;
    loop {
        // Match FastxReader: ignore empty separators before a FASTQ header.
        loop {
            if !read_fastq_line(&mut reader, &mut left[0])? {
                return Ok(count);
            }
            if !left[0].is_empty() {
                break;
            }
        }
        for line in left.iter_mut().skip(1) {
            if !read_fastq_line(&mut reader, line)? {
                return Err("truncated interleaved FASTQ record".into());
            }
        }
        for line in right.iter_mut() {
            if !read_fastq_line(&mut reader, line)? {
                return Err("interleaved paired FASTQ contains an odd number of records".into());
            }
        }
        if left[0].first() != Some(&b'@')
            || right[0].first() != Some(&b'@')
            || left[2].first() != Some(&b'+')
            || right[2].first() != Some(&b'+')
            || left[1].len() != left[3].len()
            || right[1].len() != right[3].len()
        {
            return Err("invalid interleaved FASTQ record".into());
        }
        if fragment_id(&left[0]) != fragment_id(&right[0]) {
            return Err(format!(
                "interleaved mates do not share an identifier: {} and {}",
                String::from_utf8_lossy(&left[0]),
                String::from_utf8_lossy(&right[0])
            ));
        }
        visit(&left[1], &right[1]);
        count += 1;
    }
}

pub fn read_interleaved_pairs(path: &Path) -> Result<Vec<ReadPair>, String> {
    let mut pairs = Vec::new();
    visit_interleaved_pairs(path, |first, second| {
        pairs.push(ReadPair {
            first: first.to_vec(),
            second: second.to_vec(),
        });
    })?;
    Ok(pairs)
}

fn rc(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Compare two equal-length byte slices, returning `true` when they differ by
/// at most `max_mismatch` positions. One substituted base — a single read error
/// or a benign polymorphism at the seam — no longer discards an otherwise
/// spanning read. Comparison stops as soon as the budget is exceeded.
fn within_mismatch(left: &[u8], right: &[u8], max_mismatch: usize) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let mut mismatches = 0;
    for (a, b) in left.iter().zip(right) {
        if a != b {
            mismatches += 1;
            if mismatches > max_mismatch {
                return false;
            }
        }
    }
    true
}

/// True when `needle` occurs somewhere in `haystack` within `max_mismatch`
/// substitutions. Used to reject a junction window that is not seam-specific.
fn occurs_within_mismatch(haystack: &[u8], needle: &[u8], max_mismatch: usize) -> bool {
    haystack.len() >= needle.len()
        && haystack
            .windows(needle.len())
            .any(|window| within_mismatch(window, needle, max_mismatch))
}

/// Tile the k-mers that straddle the circular seam over a small offset band.
/// A window with `before` bases taken from the end of the sequence and
/// `k - before` from the start is only genuine junction evidence when the same
/// sequence does not also occur inside the linear contig (within one mismatch,
/// so a near-repeat cannot masquerade as a unique closure). Each surviving
/// window is returned with its reverse complement so either read strand counts.
fn junction_windows(circular_sequence: &[u8], k: usize, band: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
    let len = circular_sequence.len();
    let center = k / 2;
    let low = center.saturating_sub(band).max(1);
    let high = (center + band).min(k - 1);
    let mut windows = Vec::new();
    for before in low..=high {
        let mut window = circular_sequence[len - before..].to_vec();
        window.extend_from_slice(&circular_sequence[..k - before]);
        // A window that also appears inside the linear sequence is not evidence
        // that the ends join; reads matching it could be ordinary interior reads.
        if occurs_within_mismatch(circular_sequence, &window, 1) {
            continue;
        }
        let reverse = rc(&window);
        windows.push((window, reverse));
    }
    windows
}

/// Minimum number of reads that span the circular seam, measured across a band
/// of straddling k-mers rather than a single central one. A closure is only as
/// trustworthy as its least-supported seam position, so the per-window counts
/// are reduced with `min`: one lucky (or erroneous) k-mer can no longer stand
/// in for consistent spanning coverage. Matching tolerates one mismatch on
/// either strand. Zero is returned when no window is seam-specific, keeping an
/// ambiguous repeat-bounded contig linear.
pub fn count_junction_support(
    reads: &Path,
    circular_sequence: &[u8],
    k: usize,
) -> Result<usize, String> {
    if k < 3 || circular_sequence.len() < k {
        return Ok(0);
    }
    // A narrow band keeps the off-centre windows reachable by ordinary reads
    // while still demanding the seam be crossed consistently, not at one point.
    let band = ((k - 1) / 2).saturating_sub(1).min(3);
    let windows = junction_windows(circular_sequence, k, band);
    if windows.is_empty() {
        return Ok(0);
    }
    let mut per_window = vec![0_usize; windows.len()];
    let mut reader = open_fastq(reads)?;
    while let Some(record) = read_record(&mut reader)? {
        for (index, (forward, reverse)) in windows.iter().enumerate() {
            let spans = record.sequence.windows(k).any(|read_window| {
                within_mismatch(read_window, forward, 1) || within_mismatch(read_window, reverse, 1)
            });
            if spans {
                per_window[index] += 1;
            }
        }
    }
    Ok(per_window.into_iter().min().unwrap_or(0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragment_names_normalize_mates() {
        assert_eq!(fragment_id(b"@read42/1 comment"), b"read42".to_vec());
        assert_eq!(fragment_id(b"@read42/2"), b"read42".to_vec());
    }

    #[test]
    fn streaming_pair_visitor_preserves_validated_mates() {
        let path = std::env::temp_dir().join(format!(
            "gm2_pair_visit_{}_{}.fq",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, "\n@x/1\nAAAA\n+\nFFFF\n@x/2\nTTTT\n+\nFFFF\n").unwrap();
        let mut seen = Vec::new();
        let count = visit_interleaved_pairs(&path, |first, second| {
            seen.push((first.to_vec(), second.to_vec()));
        })
        .unwrap();
        assert_eq!(count, 1);
        assert_eq!(seen, vec![(b"AAAA".to_vec(), b"TTTT".to_vec())]);
        std::fs::remove_file(path).unwrap();
    }

    fn temp_fastq(tag: &str, body: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "gm2_junction_{tag}_{}_{}.fq",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn junction_support_requires_a_spanning_read() {
        // Sequence of length 24; k = 8 tiles a small band around the seam, so a
        // supporting read must cross the join by enough on each side.
        let sequence = b"AAAACCCCGGGGTTTTACGTACGT";
        // Read spanning the seam: last 8 bases + first 8 bases of the circle.
        let path = temp_fastq(
            "span",
            "@x/1\nACGTACGTAAAACCCC\n+\nFFFFFFFFFFFFFFFF\n@x/2\nCCCCGGGG\n+\nFFFFFFFF\n",
        );
        assert_eq!(count_junction_support(&path, sequence, 8).unwrap(), 1);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn junction_support_tolerates_a_single_mismatch() {
        let sequence = b"AAAACCCCGGGGTTTTACGTACGT";
        // Same spanning read as above with one substituted base at the seam.
        let path = temp_fastq("mismatch", "@x/1\nACGTACGTAATACCCC\n+\nFFFFFFFFFFFFFFFF\n");
        assert_eq!(count_junction_support(&path, sequence, 8).unwrap(), 1);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn junction_support_is_the_minimum_spanning_depth() {
        let sequence = b"AAAACCCCGGGGTTTTACGTACGT";
        // Two reads cross the seam; a third, purely interior read touches
        // neither side of the join and must not lift the weakest seam position,
        // so the reported depth stays at 2.
        let path = temp_fastq(
            "min",
            "@a/1\nACGTACGTAAAACCCC\n+\nFFFFFFFFFFFFFFFF\n\
             @b/1\nACGTACGTAAAACCCC\n+\nFFFFFFFFFFFFFFFF\n\
             @c/1\nAAAACCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n",
        );
        assert_eq!(count_junction_support(&path, sequence, 8).unwrap(), 2);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn junction_support_ignores_a_seam_that_recurs_internally() {
        // A periodic contig: every seam k-mer also appears in the interior, so
        // spanning reads are not unique evidence of closure. Report zero and
        // keep the contig linear even though reads cross the (repeated) join.
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let path = temp_fastq("internal", "@x/1\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n");
        assert_eq!(count_junction_support(&path, sequence, 8).unwrap(), 0);
        std::fs::remove_file(path).unwrap();
    }
}
