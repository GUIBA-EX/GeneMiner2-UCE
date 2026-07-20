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

pub fn count_junction_support(
    reads: &Path,
    circular_sequence: &[u8],
    k: usize,
) -> Result<usize, String> {
    if k < 3 || circular_sequence.len() < k {
        return Ok(0);
    }
    let left = k / 2;
    let mut junction = circular_sequence[circular_sequence.len() - left..].to_vec();
    junction.extend_from_slice(&circular_sequence[..k - left]);
    let reverse = rc(&junction);
    if circular_sequence
        .windows(k)
        .any(|window| window == junction || window == reverse)
    {
        return Ok(0);
    }
    let mut reader = open_fastq(reads)?;
    let mut support = 0;
    while let Some(record) = read_record(&mut reader)? {
        if record
            .sequence
            .windows(k)
            .any(|window| window == junction || window == reverse)
        {
            support += 1;
        }
    }
    Ok(support)
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

    #[test]
    fn junction_support_requires_a_spanning_read() {
        let path = std::env::temp_dir().join(format!(
            "gm2_junction_{}_{}.fq",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(
            &path,
            "@x/1\nGGTTTTAAAACC\n+\nFFFFFFFFFFFF\n@x/2\nCCCCGGGG\n+\nFFFFFFFF\n",
        )
        .unwrap();
        assert_eq!(
            count_junction_support(&path, b"AAAACCCCGGGGTTTT", 8).unwrap(),
            1
        );
        std::fs::remove_file(path).unwrap();
    }
}
