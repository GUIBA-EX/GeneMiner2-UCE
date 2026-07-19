use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: Vec<u8>,
    pub plus: String,
    pub quality: String,
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path).map_err(|error| error.to_string())?;
    let reader: Box<dyn Read> = if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
    {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

fn read_record(reader: &mut dyn BufRead) -> Result<Option<FastqRecord>, String> {
    let mut lines = [String::new(), String::new(), String::new(), String::new()];
    if reader
        .read_line(&mut lines[0])
        .map_err(|error| error.to_string())?
        == 0
    {
        return Ok(None);
    }
    for line in lines.iter_mut().skip(1) {
        if reader.read_line(line).map_err(|error| error.to_string())? == 0 {
            return Err("truncated FASTQ record".into());
        }
    }
    for line in &mut lines {
        while line.ends_with(['\n', '\r']) {
            line.pop();
        }
    }
    if !lines[0].starts_with('@') || !lines[2].starts_with('+') {
        return Err("invalid FASTQ record".into());
    }
    Ok(Some(FastqRecord {
        header: std::mem::take(&mut lines[0]),
        sequence: lines[1]
            .bytes()
            .map(|base| base.to_ascii_uppercase())
            .collect(),
        plus: std::mem::take(&mut lines[2]),
        quality: std::mem::take(&mut lines[3]),
    }))
}

fn fragment_id(header: &str) -> String {
    header
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}

pub fn read_interleaved_pairs(path: &Path) -> Result<Vec<(FastqRecord, FastqRecord)>, String> {
    let mut reader = open_fastq(path)?;
    let mut pairs = Vec::new();
    loop {
        let Some(first) = read_record(reader.as_mut())? else {
            break;
        };
        let second = read_record(reader.as_mut())?
            .ok_or("interleaved paired FASTQ contains an odd number of records")?;
        if fragment_id(&first.header) != fragment_id(&second.header) {
            return Err(format!(
                "interleaved mates do not share an identifier: {} and {}",
                first.header, second.header
            ));
        }
        pairs.push((first, second));
    }
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
    while let Some(record) = read_record(reader.as_mut())? {
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
        assert_eq!(fragment_id("@read42/1 comment"), "read42");
        assert_eq!(fragment_id("@read42/2"), "read42");
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
