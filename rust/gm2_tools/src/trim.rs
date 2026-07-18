use crate::fasta::{read_fasta, write_fasta, FastaRecord};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter};
use std::path::Path;
use std::process::Command;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrimMode {
    All,
    Longest,
    Terminal,
    Isoform,
}

#[derive(Clone, Debug, PartialEq)]
struct SequenceMatch {
    start: usize,
    end: usize,
    reverse: bool,
}

impl SequenceMatch {
    fn from_fields(fields: &[&str]) -> Option<Self> {
        if fields.len() < 12 {
            return None;
        }
        let query_start: usize = fields[6].parse().ok()?;
        let query_end: usize = fields[7].parse().ok()?;
        let subject_start: usize = fields[8].parse().ok()?;
        let subject_end: usize = fields[9].parse().ok()?;
        Some(Self {
            start: query_start.min(query_end),
            end: query_start.max(query_end),
            reverse: (query_start > query_end) != (subject_start > subject_end),
        })
    }

    fn length(&self) -> usize {
        self.end - self.start + 1
    }
}

fn read_matches(text: &str) -> Vec<SequenceMatch> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| SequenceMatch::from_fields(&line.split('\t').collect::<Vec<_>>()))
        .collect()
}

fn merge_matches(mut matches: Vec<SequenceMatch>) -> Vec<SequenceMatch> {
    matches.sort_by_key(|item| item.start);
    let mut merged = Vec::new();
    for reverse in [false, true] {
        let mut oriented = matches
            .iter()
            .filter(|item| item.reverse == reverse)
            .cloned();
        let Some(mut current) = oriented.next() else {
            continue;
        };
        for item in oriented {
            if item.start <= current.end {
                current.end = current.end.max(item.end);
            } else {
                merged.push(current);
                current = item;
            }
        }
        merged.push(current);
    }
    merged
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' | b'U' => b'A',
        b't' | b'u' => b'a',
        b'M' => b'K',
        b'm' => b'k',
        b'R' => b'Y',
        b'r' => b'y',
        b'W' => b'W',
        b'w' => b'w',
        b'S' => b'S',
        b's' => b's',
        b'Y' => b'R',
        b'y' => b'r',
        b'K' => b'M',
        b'k' => b'm',
        b'V' => b'B',
        b'v' => b'b',
        b'H' => b'D',
        b'h' => b'd',
        b'D' => b'H',
        b'd' => b'h',
        b'B' => b'V',
        b'b' => b'v',
        other => other,
    }
}

fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = complement(*base);
    }
}

fn median_reference_length(records: &[FastaRecord]) -> Option<f64> {
    let mut lengths: Vec<usize> = records.iter().map(|record| record.sequence.len()).collect();
    lengths.sort_unstable();
    match lengths.len() {
        0 => None,
        count if count % 2 == 1 => Some(lengths[count / 2] as f64),
        count => Some((lengths[count / 2 - 1] + lengths[count / 2]) as f64 / 2.0),
    }
}

pub fn trim_record(
    query: &FastaRecord,
    reference_records: &[FastaRecord],
    blast_output: &str,
    retention_percentage: f64,
    mode: TrimMode,
) -> Option<FastaRecord> {
    let median_length = median_reference_length(reference_records)?;
    let mut matches = read_matches(blast_output);
    if mode == TrimMode::All {
        matches = merge_matches(matches);
    }
    if matches.is_empty() {
        return None;
    }

    // 正反向 HSP 不能拼进同一条序列；总覆盖长度相同的时候保持正向。
    let forward: usize = matches
        .iter()
        .filter(|item| !item.reverse)
        .map(SequenceMatch::length)
        .sum();
    let reverse: usize = matches
        .iter()
        .filter(|item| item.reverse)
        .map(SequenceMatch::length)
        .sum();
    let selected_reverse = reverse > forward;
    matches.retain(|item| item.reverse == selected_reverse);

    match mode {
        TrimMode::Longest | TrimMode::Isoform => {
            let longest = matches.into_iter().max_by_key(SequenceMatch::length)?;
            matches = vec![longest];
        }
        TrimMode::Terminal => {
            matches = vec![SequenceMatch {
                start: matches.iter().map(|item| item.start).min()?,
                end: matches.iter().map(|item| item.end).max()?,
                reverse: selected_reverse,
            }];
        }
        TrimMode::All => {}
    }
    matches.sort_by_key(|item| item.start);

    let mut sequence = Vec::new();
    for item in matches {
        if item.start == 0 || item.end > query.sequence.len() {
            return None;
        }
        sequence.extend_from_slice(&query.sequence[item.start - 1..item.end]);
    }
    if selected_reverse {
        reverse_complement(&mut sequence);
    }
    if sequence.len() as f64 / median_length * 100.0 <= retention_percentage {
        return None;
    }
    Some(FastaRecord {
        name: query.name.clone(),
        sequence,
    })
}

pub fn run_trim(
    query_path: &Path,
    reference_path: &Path,
    output_path: &Path,
    database: &Path,
    executable: &Path,
    retention_percentage: f64,
    mode: TrimMode,
) -> io::Result<()> {
    match fs::remove_file(output_path) {
        Ok(()) => {}
        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
        Err(error) => return Err(error),
    }
    let query_records = read_fasta(BufReader::new(File::open(query_path)?))?;
    let Some(query) = query_records.first() else {
        return Ok(());
    };
    let references = read_fasta(BufReader::new(File::open(reference_path)?))?;
    if references.is_empty() {
        return Ok(());
    }

    let mut command = Command::new(executable);
    command
        .env("BLAST_USAGE_REPORT", "0")
        .env("DO_NOT_TRACK", "1");
    if mode == TrimMode::Isoform {
        command.args([
            "-query",
            query_path.to_string_lossy().as_ref(),
            "-db",
            database.to_string_lossy().as_ref(),
            "-outfmt",
            "tabular",
            "-word_size",
            "13",
            "-score",
            "20",
            "-limit_lookup",
            "F",
            "-penalty",
            "-2",
        ]);
    } else {
        command.args([
            "-query",
            query_path.to_string_lossy().as_ref(),
            "-db",
            database.to_string_lossy().as_ref(),
            "-outfmt",
            "6",
            "-word_size",
            "20",
            "-min_raw_gapped_score",
            "20",
        ]);
    }
    let output = command.output()?;
    if !output.status.success() {
        return Err(io::Error::other(format!(
            "alignment program exited with {}",
            output.status
        )));
    }
    let blast_output = String::from_utf8_lossy(&output.stdout);
    if let Some(record) = trim_record(
        query,
        &references,
        &blast_output,
        retention_percentage,
        mode,
    ) {
        write_fasta(BufWriter::new(File::create(output_path)?), &[record])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(sequence: &[u8]) -> FastaRecord {
        FastaRecord {
            name: "query".to_string(),
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn normalizes_reverse_coordinates_and_reverse_complements() {
        let blast = "q\ts\t100\t4\t0\t0\t8\t5\t1\t4\t0\t100\n";
        let result = trim_record(
            &record(b"AACCGGTT"),
            &[record(b"AACCGGTT")],
            blast,
            0.0,
            TrimMode::Longest,
        )
        .unwrap();
        assert_eq!(result.sequence, b"AACC");
    }

    #[test]
    fn merges_overlapping_hsps_in_all_mode() {
        let blast =
            "q\ts\t100\t4\t0\t0\t1\t4\t1\t4\t0\t100\nq\ts\t100\t4\t0\t0\t3\t6\t3\t6\t0\t100\n";
        let result = trim_record(
            &record(b"AACCGGTT"),
            &[record(b"AACCGGTT")],
            blast,
            0.0,
            TrimMode::All,
        )
        .unwrap();
        assert_eq!(result.sequence, b"AACCGG");
    }
}
