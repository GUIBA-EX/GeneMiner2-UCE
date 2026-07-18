use std::io::{self, BufRead, Write};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: Vec<u8>,
}

/// 读取普通 FASTA。序列中的空白会被移除，与 Biopython SimpleFastaParser 一致。
pub fn read_fasta(reader: impl BufRead) -> io::Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut name: Option<String> = None;
    let mut sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = name.replace(header.trim().to_string()) {
                records.push(FastaRecord {
                    name: previous,
                    sequence: std::mem::take(&mut sequence),
                });
            }
        } else if name.is_some() {
            sequence.extend(line.bytes().filter(|base| !base.is_ascii_whitespace()));
        }
    }

    if let Some(name) = name {
        records.push(FastaRecord { name, sequence });
    }
    Ok(records)
}

pub fn write_fasta(mut writer: impl Write, records: &[FastaRecord]) -> io::Result<()> {
    for record in records {
        writer.write_all(b">")?;
        writer.write_all(record.name.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&record.sequence)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parses_multiline_fasta_and_preserves_full_title() {
        let input = Cursor::new(b"ignored\n>sample one  \nAC G\nT\tA\n>sample2\nNN\n");
        assert_eq!(
            read_fasta(input).unwrap(),
            vec![
                FastaRecord {
                    name: "sample one".to_string(),
                    sequence: b"ACGTA".to_vec(),
                },
                FastaRecord {
                    name: "sample2".to_string(),
                    sequence: b"NN".to_vec(),
                },
            ]
        );
    }
}
