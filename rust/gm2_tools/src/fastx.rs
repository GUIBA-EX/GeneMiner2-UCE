//! Buffered byte-level FASTA/FASTQ readers shared by GeneMiner2-UCE tools.
//! Text is kept as bytes: sequence tools do not need UTF-8 validation or a
//! temporary `String` for every FASTQ line.
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

pub const READ_BUFFER_SIZE: usize = 1024 * 1024;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FastxRecord {
    /// Includes the leading `>` or `@` exactly as it appeared in the input.
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
    /// Empty for FASTA records; includes the leading `+` for FASTQ records.
    pub plus: Vec<u8>,
    /// Empty for FASTA records.
    pub quality: Vec<u8>,
}

pub struct FastxReader {
    input: BufReader<Box<dyn Read>>,
    format: FastxFormat,
    pending_header: Option<Vec<u8>>,
    scratch: Vec<u8>,
    finished: bool,
}

fn is_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
}

pub fn open_input(path: &Path) -> io::Result<BufReader<Box<dyn Read>>> {
    let file = File::open(path)?;
    let input: Box<dyn Read> = if is_gzip(path) {
        // Buffer both compressed input and decompressed output. This preserves
        // flate2's concatenated-gzip support while avoiding 8 KiB read cycles.
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            READ_BUFFER_SIZE,
            file,
        )))
    } else {
        Box::new(file)
    };
    Ok(BufReader::with_capacity(READ_BUFFER_SIZE, input))
}

impl FastxReader {
    pub fn open(path: &Path, format: FastxFormat) -> io::Result<Self> {
        Ok(Self {
            input: open_input(path)?,
            format,
            pending_header: None,
            scratch: Vec::with_capacity(512),
            finished: false,
        })
    }

    fn read_line(&mut self) -> io::Result<Option<Vec<u8>>> {
        self.scratch.clear();
        if self.input.read_until(b'\n', &mut self.scratch)? == 0 {
            return Ok(None);
        }
        while matches!(self.scratch.last(), Some(b'\n' | b'\r')) {
            self.scratch.pop();
        }
        Ok(Some(std::mem::take(&mut self.scratch)))
    }

    pub fn next_record(&mut self) -> io::Result<Option<FastxRecord>> {
        match self.format {
            FastxFormat::Fasta => self.next_fasta(),
            FastxFormat::Fastq => self.next_fastq(),
        }
    }

    fn next_fastq(&mut self) -> io::Result<Option<FastxRecord>> {
        let header = loop {
            let Some(line) = self.read_line()? else {
                return Ok(None);
            };
            if !line.is_empty() {
                break line;
            }
        };
        if header.first() != Some(&b'@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed FASTQ record",
            ));
        }
        let sequence = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ sequence")
        })?;
        let plus = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ plus line")
        })?;
        if plus.first() != Some(&b'+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed FASTQ plus line",
            ));
        }
        let quality = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ quality")
        })?;
        if quality.len() != sequence.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ sequence and quality lengths differ",
            ));
        }
        Ok(Some(FastxRecord {
            header,
            sequence,
            plus,
            quality,
        }))
    }

    fn next_fasta(&mut self) -> io::Result<Option<FastxRecord>> {
        if self.finished {
            return Ok(None);
        }
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                let Some(line) = self.read_line()? else {
                    self.finished = true;
                    return Ok(None);
                };
                if line.first() == Some(&b'>') {
                    break line;
                }
                if !line.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "FASTA sequence encountered before a header",
                    ));
                }
            }
        };
        let mut sequence = Vec::new();
        loop {
            let Some(line) = self.read_line()? else {
                self.finished = true;
                break;
            };
            if line.first() == Some(&b'>') {
                self.pending_header = Some(line);
                break;
            }
            sequence.extend_from_slice(&line);
        }
        Ok(Some(FastxRecord {
            header,
            sequence,
            plus: Vec::new(),
            quality: Vec::new(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_fastq_bytes_and_crlf() {
        let path = std::env::temp_dir().join(format!("gm2-fastx-{}", std::process::id()));
        let mut file = File::create(&path).unwrap();
        file.write_all(b"@read/1\r\nAcGT\r\n+\r\n!!!!\r\n").unwrap();
        drop(file);
        let mut reader = FastxReader::open(&path, FastxFormat::Fastq).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.header, b"@read/1");
        assert_eq!(record.sequence, b"AcGT");
        assert_eq!(record.quality, b"!!!!");
        assert!(reader.next_record().unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn parses_gzip_fastq() {
        use flate2::{write::GzEncoder, Compression};
        let path =
            std::env::temp_dir().join(format!("gm2-fastx-gzip-{}.fq.gz", std::process::id()));
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(b"@read\nACGT\n+\n!!!!\n").unwrap();
        encoder.finish().unwrap();
        let mut reader = FastxReader::open(&path, FastxFormat::Fastq).unwrap();
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"ACGT");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn parses_multiline_fasta() {
        let path = std::env::temp_dir().join(format!("gm2-fastx-fasta-{}", std::process::id()));
        let mut file = File::create(&path).unwrap();
        file.write_all(b">one\nAC\nGT\n>two\nNN\n").unwrap();
        drop(file);
        let mut reader = FastxReader::open(&path, FastxFormat::Fasta).unwrap();
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"ACGT");
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"NN");
        assert!(reader.next_record().unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }
}
