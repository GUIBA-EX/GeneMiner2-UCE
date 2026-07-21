use gm2_tools::fastx::FastxRecord;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub type FragmentId = u32;
pub type LocusId = u32;

#[derive(Clone, Debug)]
pub struct Fragment {
    pub ordinal: u64,
    pub r1: FastxRecord,
    pub r2: FastxRecord,
}

impl Fragment {
    pub fn bases(&self) -> u64 {
        (self.r1.sequence.len() + self.r2.sequence.len()) as u64
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Candidate {
    pub fragment_id: FragmentId,
    pub ordinal: u64,
    pub fragment_bases: u32,
    pub max_exact: u16,
    /// Normalized coverage of the verification reference in 64 bins.
    pub covered_bins: u64,
    /// Bit 0/1: exact evidence extends the left/right end of a contig reference.
    pub terminal_mask: u8,
    /// Maximum exact-seed overhang beyond the left/right reference edge.
    pub left_extension: u16,
    pub right_extension: u16,
    /// Number of mates with an exact seed on this locus.
    pub aligned_mates: u8,
    /// Number of loci to which this fragment passed run-k verification.
    pub locus_count: u16,
}

#[derive(Clone, Debug)]
pub struct Locus {
    pub name: String,
    pub effective_length: f64,
}

pub struct FragmentBank {
    memory: Vec<Fragment>,
    memory_bytes: u64,
    memory_limit: u64,
    spill: Option<SpillStore>,
    spill_path: PathBuf,
    total: usize,
    sequence_bases: u64,
}

struct SpillStore {
    writer: Option<BufWriter<File>>,
    count: usize,
    encoded_bytes: u64,
}

impl FragmentBank {
    pub fn new(memory_limit: u64, spill_path: PathBuf) -> Self {
        Self {
            memory: Vec::new(),
            memory_bytes: 0,
            memory_limit,
            spill: None,
            spill_path,
            total: 0,
            sequence_bases: 0,
        }
    }

    pub fn insert(&mut self, fragment: Fragment) -> Result<FragmentId, String> {
        let id =
            u32::try_from(self.total).map_err(|_| "too many retained fragments".to_string())?;
        let storage_bytes = fragment_storage_bytes(&fragment);
        self.sequence_bases += fragment.bases();
        if self.spill.is_none()
            && self.memory_bytes.saturating_add(storage_bytes) <= self.memory_limit
        {
            self.memory_bytes += storage_bytes;
            self.memory.push(fragment);
        } else {
            self.ensure_spill()?;
            let spill = self.spill.as_mut().expect("spill initialized");
            let writer = spill.writer.as_mut().expect("spill writer open");
            spill.encoded_bytes += write_fragment(writer, &fragment)?;
            spill.count += 1;
        }
        self.total += 1;
        Ok(id)
    }

    fn ensure_spill(&mut self) -> Result<(), String> {
        if self.spill.is_some() {
            return Ok(());
        }
        if let Some(parent) = self.spill_path.parent() {
            fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        let writer = BufWriter::with_capacity(
            4 * 1024 * 1024,
            File::create(&self.spill_path).map_err(|e| e.to_string())?,
        );
        self.spill = Some(SpillStore {
            writer: Some(writer),
            count: 0,
            encoded_bytes: 0,
        });
        Ok(())
    }

    pub fn stream_in_order(
        &mut self,
        mut visit: impl FnMut(FragmentId, &Fragment) -> Result<(), String>,
    ) -> Result<(), String> {
        for (id, fragment) in self.memory.iter().enumerate() {
            visit(id as FragmentId, fragment)?;
        }
        let Some(spill) = self.spill.as_mut() else {
            return Ok(());
        };
        if let Some(mut writer) = spill.writer.take() {
            writer.flush().map_err(|e| e.to_string())?;
        }
        let mut reader = BufReader::with_capacity(
            4 * 1024 * 1024,
            File::open(&self.spill_path).map_err(|e| e.to_string())?,
        );
        let first_id = self.memory.len();
        for offset in 0..spill.count {
            let fragment = read_fragment(&mut reader)?;
            visit((first_id + offset) as FragmentId, &fragment)?;
        }
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.total
    }

    pub fn is_empty(&self) -> bool {
        self.total == 0
    }

    pub fn bases(&self) -> u64 {
        self.sequence_bases
    }

    pub fn memory_bytes(&self) -> u64 {
        self.memory_bytes
    }

    pub fn spill_bytes(&self) -> u64 {
        self.spill.as_ref().map_or(0, |spill| spill.encoded_bytes)
    }
}

impl Drop for FragmentBank {
    fn drop(&mut self) {
        if let Some(mut spill) = self.spill.take() {
            drop(spill.writer.take());
            let _ = fs::remove_file(&self.spill_path);
        }
    }
}

fn fragment_storage_bytes(fragment: &Fragment) -> u64 {
    fn record_bytes(record: &FastxRecord) -> u64 {
        (record.header.len() + record.sequence.len() + record.plus.len() + record.quality.len())
            as u64
    }
    8 + record_bytes(&fragment.r1) + record_bytes(&fragment.r2)
}

fn write_u32(out: &mut impl Write, value: usize) -> Result<u64, String> {
    let value = u32::try_from(value).map_err(|_| "fragment field is too large".to_string())?;
    out.write_all(&value.to_le_bytes())
        .map_err(|e| e.to_string())?;
    Ok(4)
}

fn write_bytes(out: &mut impl Write, value: &[u8]) -> Result<u64, String> {
    let mut written = write_u32(out, value.len())?;
    out.write_all(value).map_err(|e| e.to_string())?;
    written += value.len() as u64;
    Ok(written)
}

fn write_record(out: &mut impl Write, record: &FastxRecord) -> Result<u64, String> {
    Ok(write_bytes(out, &record.header)?
        + write_bytes(out, &record.sequence)?
        + write_bytes(out, &record.plus)?
        + write_bytes(out, &record.quality)?)
}

fn write_fragment(out: &mut impl Write, fragment: &Fragment) -> Result<u64, String> {
    out.write_all(&fragment.ordinal.to_le_bytes())
        .map_err(|e| e.to_string())?;
    Ok(8 + write_record(out, &fragment.r1)? + write_record(out, &fragment.r2)?)
}

fn read_u32(input: &mut impl Read) -> Result<usize, String> {
    let mut bytes = [0_u8; 4];
    input.read_exact(&mut bytes).map_err(|e| e.to_string())?;
    Ok(u32::from_le_bytes(bytes) as usize)
}

fn read_bytes(input: &mut impl Read) -> Result<Vec<u8>, String> {
    let length = read_u32(input)?;
    let mut value = vec![0_u8; length];
    input.read_exact(&mut value).map_err(|e| e.to_string())?;
    Ok(value)
}

fn read_record(input: &mut impl Read) -> Result<FastxRecord, String> {
    Ok(FastxRecord {
        header: read_bytes(input)?,
        sequence: read_bytes(input)?,
        plus: read_bytes(input)?,
        quality: read_bytes(input)?,
    })
}

fn read_fragment(input: &mut impl Read) -> Result<Fragment, String> {
    let mut ordinal = [0_u8; 8];
    input.read_exact(&mut ordinal).map_err(|e| e.to_string())?;
    Ok(Fragment {
        ordinal: u64::from_le_bytes(ordinal),
        r1: read_record(input)?,
        r2: read_record(input)?,
    })
}

pub fn default_spill_path(output: &Path) -> PathBuf {
    output.join(".uce_filter.fragments.spool")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(name: &[u8], sequence: &[u8]) -> FastxRecord {
        FastxRecord {
            header: name.to_vec(),
            sequence: sequence.to_vec(),
            plus: b"+".to_vec(),
            quality: vec![b'I'; sequence.len()],
        }
    }

    #[test]
    fn hybrid_store_streams_memory_then_spill_in_order() {
        let path = std::env::temp_dir().join(format!("uce-fragment-spill-{}", std::process::id()));
        let _ = fs::remove_file(&path);
        let mut bank = FragmentBank::new(1, path.clone());
        for ordinal in 0..3 {
            bank.insert(Fragment {
                ordinal,
                r1: record(b"@r/1", b"ACGT"),
                r2: record(b"@r/2", b"TGCA"),
            })
            .unwrap();
        }
        let mut observed = Vec::new();
        bank.stream_in_order(|id, fragment| {
            observed.push((id, fragment.ordinal, fragment.r1.sequence.clone()));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            observed,
            vec![
                (0, 0, b"ACGT".to_vec()),
                (1, 1, b"ACGT".to_vec()),
                (2, 2, b"ACGT".to_vec()),
            ]
        );
        drop(bank);
        assert!(!path.exists());
    }
}
