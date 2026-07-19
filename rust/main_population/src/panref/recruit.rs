use super::bait::BaitCatalog;
use super::bait_index::BaitIndex;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::{BTreeSet, HashMap, VecDeque};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::{Child, Command, Stdio};

// UCE libraries contain many non-target reads, so recruitment must be a cheap
// exact test first.  A 23-mer leaves multiple independent anchors in a 150 bp
// read even when bait and sample are moderately divergent.
const SEED_LEN: usize = 23;
const MAX_OPEN_LOCUS_FILES: usize = 96;
const RECRUIT_BATCH_PAIRS: usize = 4096;

#[derive(Debug)]
pub(crate) struct RecruitStats {
    pub(crate) per_locus: Vec<u64>,
    pub(crate) strong_pairs: u64,
    pub(crate) rescued_pairs: u64,
    pub(crate) ambiguous_pairs: u64,
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct PairAssignment {
    pub(crate) shared_loci: Vec<u32>,
    pub(crate) rescued_loci: Vec<u32>,
    pub(crate) ambiguous: bool,
}

impl PairAssignment {
    pub(crate) fn all_loci(&self) -> impl Iterator<Item = u32> + '_ {
        self.shared_loci
            .iter()
            .chain(self.rescued_loci.iter())
            .copied()
    }
}

#[derive(Clone, Debug)]
struct FastqRecord {
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

struct FastqInput {
    reader: Box<dyn BufRead>,
    child: Option<Child>,
}

impl FastqInput {
    fn finish(mut self) -> Result<(), String> {
        drop(self.reader);
        if let Some(mut child) = self.child.take() {
            let status = child.wait().map_err(|error| error.to_string())?;
            if !status.success() {
                return Err(format!("pigz decompression failed with {status}"));
            }
        }
        Ok(())
    }
}

struct LocusWriters {
    directory: std::path::PathBuf,
    writers: HashMap<u32, BufWriter<File>>,
    order: VecDeque<u32>,
}

impl LocusWriters {
    fn new(directory: &Path) -> Result<Self, String> {
        fs::create_dir_all(directory).map_err(|e| e.to_string())?;
        Ok(Self {
            directory: directory.to_path_buf(),
            writers: HashMap::new(),
            order: VecDeque::new(),
        })
    }

    fn write_pair(
        &mut self,
        locus: u32,
        sample_id: &str,
        pair_id: u64,
        first: &FastqRecord,
        second: Option<&FastqRecord>,
    ) -> Result<(), String> {
        if !self.writers.contains_key(&locus) {
            if self.writers.len() == MAX_OPEN_LOCUS_FILES {
                let old = self.order.pop_front().expect("nonempty writer cache");
                self.writers
                    .remove(&old)
                    .expect("cached writer")
                    .flush()
                    .map_err(|e| e.to_string())?;
            }
            let path = self
                .directory
                .join(format!("locus_{locus:05}.interleaved.fq"));
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .map_err(|e| e.to_string())?;
            self.writers.insert(locus, BufWriter::new(file));
            self.order.push_back(locus);
        }
        let writer = self.writers.get_mut(&locus).expect("writer just inserted");
        if let Some(second) = second {
            write_record(writer, sample_id, pair_id, 1, first)?;
            write_record(writer, sample_id, pair_id, 2, second)?;
        } else {
            write_record(writer, sample_id, pair_id, 0, first)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(), String> {
        for (_, mut writer) in self.writers.drain() {
            writer.flush().map_err(|e| e.to_string())?;
        }
        Ok(())
    }
}

fn write_record(
    writer: &mut BufWriter<File>,
    sample_id: &str,
    pair_id: u64,
    mate: u8,
    record: &FastqRecord,
) -> Result<(), String> {
    writeln!(writer, "@{sample_id}:panref_{pair_id}/{mate}").map_err(|e| e.to_string())?;
    writer
        .write_all(&record.sequence)
        .map_err(|e| e.to_string())?;
    writer.write_all(b"\n+\n").map_err(|e| e.to_string())?;
    writer
        .write_all(&record.quality)
        .map_err(|e| e.to_string())?;
    writer.write_all(b"\n").map_err(|e| e.to_string())
}

pub(crate) fn recruit_pair(index: &BaitIndex, read1: &[u8], read2: &[u8]) -> PairAssignment {
    let left = candidate_loci(index, read1);
    let right = candidate_loci(index, read2);
    let shared_loci: Vec<u32> = left.intersection(&right).copied().collect();
    let rescued_loci = if shared_loci.is_empty() {
        if left.len() == 1 && right.is_empty() {
            left.iter().copied().collect()
        } else if right.len() == 1 && left.is_empty() {
            right.iter().copied().collect()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    // A pair shared by multiple UCE loci is ambiguous, not evidence for each locus.
    let ambiguous = shared_loci.len() > 1
        || (shared_loci.is_empty()
            && rescued_loci.is_empty()
            && (!left.is_empty() || !right.is_empty()));
    PairAssignment {
        shared_loci: if shared_loci.len() == 1 {
            shared_loci
        } else {
            Vec::new()
        },
        rescued_loci,
        ambiguous,
    }
}

struct ClassifiedPair {
    pair_id: u64,
    locus: Option<u32>,
    strong: bool,
    rescued: bool,
    ambiguous: bool,
    first: FastqRecord,
    second: Option<FastqRecord>,
}

/// Stream bait-oriented, uniquely assigned read pairs without materialising a
/// per-locus FASTQ. Callers decide whether strict PE-supported pairs only or
/// single-mate rescue pairs are admissible.
pub(crate) fn stream_recruited_pairs(
    index: &BaitIndex,
    catalog: &BaitCatalog,
    read1: &Path,
    read2: &Path,
    threads: usize,
    include_rescued: bool,
    mut consume: impl FnMut(u32, bool, &[u8], Option<&[u8]>),
) -> Result<RecruitStats, String> {
    let mut first = open_fastq(read1, threads)?;
    let mut second = if read1 == read2 {
        None
    } else {
        Some(open_fastq(read2, threads)?)
    };
    let mut counts = vec![0; index.locus_count()];
    let mut strong_pairs = 0_u64;
    let mut rescued_pairs = 0_u64;
    let mut ambiguous_pairs = 0_u64;
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|error| error.to_string())?;
    let mut pair_id = 0_u64;
    loop {
        let mut batch = Vec::with_capacity(RECRUIT_BATCH_PAIRS);
        for _ in 0..RECRUIT_BATCH_PAIRS {
            let Some(record1) = read_record(first.reader.as_mut())? else {
                if let Some(reader) = second.as_mut() {
                    if read_record(reader.reader.as_mut())?.is_some() {
                        return Err("paired FASTQ files have different record counts".into());
                    }
                }
                break;
            };
            let record2 = match second.as_mut() {
                Some(reader) => Some(
                    read_record(reader.reader.as_mut())?
                        .ok_or("paired FASTQ files have different record counts")?,
                ),
                None => None,
            };
            pair_id += 1;
            batch.push((pair_id, record1, record2));
        }
        if batch.is_empty() {
            break;
        }
        let classified = pool.install(|| {
            batch
                .into_par_iter()
                .map(|(id, left, right)| classify_pair(index, catalog, id, left, right))
                .collect::<Vec<_>>()
        });
        for pair in classified {
            strong_pairs += u64::from(pair.strong);
            rescued_pairs += u64::from(pair.rescued);
            ambiguous_pairs += u64::from(pair.ambiguous);
            if let Some(locus) = pair.locus {
                if pair.strong || (include_rescued && pair.rescued) {
                    counts[locus as usize] += 1;
                    consume(
                        locus,
                        pair.strong,
                        &pair.first.sequence,
                        pair.second.as_ref().map(|x| x.sequence.as_slice()),
                    );
                }
            }
        }
    }
    first.finish()?;
    if let Some(second) = second {
        second.finish()?;
    }
    Ok(RecruitStats {
        per_locus: counts,
        strong_pairs,
        rescued_pairs,
        ambiguous_pairs,
    })
}

fn classify_pair(
    index: &BaitIndex,
    catalog: &BaitCatalog,
    pair_id: u64,
    record1: FastqRecord,
    record2: Option<FastqRecord>,
) -> ClassifiedPair {
    let assignment = recruit_pair(
        index,
        &record1.sequence,
        record2.as_ref().map_or(&[], |record| &record.sequence),
    );
    let strong = !assignment.shared_loci.is_empty();
    let rescued = !assignment.rescued_loci.is_empty();
    let ambiguous = assignment.ambiguous;
    let locus = assignment.all_loci().next();
    if let Some(locus) = locus {
        let baits = &catalog.loci[locus as usize].records;
        let first = orient_to_bait(&record1, baits);
        let second = record2.as_ref().map(|record| orient_to_bait(record, baits));
        ClassifiedPair {
            pair_id,
            locus: Some(locus),
            strong,
            rescued,
            ambiguous,
            first,
            second,
        }
    } else {
        ClassifiedPair {
            pair_id,
            locus: None,
            strong,
            rescued,
            ambiguous,
            first: record1,
            second: record2,
        }
    }
}

pub(crate) fn recruit_pairs_to_fastq(
    index: &BaitIndex,
    catalog: &BaitCatalog,
    read1: &Path,
    read2: &Path,
    output: &Path,
    sample_id: &str,
    threads: usize,
) -> Result<RecruitStats, String> {
    let mut first = open_fastq(read1, threads)?;
    let mut second = if read1 == read2 {
        None
    } else {
        Some(open_fastq(read2, threads)?)
    };
    let mut writers = LocusWriters::new(output)?;
    let mut counts = vec![0; index.locus_count()];
    let mut strong_pairs = 0_u64;
    let mut rescued_pairs = 0_u64;
    let mut ambiguous_pairs = 0_u64;
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|error| error.to_string())?;
    let mut pair_id = 0_u64;
    loop {
        let mut batch = Vec::with_capacity(RECRUIT_BATCH_PAIRS);
        for _ in 0..RECRUIT_BATCH_PAIRS {
            let Some(record1) = read_record(first.reader.as_mut())? else {
                if let Some(reader) = second.as_mut() {
                    if read_record(reader.reader.as_mut())?.is_some() {
                        return Err("paired FASTQ files have different record counts".into());
                    }
                }
                break;
            };
            let record2 = match second.as_mut() {
                Some(reader) => Some(
                    read_record(reader.reader.as_mut())?
                        .ok_or("paired FASTQ files have different record counts")?,
                ),
                None => None,
            };
            pair_id += 1;
            batch.push((pair_id, record1, record2));
        }
        if batch.is_empty() {
            break;
        }
        let classified = pool.install(|| {
            batch
                .into_par_iter()
                .map(|(id, first, second)| classify_pair(index, catalog, id, first, second))
                .collect::<Vec<_>>()
        });
        for pair in classified {
            strong_pairs += u64::from(pair.strong);
            rescued_pairs += u64::from(pair.rescued);
            ambiguous_pairs += u64::from(pair.ambiguous);
            if let Some(locus) = pair.locus {
                writers.write_pair(
                    locus,
                    sample_id,
                    pair.pair_id,
                    &pair.first,
                    pair.second.as_ref(),
                )?;
                counts[locus as usize] += 1;
            }
        }
    }
    writers.finish()?;
    first.finish()?;
    if let Some(second) = second {
        second.finish()?;
    }
    Ok(RecruitStats {
        per_locus: counts,
        strong_pairs,
        rescued_pairs,
        ambiguous_pairs,
    })
}

fn orient_to_bait(record: &FastqRecord, baits: &[Vec<u8>]) -> FastqRecord {
    let forward = bait_seed_hits(&record.sequence, baits);
    let reverse = bait_seed_hits(&reverse_complement(&record.sequence), baits);
    if reverse > forward {
        FastqRecord {
            sequence: reverse_complement(&record.sequence),
            quality: record.quality.iter().rev().copied().collect(),
        }
    } else {
        record.clone()
    }
}

fn bait_seed_hits(sequence: &[u8], baits: &[Vec<u8>]) -> usize {
    if sequence.len() < SEED_LEN {
        return 0;
    }
    baits
        .iter()
        .map(|bait| {
            sequence
                .windows(SEED_LEN)
                .step_by(SEED_LEN)
                .filter(|seed| bait.windows(SEED_LEN).any(|target| target == *seed))
                .count()
        })
        .sum()
}

fn candidate_loci(index: &BaitIndex, read: &[u8]) -> BTreeSet<u32> {
    if read.len() < SEED_LEN {
        return BTreeSet::new();
    }
    let mut loci = BTreeSet::new();
    for sequence in [read.to_ascii_uppercase(), reverse_complement(read)] {
        loci.extend(index.minimizer_loci(&sequence));
        if loci.len() > 1 {
            return loci;
        }
    }
    loci
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

fn open_fastq(path: &Path, threads: usize) -> Result<FastqInput, String> {
    let compressed = path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"));
    if compressed {
        let workers = threads.max(1).div_ceil(2);
        if let Ok(mut child) = Command::new("pigz")
            .arg("-dc")
            .arg("-p")
            .arg(workers.to_string())
            .arg(path)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
        {
            let stdout = child.stdout.take().ok_or("cannot capture pigz stdout")?;
            return Ok(FastqInput {
                reader: Box::new(BufReader::new(stdout)),
                child: Some(child),
            });
        }
    }
    let file =
        File::open(path).map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    let input: Box<dyn Read> = if compressed {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(FastqInput {
        reader: Box::new(BufReader::new(input)),
        child: None,
    })
}

fn read_record(reader: &mut dyn BufRead) -> Result<Option<FastqRecord>, String> {
    let mut header = String::new();
    if reader
        .read_line(&mut header)
        .map_err(|error| error.to_string())?
        == 0
    {
        return Ok(None);
    }
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    for line in [&mut sequence, &mut plus, &mut quality] {
        if reader.read_line(line).map_err(|error| error.to_string())? == 0 {
            return Err("truncated FASTQ record".into());
        }
    }
    if !header.starts_with('@') || !plus.starts_with('+') {
        return Err("invalid FASTQ record".into());
    }
    let sequence = sequence.trim().as_bytes().to_ascii_uppercase();
    let quality = quality.trim().as_bytes().to_vec();
    if sequence.len() != quality.len() {
        return Err("FASTQ sequence and quality lengths differ".into());
    }
    Ok(Some(FastqRecord { sequence, quality }))
}

#[cfg(test)]
mod tests {
    use super::{orient_to_bait, recruit_pair, FastqRecord, PairAssignment};
    use crate::panref::bait_index::{BaitIndex, BaitRecord};
    #[test]
    fn paired_evidence_beats_single_mate_rescue() {
        let index = BaitIndex::build(&[BaitRecord {
            locus: 0,
            sequence: b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),
        }])
        .unwrap();
        assert_eq!(
            recruit_pair(
                &index,
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            ),
            PairAssignment {
                shared_loci: vec![0],
                rescued_loci: vec![],
                ambiguous: false
            }
        );
        assert_eq!(
            recruit_pair(
                &index,
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                b"ACGTACGTACGTACGTACGTACGTACGTACG"
            ),
            PairAssignment {
                shared_loci: vec![],
                rescued_loci: vec![0],
                ambiguous: false
            }
        );
    }

    #[test]
    fn bait_orientation_normalizes_reverse_reads() {
        let bait = vec![b"ACGTCCTAGGATTCGACCTGTAAGCGTACCAA".to_vec()];
        let record = FastqRecord {
            sequence: super::reverse_complement(&bait[0]),
            quality: vec![b'F'; bait[0].len()],
        };
        let oriented = orient_to_bait(&record, &bait);
        assert_eq!(oriented.sequence, bait[0]);
        assert_eq!(
            oriented.quality,
            record.quality.iter().rev().copied().collect::<Vec<_>>()
        );
    }
}
