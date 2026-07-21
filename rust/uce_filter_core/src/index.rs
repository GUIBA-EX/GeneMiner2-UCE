use crate::model::{Locus, LocusId};
use ahash::AHashMap;
use gm2_tools::fastx::{FastxFormat, FastxReader};
use std::fs;
use std::path::{Path, PathBuf};

const INVALID: u8 = u8::MAX;

const fn base_table() -> [u8; 256] {
    let mut table = [INVALID; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
}

pub const BASE_CODE: [u8; 256] = base_table();

#[inline(always)]
pub fn code(base: u8) -> Option<u8> {
    let value = BASE_CODE[base as usize];
    (value != INVALID).then_some(value)
}

pub fn valid_dna(sequence: &[u8]) -> bool {
    sequence.iter().all(|&base| code(base).is_some())
}

#[derive(Clone, Copy, Debug)]
pub struct AnchorOccurrence {
    pub locus: LocusId,
    pub sequence: u32,
    pub position: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExactSeed {
    pub sequence: u32,
    pub read_start: u16,
    pub read_end: u16,
    pub reference_start: u32,
    pub reference_end: u32,
}

impl ExactSeed {
    pub fn len(self) -> usize {
        self.read_end as usize - self.read_start as usize
    }

    pub fn is_empty(self) -> bool {
        self.read_start == self.read_end
    }
}

#[derive(Debug)]
enum LocusHits {
    One(LocusId),
    Many(Vec<LocusId>),
}

impl LocusHits {
    fn insert(&mut self, locus: LocusId) {
        match self {
            Self::One(existing) if *existing != locus => {
                *self = Self::Many(vec![*existing, locus]);
            }
            Self::Many(values) if !values.contains(&locus) => values.push(locus),
            _ => {}
        }
    }

    fn values(&self) -> &[LocusId] {
        match self {
            Self::One(value) => std::slice::from_ref(value),
            Self::Many(values) => values,
        }
    }
}

#[derive(Debug)]
enum AnchorHits {
    One(AnchorOccurrence),
    Many(Vec<AnchorOccurrence>),
}

impl AnchorHits {
    fn push(&mut self, occurrence: AnchorOccurrence) {
        match self {
            Self::One(existing) => {
                *self = Self::Many(vec![*existing, occurrence]);
            }
            Self::Many(values) => values.push(occurrence),
        }
    }

    fn values(&self) -> &[AnchorOccurrence] {
        match self {
            Self::One(value) => std::slice::from_ref(value),
            Self::Many(values) => values,
        }
    }
}

#[derive(Clone, Debug)]
pub struct OrientedReference {
    pub locus: LocusId,
    pub strand: u8,
    pub bases: Vec<u8>,
}

#[derive(Debug)]
pub struct UceIndex {
    pub k: usize,
    pub run_k: usize,
    pub loci: Vec<Locus>,
    pub references: Vec<OrientedReference>,
    recruit: AHashMap<u128, LocusHits>,
    anchors: AHashMap<u128, AnchorHits>,
}

fn stripped_extension(path: &Path) -> Option<String> {
    let base = if path
        .extension()
        .and_then(|v| v.to_str())
        .is_some_and(|v| v.eq_ignore_ascii_case("gz"))
    {
        PathBuf::from(path.file_stem()?)
    } else {
        path.to_path_buf()
    };
    base.extension()?.to_str().map(|v| v.to_ascii_lowercase())
}

fn reference_name(path: &Path) -> Result<String, String> {
    let base = if path
        .extension()
        .and_then(|v| v.to_str())
        .is_some_and(|v| v.eq_ignore_ascii_case("gz"))
    {
        PathBuf::from(
            path.file_stem()
                .ok_or_else(|| "invalid reference path".to_string())?,
        )
    } else {
        path.to_path_buf()
    };
    base.file_stem()
        .and_then(|v| v.to_str())
        .map(str::to_owned)
        .ok_or_else(|| format!("invalid reference name: {}", path.display()))
}

fn reference_paths(reference: &Path) -> Result<Vec<PathBuf>, String> {
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if path.is_file()
                && matches!(
                    stripped_extension(&path).as_deref(),
                    Some("fa" | "fas" | "fasta")
                )
            {
                paths.push(path);
            }
        }
    } else if reference.is_file() {
        paths.push(reference.to_path_buf());
    }
    paths.sort();
    if paths.is_empty() {
        return Err("no reference FASTA file found".to_string());
    }
    Ok(paths)
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match code(base) {
            Some(value) => b"TGCA"[value as usize],
            None => base,
        })
        .collect()
}

pub fn scan_kmers(
    sequence: &[u8],
    k: usize,
    step: usize,
    canonical: bool,
    mut visit: impl FnMut(u128, usize),
) {
    if k == 0 || k > 64 || sequence.len() < k {
        return;
    }
    let mask = if k == 64 {
        u128::MAX
    } else {
        (1_u128 << (2 * k)) - 1
    };
    let reverse_shift = 2 * (k - 1);
    let tail = sequence.len() - k;
    let (mut forward, mut reverse, mut valid, mut next_probe) = (0_u128, 0_u128, 0_usize, 0_usize);
    for (end, &base) in sequence.iter().enumerate() {
        if let Some(value) = code(base) {
            forward = ((forward << 2) | value as u128) & mask;
            reverse = (reverse >> 2) | (((3 - value) as u128) << reverse_shift);
            valid += 1;
        } else {
            forward = 0;
            reverse = 0;
            valid = 0;
        }
        if end + 1 < k {
            continue;
        }
        let start = end + 1 - k;
        let sampled = start == next_probe;
        if sampled {
            next_probe = next_probe.saturating_add(step.max(1));
        }
        if valid >= k && (sampled || start == tail) {
            visit(
                if canonical {
                    forward.min(reverse)
                } else {
                    forward
                },
                start,
            );
        }
    }
}

impl UceIndex {
    pub fn build(reference: &Path, k: usize) -> Result<Self, String> {
        Self::build_split(reference, reference, k)
    }

    pub fn build_split(
        recruit_reference: &Path,
        verify_reference: &Path,
        k: usize,
    ) -> Result<Self, String> {
        if !(1..=64).contains(&k) {
            return Err("UCEFilter currently supports k-mer sizes 1..=64".to_string());
        }
        let run_k = std::cmp::max(k / 2, k.saturating_sub(13)) | 1;
        let mut index = Self {
            k,
            run_k,
            loci: Vec::new(),
            references: Vec::new(),
            recruit: AHashMap::new(),
            anchors: AHashMap::new(),
        };
        for path in reference_paths(verify_reference)? {
            let locus = index.loci.len() as LocusId;
            let name = reference_name(&path)?;
            let mut reader =
                FastxReader::open(&path, FastxFormat::Fasta).map_err(|e| e.to_string())?;
            let mut originals = Vec::new();
            while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
                if record.sequence.len() >= k {
                    originals.push(record.sequence);
                }
            }
            let max_len = originals.iter().map(Vec::len).max().unwrap_or(0) as f64;
            let effective_length = if originals.is_empty() {
                0.0
            } else {
                (max_len * ((originals.len() as f64).log10() + 1.0)).trunc()
            };
            index.loci.push(Locus {
                name,
                effective_length,
            });
            for original in originals {
                scan_kmers(&original, k, 1, true, |key, _| {
                    index
                        .recruit
                        .entry(key)
                        .and_modify(|hits| hits.insert(locus))
                        .or_insert(LocusHits::One(locus));
                });
                index.add_oriented(locus, 1, original.clone());
                index.add_oriented(locus, 2, reverse_complement(&original));
            }
        }
        if recruit_reference != verify_reference {
            index.recruit.clear();
            let locus_by_name: AHashMap<_, _> = index
                .loci
                .iter()
                .enumerate()
                .map(|(id, locus)| (locus.name.clone(), id as LocusId))
                .collect();
            for path in reference_paths(recruit_reference)? {
                let name = reference_name(&path)?;
                let Some(&locus) = locus_by_name.get(&name) else {
                    continue;
                };
                let mut reader =
                    FastxReader::open(&path, FastxFormat::Fasta).map_err(|e| e.to_string())?;
                while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
                    scan_kmers(&record.sequence, k, 1, true, |key, _| {
                        index
                            .recruit
                            .entry(key)
                            .and_modify(|hits| hits.insert(locus))
                            .or_insert(LocusHits::One(locus));
                    });
                }
            }
        }
        Ok(index)
    }

    fn add_oriented(&mut self, locus: LocusId, strand: u8, bases: Vec<u8>) {
        let sequence = self.references.len() as u32;
        scan_kmers(&bases, self.run_k, 1, false, |key, position| {
            let occurrence = AnchorOccurrence {
                locus,
                sequence,
                position: position as u32,
            };
            self.anchors
                .entry(key)
                .and_modify(|hits| hits.push(occurrence))
                .or_insert(AnchorHits::One(occurrence));
        });
        self.references.push(OrientedReference {
            locus,
            strand,
            bases,
        });
    }

    pub fn recruit(&self, sequence: &[u8], step: usize, hits: &mut Vec<LocusId>) {
        scan_kmers(sequence, self.k, step, true, |key, _| {
            if let Some(loci) = self.recruit.get(&key) {
                for &locus in loci.values() {
                    if !hits.contains(&locus) {
                        hits.push(locus);
                    }
                }
            }
        });
    }

    pub fn orientation_events(&self, sequence: &[u8], candidates: &[LocusId]) -> Vec<Vec<u8>> {
        let windows = sequence.len().saturating_sub(self.run_k).saturating_add(1);
        let mut result = vec![vec![0_u8; windows]; candidates.len()];
        if !valid_dna(sequence) || windows == 0 {
            return result;
        }
        scan_kmers(sequence, self.run_k, 1, false, |key, position| {
            if let Some(entries) = self.anchors.get(&key) {
                for occurrence in entries.values() {
                    let locus = occurrence.locus;
                    let mask = self.references[occurrence.sequence as usize].strand;
                    if let Some(i) = candidates.iter().position(|&value| value == locus) {
                        result[i][position] |= mask;
                    }
                }
            }
        });
        result
    }

    pub fn best_exact(&self, read: &[u8], locus: LocusId) -> Option<ExactSeed> {
        if !valid_dna(read) || read.len() < self.run_k {
            return None;
        }
        let mut best = None::<ExactSeed>;
        let mut covered: AHashMap<(u32, isize), usize> = AHashMap::new();
        scan_kmers(read, self.run_k, 1, false, |key, read_pos| {
            let Some(occurrences) = self.anchors.get(&key) else {
                return;
            };
            for occurrence in occurrences
                .values()
                .iter()
                .filter(|entry| entry.locus == locus)
            {
                let reference = &self.references[occurrence.sequence as usize].bases;
                let ref_pos = occurrence.position as usize;
                let diagonal = ref_pos as isize - read_pos as isize;
                if covered
                    .get(&(occurrence.sequence, diagonal))
                    .is_some_and(|&end| read_pos < end)
                {
                    continue;
                }
                let mut left = 0_usize;
                while left < read_pos
                    && left < ref_pos
                    && read[read_pos - left - 1] == reference[ref_pos - left - 1]
                {
                    left += 1;
                }
                let mut right = self.run_k;
                while read_pos + right < read.len()
                    && ref_pos + right < reference.len()
                    && read[read_pos + right] == reference[ref_pos + right]
                {
                    right += 1;
                }
                let candidate = ExactSeed {
                    sequence: occurrence.sequence,
                    read_start: (read_pos - left).min(u16::MAX as usize) as u16,
                    read_end: (read_pos + right).min(u16::MAX as usize) as u16,
                    reference_start: (ref_pos - left).min(u32::MAX as usize) as u32,
                    reference_end: (ref_pos + right).min(u32::MAX as usize) as u32,
                };
                if best.is_none_or(|current| {
                    candidate.len() > current.len()
                        || (candidate.len() == current.len()
                            && (
                                candidate.sequence,
                                candidate.reference_start,
                                candidate.read_start,
                            ) < (
                                current.sequence,
                                current.reference_start,
                                current.read_start,
                            ))
                }) {
                    best = Some(candidate);
                }
                covered.insert((occurrence.sequence, diagonal), read_pos + right);
            }
        });
        best
    }

    pub fn max_exact(&self, read: &[u8], locus: LocusId) -> usize {
        self.best_exact(read, locus).map_or(0, ExactSeed::len)
    }

    pub fn anchor_entries(&self) -> usize {
        self.anchors.values().map(|hits| hits.values().len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEST_DIR: AtomicU64 = AtomicU64::new(0);

    fn test_index(k: usize, reference: &[u8]) -> UceIndex {
        let root = std::env::temp_dir().join(format!(
            "uce-filter-index-{}-{k}-{}",
            std::process::id(),
            NEXT_TEST_DIR.fetch_add(1, Ordering::Relaxed),
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let mut out = fs::File::create(root.join("locus.fa")).unwrap();
        out.write_all(b">ref\n").unwrap();
        out.write_all(reference).unwrap();
        out.write_all(b"\n").unwrap();
        drop(out);
        let index = UceIndex::build(&root, k).unwrap();
        fs::remove_dir_all(root).unwrap();
        index
    }

    fn brute_hit(read: &[u8], reference: &[u8], k: usize) -> bool {
        if !valid_dna(read) || read.len() < k || reference.len() < k {
            return false;
        }
        let reverse = reverse_complement(reference);
        read.windows(k).any(|window| {
            reference.windows(k).any(|candidate| candidate == window)
                || reverse.windows(k).any(|candidate| candidate == window)
        })
    }

    #[test]
    fn maximum_exact_reproduces_all_thresholds() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let read = b"TTTTGCAACGATTCGGTACCAGGG";
        let index = test_index(16, reference);
        let maximum = index.max_exact(read, 0);
        for k in [1, 5, 16, 17, 20, 24] {
            assert_eq!(
                brute_hit(read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }

    #[test]
    fn ambiguous_read_has_no_legacy_exact_match() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let index = test_index(16, reference);
        assert_eq!(index.max_exact(b"ACGTTGCAACGATTCNGTACC", 0), 0);
    }

    #[test]
    fn maximum_exact_handles_long_and_reverse_strand_thresholds() {
        let reference =
            b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTTAGCTACGATGCTAGGCTTACCGATGGCATTCG";
        let read = reverse_complement(&reference[8..80]);
        let index = test_index(16, reference);
        let maximum = index.max_exact(&read, 0);
        for k in [16, 31, 32, 33, 63, 64, 67] {
            assert_eq!(
                brute_hit(&read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }

    #[test]
    fn maximum_exact_restarts_after_a_mismatch_on_the_same_diagonal() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let mut read = reference[4..49].to_vec();
        read[23] = if read[23] == b'A' { b'C' } else { b'A' };
        let index = test_index(16, reference);
        let maximum = index.max_exact(&read, 0);
        for k in [16, 20, 23, 24, 31] {
            assert_eq!(
                brute_hit(&read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }
}
