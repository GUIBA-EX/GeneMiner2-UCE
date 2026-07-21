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

#[derive(Clone, Debug, Default)]
pub struct IndexProfile {
    pub recruit_probes: u64,
    pub recruit_bloom_rejected: u64,
    pub recruit_hits: u64,
    pub anchor_hit_keys: u64,
    pub anchor_occurrences: u64,
    pub exact_extensions: u64,
    pub exact_seed_bases: u64,
}

/// Four probes confined to one 64-bit word. It has no false negatives; a
/// positive is always confirmed by the exact reference hash table.
#[derive(Debug)]
struct BlockedBloom {
    blocks: Vec<u64>,
    mask: usize,
}

impl BlockedBloom {
    fn for_keys(keys: usize) -> Self {
        let blocks = keys
            .saturating_mul(12)
            .div_ceil(64)
            .max(1)
            .next_power_of_two();
        Self {
            blocks: vec![0; blocks],
            mask: blocks - 1,
        }
    }

    #[inline(always)]
    fn mix(mut value: u64) -> u64 {
        value ^= value >> 30;
        value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value ^= value >> 27;
        value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    #[inline(always)]
    fn bit_mask(hash: u64) -> u64 {
        (1_u64 << ((hash >> 16) & 63))
            | (1_u64 << ((hash >> 28) & 63))
            | (1_u64 << ((hash >> 40) & 63))
            | (1_u64 << ((hash >> 52) & 63))
    }

    #[inline(always)]
    fn insert_hash(&mut self, hash: u64) {
        let block = hash as usize & self.mask;
        self.blocks[block] |= Self::bit_mask(hash);
    }

    #[inline(always)]
    fn may_contain_hash(&self, hash: u64) -> bool {
        let block = hash as usize & self.mask;
        self.blocks[block] & Self::bit_mask(hash) == Self::bit_mask(hash)
    }
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

#[derive(Debug, Default)]
pub struct ReadEvidenceScratch {
    orientation_events: Vec<u8>,
    windows: usize,
    best_exact: Vec<Option<ExactSeed>>,
    covered_small: Vec<((usize, u32, isize), usize)>,
    covered_large: AHashMap<(usize, u32, isize), usize>,
    large_coverage: bool,
    locus_generation: Vec<u32>,
    locus_slots: Vec<usize>,
    generation: u32,
}

impl ReadEvidenceScratch {
    const SMALL_COVERAGE_LIMIT: usize = 16;

    fn reset(&mut self, candidates: &[LocusId], windows: usize, locus_count: usize) {
        self.windows = windows;
        self.orientation_events
            .resize(candidates.len() * windows, 0);
        self.orientation_events.fill(0);
        self.best_exact.resize(candidates.len(), None);
        self.best_exact.fill(None);
        self.covered_small.clear();
        self.covered_large.clear();
        self.large_coverage = false;
        self.locus_generation.resize(locus_count, 0);
        self.locus_slots.resize(locus_count, 0);
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.locus_generation.fill(0);
            self.generation = 1;
        }
        for (slot, &locus) in candidates.iter().enumerate() {
            let locus = locus as usize;
            self.locus_generation[locus] = self.generation;
            self.locus_slots[locus] = slot;
        }
    }

    #[inline(always)]
    fn slot_for(&self, locus: LocusId) -> Option<usize> {
        let locus = locus as usize;
        (self.locus_generation[locus] == self.generation).then_some(self.locus_slots[locus])
    }

    pub fn orientation(&self, candidate: usize) -> &[u8] {
        let start = candidate * self.windows;
        &self.orientation_events[start..start + self.windows]
    }

    pub fn best(&self, candidate: usize) -> Option<ExactSeed> {
        self.best_exact[candidate]
    }

    fn covered_end(&self, key: (usize, u32, isize)) -> Option<usize> {
        if self.large_coverage {
            self.covered_large.get(&key).copied()
        } else {
            self.covered_small
                .iter()
                .find_map(|&(candidate, end)| (candidate == key).then_some(end))
        }
    }

    fn record_coverage(&mut self, key: (usize, u32, isize), end: usize) {
        if self.large_coverage {
            self.covered_large.insert(key, end);
            return;
        }
        if let Some(entry) = self
            .covered_small
            .iter_mut()
            .find(|(candidate, _)| *candidate == key)
        {
            entry.1 = end;
            return;
        }
        if self.covered_small.len() < Self::SMALL_COVERAGE_LIMIT {
            self.covered_small.push((key, end));
            return;
        }
        self.covered_large.extend(self.covered_small.drain(..));
        self.large_coverage = true;
        self.covered_large.insert(key, end);
    }
}

/// Reusable, generation-stamped candidate set for one sample. It preserves
/// the previous recruit order while eliminating per-fragment allocation and
/// linear duplicate checks.
#[derive(Debug, Default)]
pub struct RecruitScratch {
    loci: Vec<LocusId>,
    seen_generation: Vec<u32>,
    generation: u32,
}

impl RecruitScratch {
    pub fn begin(&mut self, locus_count: usize) {
        self.loci.clear();
        self.seen_generation.resize(locus_count, 0);
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.seen_generation.fill(0);
            self.generation = 1;
        }
    }

    #[inline(always)]
    fn insert(&mut self, locus: LocusId) {
        let locus_index = locus as usize;
        if self.seen_generation[locus_index] != self.generation {
            self.seen_generation[locus_index] = self.generation;
            self.loci.push(locus);
        }
    }

    pub fn loci(&self) -> &[LocusId] {
        &self.loci
    }

    pub fn sort(&mut self) {
        self.loci.sort_unstable();
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

#[derive(Debug)]
enum RecruitIndex {
    Short(AHashMap<u64, LocusHits>),
    Long(AHashMap<u128, LocusHits>),
}

enum AnchorIndex {
    Short(AHashMap<u64, AnchorHits>),
    Long(AHashMap<u128, AnchorHits>),
}

impl RecruitIndex {
    fn new(k: usize) -> Self {
        if k <= 32 {
            Self::Short(AHashMap::new())
        } else {
            Self::Long(AHashMap::new())
        }
    }

    fn clear(&mut self) {
        match self {
            Self::Short(index) => index.clear(),
            Self::Long(index) => index.clear(),
        }
    }

    fn insert_sequence(&mut self, sequence: &[u8], k: usize, locus: LocusId) {
        match self {
            Self::Short(index) => scan_kmers_u64(sequence, k, 1, true, |key, _| {
                index
                    .entry(key)
                    .and_modify(|hits| hits.insert(locus))
                    .or_insert(LocusHits::One(locus));
            }),
            Self::Long(index) => scan_kmers(sequence, k, 1, true, |key, _| {
                index
                    .entry(key)
                    .and_modify(|hits| hits.insert(locus))
                    .or_insert(LocusHits::One(locus));
            }),
        }
    }

    fn bloom(&self) -> BlockedBloom {
        match self {
            Self::Short(index) => {
                let mut bloom = BlockedBloom::for_keys(index.len());
                for &key in index.keys() {
                    bloom.insert_hash(BlockedBloom::mix(key));
                }
                bloom
            }
            Self::Long(index) => {
                let mut bloom = BlockedBloom::for_keys(index.len());
                for &key in index.keys() {
                    bloom.insert_hash(BlockedBloom::mix(key as u64 ^ (key >> 64) as u64));
                }
                bloom
            }
        }
    }

    fn scan(
        &self,
        sequence: &[u8],
        k: usize,
        step: usize,
        bloom: &BlockedBloom,
        mut visit: impl FnMut(Option<&LocusHits>, bool),
    ) {
        match self {
            Self::Short(index) => scan_kmers_u64(sequence, k, step, true, |key, _| {
                let may_contain = bloom.may_contain_hash(BlockedBloom::mix(key));
                visit(may_contain.then(|| index.get(&key)).flatten(), !may_contain);
            }),
            Self::Long(index) => scan_kmers(sequence, k, step, true, |key, _| {
                let may_contain =
                    bloom.may_contain_hash(BlockedBloom::mix(key as u64 ^ (key >> 64) as u64));
                visit(may_contain.then(|| index.get(&key)).flatten(), !may_contain);
            }),
        }
    }
}

impl AnchorIndex {
    fn new(k: usize) -> Self {
        if k <= 32 {
            Self::Short(AHashMap::new())
        } else {
            Self::Long(AHashMap::new())
        }
    }

    fn insert_sequence(&mut self, bases: &[u8], k: usize, locus: LocusId, sequence: u32) {
        let insert = |position: usize, hits: &mut AnchorHits| {
            hits.push(AnchorOccurrence {
                locus,
                sequence,
                position: position as u32,
            });
        };
        match self {
            Self::Short(index) => scan_kmers_u64(bases, k, 1, false, |key, position| {
                index
                    .entry(key)
                    .and_modify(|hits| insert(position, hits))
                    .or_insert(AnchorHits::One(AnchorOccurrence {
                        locus,
                        sequence,
                        position: position as u32,
                    }));
            }),
            Self::Long(index) => scan_kmers(bases, k, 1, false, |key, position| {
                index
                    .entry(key)
                    .and_modify(|hits| insert(position, hits))
                    .or_insert(AnchorHits::One(AnchorOccurrence {
                        locus,
                        sequence,
                        position: position as u32,
                    }));
            }),
        }
    }

    fn scan(&self, sequence: &[u8], k: usize, mut visit: impl FnMut(&[AnchorOccurrence], usize)) {
        match self {
            Self::Short(index) => scan_kmers_u64(sequence, k, 1, false, |key, position| {
                if let Some(hits) = index.get(&key) {
                    visit(hits.values(), position);
                }
            }),
            Self::Long(index) => scan_kmers(sequence, k, 1, false, |key, position| {
                if let Some(hits) = index.get(&key) {
                    visit(hits.values(), position);
                }
            }),
        }
    }

    fn entries(&self) -> usize {
        match self {
            Self::Short(index) => index.values().map(|hits| hits.values().len()).sum(),
            Self::Long(index) => index.values().map(|hits| hits.values().len()).sum(),
        }
    }
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

pub struct UceIndex {
    pub k: usize,
    pub run_k: usize,
    pub loci: Vec<Locus>,
    pub references: Vec<OrientedReference>,
    recruit: RecruitIndex,
    recruit_bloom: BlockedBloom,
    anchors: AnchorIndex,
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

pub fn scan_kmers_u64(
    sequence: &[u8],
    k: usize,
    step: usize,
    canonical: bool,
    mut visit: impl FnMut(u64, usize),
) {
    debug_assert!((1..=32).contains(&k));
    if k == 0 || k > 32 || sequence.len() < k {
        return;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * k)) - 1
    };
    let reverse_shift = 2 * (k - 1);
    let tail = sequence.len() - k;
    let (mut forward, mut reverse, mut valid, mut next_probe) = (0_u64, 0_u64, 0_usize, 0_usize);
    for (end, &base) in sequence.iter().enumerate() {
        if let Some(value) = code(base) {
            forward = ((forward << 2) | value as u64) & mask;
            reverse = (reverse >> 2) | (((3 - value) as u64) << reverse_shift);
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
            recruit: RecruitIndex::new(k),
            recruit_bloom: BlockedBloom::for_keys(1),
            anchors: AnchorIndex::new(run_k),
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
                index.recruit.insert_sequence(&original, k, locus);
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
                    index.recruit.insert_sequence(&record.sequence, k, locus);
                }
            }
        }
        index.recruit_bloom = index.recruit.bloom();
        Ok(index)
    }

    fn add_oriented(&mut self, locus: LocusId, strand: u8, bases: Vec<u8>) {
        let sequence = self.references.len() as u32;
        self.anchors
            .insert_sequence(&bases, self.run_k, locus, sequence);
        self.references.push(OrientedReference {
            locus,
            strand,
            bases,
        });
    }

    pub fn recruit(
        &self,
        sequence: &[u8],
        step: usize,
        hits: &mut RecruitScratch,
        profile: Option<&mut IndexProfile>,
    ) {
        let mut profile = profile;
        self.recruit.scan(
            sequence,
            self.k,
            step,
            &self.recruit_bloom,
            |loci, bloom_rejected| {
                if let Some(profile) = profile.as_deref_mut() {
                    profile.recruit_probes += 1;
                    profile.recruit_bloom_rejected += u64::from(bloom_rejected);
                    profile.recruit_hits += u64::from(loci.is_some());
                }
                let Some(loci) = loci else { return };
                for &locus in loci.values() {
                    hits.insert(locus);
                }
            },
        );
    }

    pub fn orientation_events(&self, sequence: &[u8], candidates: &[LocusId]) -> Vec<Vec<u8>> {
        let windows = sequence.len().saturating_sub(self.run_k).saturating_add(1);
        let mut result = vec![vec![0_u8; windows]; candidates.len()];
        if !valid_dna(sequence) || windows == 0 {
            return result;
        }
        self.anchors
            .scan(sequence, self.run_k, |entries, position| {
                for occurrence in entries {
                    let locus = occurrence.locus;
                    let mask = self.references[occurrence.sequence as usize].strand;
                    if let Ok(i) = candidates.binary_search(&locus) {
                        result[i][position] |= mask;
                    }
                }
            });
        result
    }

    /// Collects run-k orientation events and the best exact seed for every
    /// candidate locus in one anchor-index traversal.
    pub fn read_evidence(
        &self,
        read: &[u8],
        candidates: &[LocusId],
        result: &mut ReadEvidenceScratch,
        profile: Option<&mut IndexProfile>,
    ) {
        let windows = read.len().saturating_sub(self.run_k).saturating_add(1);
        result.reset(candidates, windows, self.loci.len());
        if !valid_dna(read) || windows == 0 || candidates.is_empty() {
            return;
        }
        let mut profile = profile;
        self.anchors
            .scan(read, self.run_k, |occurrences, read_pos| {
                if let Some(profile) = profile.as_deref_mut() {
                    profile.anchor_hit_keys += 1;
                    profile.anchor_occurrences += occurrences.len() as u64;
                }
                for occurrence in occurrences {
                    let Some(slot) = result.slot_for(occurrence.locus) else {
                        continue;
                    };
                    let reference = &self.references[occurrence.sequence as usize];
                    result.orientation_events[slot * windows + read_pos] |= reference.strand;
                    let ref_pos = occurrence.position as usize;
                    let diagonal = ref_pos as isize - read_pos as isize;
                    let coverage_key = (slot, occurrence.sequence, diagonal);
                    if result
                        .covered_end(coverage_key)
                        .is_some_and(|end| read_pos < end)
                    {
                        continue;
                    }
                    let mut left = 0_usize;
                    while left < read_pos
                        && left < ref_pos
                        && read[read_pos - left - 1] == reference.bases[ref_pos - left - 1]
                    {
                        left += 1;
                    }
                    let mut right = self.run_k;
                    while read_pos + right < read.len()
                        && ref_pos + right < reference.bases.len()
                        && read[read_pos + right] == reference.bases[ref_pos + right]
                    {
                        right += 1;
                    }
                    if let Some(profile) = profile.as_deref_mut() {
                        profile.exact_extensions += 1;
                        profile.exact_seed_bases += (left + right) as u64;
                    }
                    let candidate = ExactSeed {
                        sequence: occurrence.sequence,
                        read_start: (read_pos - left).min(u16::MAX as usize) as u16,
                        read_end: (read_pos + right).min(u16::MAX as usize) as u16,
                        reference_start: (ref_pos - left).min(u32::MAX as usize) as u32,
                        reference_end: (ref_pos + right).min(u32::MAX as usize) as u32,
                    };
                    let best = &mut result.best_exact[slot];
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
                        *best = Some(candidate);
                    }
                    result.record_coverage(coverage_key, read_pos + right);
                }
            });
    }

    pub fn best_exact(&self, read: &[u8], locus: LocusId) -> Option<ExactSeed> {
        if !valid_dna(read) || read.len() < self.run_k {
            return None;
        }
        let mut best = None::<ExactSeed>;
        let mut covered: AHashMap<(u32, isize), usize> = AHashMap::new();
        self.anchors
            .scan(read, self.run_k, |occurrences, read_pos| {
                for occurrence in occurrences.iter().filter(|entry| entry.locus == locus) {
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
        self.anchors.entries()
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
    fn evidence_scratch_keeps_slots_beyond_u16_range() {
        let candidates: Vec<LocusId> = (0..=u16::MAX as u32 + 1).collect();
        let mut scratch = ReadEvidenceScratch::default();
        scratch.reset(&candidates, 1, candidates.len());
        assert_eq!(
            scratch.slot_for(u16::MAX as u32 + 1),
            Some(u16::MAX as usize + 1)
        );
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
