//! Per-locus FM indexes for exact-match evidence.
//!
//! Recruitment remains a global rolling-k-mer lookup. Only recruited loci are
//! queried here. Each locus stores its forward and reverse-complement reference
//! segments in one FM index, separated by sentinels. Runtime state is limited to
//! symbol rank bitvectors and sparse suffix-array samples: no positional anchor
//! hash table and no full suffix array or RMQ remain resident.

use crate::index::{code, OrientedReference};

const ALPHABET_SIZE: usize = 6;
const SA_SAMPLE_RATE: usize = 16;
const SENTINEL: u8 = 0;
const SEPARATOR: u8 = 5;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct MemSeed {
    pub sequence: u32,
    pub read_start: usize,
    pub read_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
}

impl MemSeed {
    pub fn len(self) -> usize {
        self.read_end - self.read_start
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct MemQueryProfile {
    pub index_queries: u64,
    pub run_windows: u64,
    pub matching_windows: u64,
    pub mem_starts: u64,
    pub mem_bases: u64,
}

#[derive(Clone, Copy, Debug)]
struct ReferenceSpan {
    text_start: u32,
    text_end: u32,
    sequence: u32,
    reference_start: u32,
    reference_end: u32,
    strand: u8,
}

#[derive(Debug, Default)]
pub(crate) struct LocusMemIndex {
    length: usize,
    spans: Vec<ReferenceSpan>,
    cumulative: [u32; ALPHABET_SIZE],
    symbol_bits: [Vec<u64>; ALPHABET_SIZE],
    symbol_prefix: [Vec<u32>; ALPHABET_SIZE],
    sample_bits: Vec<u64>,
    sample_prefix: Vec<u32>,
    sample_values: Vec<u32>,
}

impl LocusMemIndex {
    fn build(&mut self, text: &[u8]) -> Result<(), String> {
        if text.is_empty() {
            return Ok(());
        }
        let suffixes = build_suffix_array(text)?;
        self.length = text.len();
        let words = self.length.div_ceil(64);
        self.symbol_bits = std::array::from_fn(|_| vec![0_u64; words]);
        self.sample_bits = vec![0_u64; words];
        self.sample_values.clear();

        let mut frequencies = [0_u32; ALPHABET_SIZE];
        for &symbol in text {
            frequencies[symbol as usize] += 1;
        }
        let mut preceding = 0_u32;
        for (symbol, frequency) in frequencies.into_iter().enumerate() {
            self.cumulative[symbol] = preceding;
            preceding += frequency;
        }

        for (row, &position) in suffixes.iter().enumerate() {
            let symbol = if position == 0 {
                SENTINEL
            } else {
                text[position - 1]
            };
            self.symbol_bits[symbol as usize][row / 64] |= 1_u64 << (row % 64);
            if position % SA_SAMPLE_RATE == 0 {
                self.sample_bits[row / 64] |= 1_u64 << (row % 64);
                self.sample_values.push(position as u32);
            }
        }
        self.symbol_prefix = std::array::from_fn(|symbol| {
            let mut prefix = Vec::with_capacity(words + 1);
            prefix.push(0_u32);
            for &word in &self.symbol_bits[symbol] {
                prefix.push(prefix.last().copied().unwrap() + word.count_ones());
            }
            prefix
        });
        self.sample_prefix = Vec::with_capacity(words + 1);
        self.sample_prefix.push(0_u32);
        for &word in &self.sample_bits {
            self.sample_prefix
                .push(self.sample_prefix.last().copied().unwrap() + word.count_ones());
        }
        Ok(())
    }

    #[inline(always)]
    fn rank(&self, symbol: u8, end: usize) -> usize {
        let word = end / 64;
        let offset = end % 64;
        let mut count = self.symbol_prefix[symbol as usize][word] as usize;
        if offset > 0 {
            let mask = (1_u64 << offset) - 1;
            count += (self.symbol_bits[symbol as usize][word] & mask).count_ones() as usize;
        }
        count
    }

    #[inline(always)]
    fn extend(&self, interval: (usize, usize), symbol: u8) -> (usize, usize) {
        let base = self.cumulative[symbol as usize] as usize;
        (
            base + self.rank(symbol, interval.0),
            base + self.rank(symbol, interval.1),
        )
    }

    fn symbol_at(&self, row: usize) -> u8 {
        let mask = 1_u64 << (row % 64);
        for symbol in 0..ALPHABET_SIZE {
            if self.symbol_bits[symbol][row / 64] & mask != 0 {
                return symbol as u8;
            }
        }
        unreachable!("each BWT row has exactly one symbol")
    }

    fn sampled_sa(&self, row: usize) -> Option<usize> {
        let word = row / 64;
        let offset = row % 64;
        let bit = 1_u64 << offset;
        if self.sample_bits[word] & bit == 0 {
            return None;
        }
        let preceding_mask = bit - 1;
        let sample = self.sample_prefix[word] as usize
            + (self.sample_bits[word] & preceding_mask).count_ones() as usize;
        Some(self.sample_values[sample] as usize)
    }

    fn locate(&self, mut row: usize) -> usize {
        let mut steps = 0_usize;
        loop {
            if let Some(position) = self.sampled_sa(row) {
                return (position + steps) % self.length;
            }
            let symbol = self.symbol_at(row);
            row = self.cumulative[symbol as usize] as usize + self.rank(symbol, row);
            steps += 1;
            debug_assert!(steps < SA_SAMPLE_RATE);
        }
    }

    fn span_at(&self, position: usize) -> Option<ReferenceSpan> {
        let index = self
            .spans
            .partition_point(|span| span.text_end as usize <= position);
        self.spans.get(index).copied().filter(|span| {
            span.text_start as usize <= position && position < span.text_end as usize
        })
    }

    #[cfg(test)]
    fn interval(&self, pattern: &[u8]) -> (usize, usize) {
        let mut interval = (0_usize, self.length);
        for &base in pattern {
            interval = self.extend(interval, base_symbol(base));
            if interval.0 == interval.1 {
                break;
            }
        }
        interval
    }

    fn occurrence(&self, row: usize, length: usize) -> Option<(u32, usize, usize, u8)> {
        let text_start = self.locate(row);
        let span = self.span_at(text_start)?;
        if text_start + length > span.text_end as usize {
            return None;
        }
        let reference_end = span.reference_end as usize - (text_start - span.text_start as usize);
        let reference_start = reference_end.checked_sub(length)?;
        debug_assert!(reference_start >= span.reference_start as usize);
        Some((span.sequence, reference_start, reference_end, span.strand))
    }

    fn best_occurrence(
        &self,
        interval: (usize, usize),
        length: usize,
    ) -> Option<(u32, usize, usize)> {
        (interval.0..interval.1)
            .filter_map(|row| self.occurrence(row, length))
            .map(|(sequence, start, end, _)| (sequence, start, end))
            .min_by_key(|&(sequence, start, _)| (sequence, start))
    }

    pub fn collect(
        &self,
        read: &[u8],
        run_k: usize,
        occurrence_counts: &mut [u32],
        strand_masks: &mut [u8],
        profile: &mut MemQueryProfile,
    ) -> Option<MemSeed> {
        if self.length == 0 {
            return None;
        }
        profile.index_queries += 1;
        let mut best = None::<MemSeed>;
        for read_start in 0..occurrence_counts.len() {
            profile.run_windows += 1;
            let mut interval = (0_usize, self.length);
            let mut length = 0_usize;
            for &base in &read[read_start..] {
                let next = self.extend(interval, base_symbol(base));
                if next.0 == next.1 {
                    break;
                }
                interval = next;
                length += 1;
                if length == run_k {
                    let occurrences = interval.1 - interval.0;
                    profile.matching_windows += 1;
                    occurrence_counts[read_start] = occurrences.min(u32::MAX as usize) as u32;
                    if occurrences == 1 {
                        if let Some((_, _, _, strand)) = self.occurrence(interval.0, run_k) {
                            strand_masks[read_start] = strand;
                        }
                    }
                }
            }
            if length < run_k {
                continue;
            }
            profile.mem_starts += 1;
            profile.mem_bases += length as u64;
            let Some((sequence, reference_start, reference_end)) =
                self.best_occurrence(interval, length)
            else {
                continue;
            };
            let candidate = MemSeed {
                sequence,
                read_start,
                read_end: read_start + length,
                reference_start,
                reference_end,
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
        }
        best
    }

    pub fn symbols(&self) -> usize {
        self.length
    }
}

pub(crate) fn build_locus_indexes(
    references: &[OrientedReference],
    locus_count: usize,
    run_k: usize,
) -> Result<Vec<LocusMemIndex>, String> {
    let mut loci = (0..locus_count)
        .map(|_| LocusMemIndex::default())
        .collect::<Vec<_>>();
    let mut texts = (0..locus_count).map(|_| Vec::new()).collect::<Vec<_>>();
    for (sequence, reference) in references.iter().enumerate() {
        let sequence = u32::try_from(sequence)
            .map_err(|_| "reference sequence count exceeds u32".to_string())?;
        let locus_id = reference.locus as usize;
        let locus = &mut loci[locus_id];
        let text = &mut texts[locus_id];
        let mut cursor = 0_usize;
        while cursor < reference.bases.len() {
            while cursor < reference.bases.len() && code(reference.bases[cursor]).is_none() {
                cursor += 1;
            }
            let reference_start = cursor;
            while cursor < reference.bases.len() && code(reference.bases[cursor]).is_some() {
                cursor += 1;
            }
            let reference_end = cursor;
            if reference_end - reference_start < run_k {
                continue;
            }
            if !text.is_empty() {
                text.push(SEPARATOR);
            }
            let text_start = text.len();
            text.extend(
                reference.bases[reference_start..reference_end]
                    .iter()
                    .rev()
                    .map(|&base| base_symbol(base)),
            );
            locus.spans.push(ReferenceSpan {
                text_start: u32::try_from(text_start)
                    .map_err(|_| "locus text coordinate exceeds u32".to_string())?,
                text_end: u32::try_from(text.len())
                    .map_err(|_| "locus text coordinate exceeds u32".to_string())?,
                sequence,
                reference_start: u32::try_from(reference_start)
                    .map_err(|_| "reference coordinate exceeds u32".to_string())?,
                reference_end: u32::try_from(reference_end)
                    .map_err(|_| "reference coordinate exceeds u32".to_string())?,
                strand: reference.strand,
            });
        }
    }
    for (locus, text) in loci.iter_mut().zip(&mut texts) {
        if !text.is_empty() {
            text.push(SENTINEL);
            locus.build(text)?;
        }
    }
    Ok(loci)
}

#[inline(always)]
fn base_symbol(base: u8) -> u8 {
    code(base).map_or(SEPARATOR, |value| value + 1)
}

fn build_suffix_array(text: &[u8]) -> Result<Vec<usize>, String> {
    if text.len() > u32::MAX as usize {
        return Err("locus FM-index text length exceeds u32".to_string());
    }
    let length = text.len();
    let mut suffixes = (0..length).collect::<Vec<_>>();
    let mut ranks = text.iter().map(|&symbol| symbol as i64).collect::<Vec<_>>();
    let mut next_ranks = vec![0_i64; length];
    let mut width = 1_usize;
    while width < length {
        suffixes.sort_unstable_by_key(|&position| {
            (
                ranks[position],
                ranks.get(position + width).copied().unwrap_or(-1),
            )
        });
        if let Some(&first) = suffixes.first() {
            next_ranks[first] = 0;
        }
        for pair in suffixes.windows(2) {
            let previous = pair[0];
            let current = pair[1];
            let previous_key = (
                ranks[previous],
                ranks.get(previous + width).copied().unwrap_or(-1),
            );
            let current_key = (
                ranks[current],
                ranks.get(current + width).copied().unwrap_or(-1),
            );
            next_ranks[current] = next_ranks[previous] + i64::from(previous_key != current_key);
        }
        std::mem::swap(&mut ranks, &mut next_ranks);
        if suffixes
            .last()
            .is_some_and(|&position| ranks[position] as usize + 1 == length)
        {
            break;
        }
        width = width.saturating_mul(2);
    }
    Ok(suffixes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_index(texts: &[(&[u8], u8)]) -> LocusMemIndex {
        let mut index = LocusMemIndex::default();
        let mut text = Vec::new();
        for (sequence, &(bases, strand)) in texts.iter().enumerate() {
            if !text.is_empty() {
                text.push(SEPARATOR);
            }
            let start = text.len();
            text.extend(bases.iter().rev().map(|&base| base_symbol(base)));
            index.spans.push(ReferenceSpan {
                text_start: start as u32,
                text_end: text.len() as u32,
                sequence: sequence as u32,
                reference_start: 0,
                reference_end: bases.len() as u32,
                strand,
            });
        }
        text.push(SENTINEL);
        index.build(&text).unwrap();
        index
    }

    #[test]
    fn fm_intervals_count_repeats_and_locate_leftmost() {
        let index = text_index(&[(b"ACGTACGTACGA", 1)]);
        let interval = index.interval(b"ACG");
        assert_eq!(interval.1 - interval.0, 3);
        assert_eq!(index.best_occurrence(interval, 3), Some((0, 0, 3)));
    }

    #[test]
    fn longest_prefix_returns_leftmost_equal_occurrence() {
        let index = text_index(&[(b"TTACGTTTACGT", 1)]);
        let read = b"ACGTAA";
        let mut counts = vec![0; read.len() - 3 + 1];
        let mut strands = vec![0; counts.len()];
        let seed = index
            .collect(
                read,
                3,
                &mut counts,
                &mut strands,
                &mut MemQueryProfile::default(),
            )
            .unwrap();
        assert_eq!((seed.read_start, seed.read_end), (0, 4));
        assert_eq!((seed.reference_start, seed.reference_end), (2, 6));
    }

    #[test]
    fn sentinels_prevent_cross_sequence_matches() {
        let index = text_index(&[(b"AAAAC", 1), (b"GTTTT", 2)]);
        let interval = index.interval(b"ACG");
        assert_eq!(interval.0, interval.1);
    }
}
