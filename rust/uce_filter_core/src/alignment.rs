use crate::index::{ExactSeed, UceIndex};
use crate::model::LocusId;

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -1;
const STOP: u8 = 3;

#[derive(Clone, Copy, Debug)]
pub struct AlignmentEvidence {
    pub sequence: u32,
    pub strand: u8,
    pub score: i32,
    pub exact_seed_length: u16,
    pub query_start: u16,
    pub query_end: u16,
    pub reference_start: u32,
    pub reference_end: u32,
    pub matches: u16,
    pub mismatches: u16,
    pub gap_bases: u16,
}

impl AlignmentEvidence {
    pub fn aligned_columns(self) -> usize {
        self.matches as usize + self.mismatches as usize + self.gap_bases as usize
    }

    pub fn reference_overlap(self) -> usize {
        self.reference_end as usize - self.reference_start as usize
    }

    pub fn identity(self) -> f64 {
        let columns = self.aligned_columns();
        if columns == 0 {
            0.0
        } else {
            self.matches as f64 / columns as f64
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct TerminalEvidence {
    pub effective_window: u32,
    pub original_reference_start: u32,
    pub original_reference_end: u32,
    pub left_overhang: u16,
    pub right_overhang: u16,
    pub near_left_terminal: bool,
    pub near_right_terminal: bool,
    pub extends_left_terminal: bool,
    pub extends_right_terminal: bool,
}

pub fn terminal_evidence(
    alignment: AlignmentEvidence,
    query_length: usize,
    reference_length: usize,
    terminal_window: usize,
) -> TerminalEvidence {
    let query_left = alignment.query_start as usize;
    let query_right = query_length.saturating_sub(alignment.query_end as usize);
    let (reference_start, reference_end, left_overhang, right_overhang) = if alignment.strand == 1 {
        (
            alignment.reference_start as usize,
            alignment.reference_end as usize,
            query_left,
            query_right,
        )
    } else {
        (
            reference_length.saturating_sub(alignment.reference_end as usize),
            reference_length.saturating_sub(alignment.reference_start as usize),
            query_right,
            query_left,
        )
    };
    let effective_window = terminal_window.min((reference_length / 5).max(1));
    let near_left_terminal = reference_start <= effective_window;
    let near_right_terminal = reference_length.saturating_sub(reference_end) <= effective_window;
    TerminalEvidence {
        effective_window: effective_window.min(u32::MAX as usize) as u32,
        original_reference_start: reference_start.min(u32::MAX as usize) as u32,
        original_reference_end: reference_end.min(u32::MAX as usize) as u32,
        left_overhang: left_overhang.min(u16::MAX as usize) as u16,
        right_overhang: right_overhang.min(u16::MAX as usize) as u16,
        near_left_terminal,
        near_right_terminal,
        extends_left_terminal: near_left_terminal && left_overhang > 0,
        extends_right_terminal: near_right_terminal && right_overhang > 0,
    }
}

fn best_transition(candidates: &[(i32, u8)]) -> (i32, u8) {
    let mut best = (0, STOP);
    for &(score, state) in candidates {
        if score > best.0 || (score == best.0 && score > 0 && state < best.1) {
            best = (score, state);
        }
    }
    best
}

fn seeded_local_affine(
    query: &[u8],
    reference: &[u8],
    seed: ExactSeed,
    strand: u8,
    band: usize,
) -> Option<AlignmentEvidence> {
    if query.is_empty() || reference.is_empty() || seed.is_empty() {
        return None;
    }
    let expected_start = (seed.reference_start as usize).saturating_sub(seed.read_start as usize);
    let slice_start = expected_start.saturating_sub(band);
    let slice_end = (expected_start + query.len() + band).min(reference.len());
    if slice_start >= slice_end {
        return None;
    }
    let target = &reference[slice_start..slice_end];
    let rows = query.len() + 1;
    let cols = target.len() + 1;
    let cells = rows.checked_mul(cols)?;
    let mut scores = [vec![0_i32; cells], vec![0_i32; cells], vec![0_i32; cells]];
    let mut trace = [vec![STOP; cells], vec![STOP; cells], vec![STOP; cells]];
    let mut best = (0_i32, 0_usize, 0_usize, 0_u8);

    for i in 1..rows {
        for j in 1..cols {
            let idx = i * cols + j;
            let diagonal = (i - 1) * cols + j - 1;
            let above = (i - 1) * cols + j;
            let left = i * cols + j - 1;
            let substitution = if query[i - 1] == target[j - 1] {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            let matched = best_transition(&[
                (scores[0][diagonal] + substitution, 0),
                (scores[1][diagonal] + substitution, 1),
                (scores[2][diagonal] + substitution, 2),
            ]);
            scores[0][idx] = matched.0;
            trace[0][idx] = matched.1;

            let inserted = best_transition(&[
                (scores[0][above] + GAP_OPEN, 0),
                (scores[1][above] + GAP_EXTEND, 1),
                (scores[2][above] + GAP_OPEN, 2),
            ]);
            scores[1][idx] = inserted.0;
            trace[1][idx] = inserted.1;

            let deleted = best_transition(&[
                (scores[0][left] + GAP_OPEN, 0),
                (scores[1][left] + GAP_OPEN, 1),
                (scores[2][left] + GAP_EXTEND, 2),
            ]);
            scores[2][idx] = deleted.0;
            trace[2][idx] = deleted.1;

            for state in 0..3_u8 {
                let score = scores[state as usize][idx];
                if score > best.0
                    || (score == best.0 && score > 0 && (i, j, state) < (best.1, best.2, best.3))
                {
                    best = (score, i, j, state);
                }
            }
        }
    }
    if best.0 <= 0 {
        return None;
    }

    let (score, query_end, target_end, mut state) = best;
    let (mut i, mut j) = (query_end, target_end);
    let (mut matches, mut mismatches, mut gaps) = (0_usize, 0_usize, 0_usize);
    loop {
        let idx = i * cols + j;
        if state == STOP || scores[state as usize][idx] == 0 {
            break;
        }
        let previous = trace[state as usize][idx];
        match state {
            0 => {
                if query[i - 1] == target[j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                gaps += 1;
                i -= 1;
            }
            2 => {
                gaps += 1;
                j -= 1;
            }
            _ => break,
        }
        state = previous;
    }
    Some(AlignmentEvidence {
        sequence: seed.sequence,
        strand,
        score,
        exact_seed_length: seed.len().min(u16::MAX as usize) as u16,
        query_start: i.min(u16::MAX as usize) as u16,
        query_end: query_end.min(u16::MAX as usize) as u16,
        reference_start: (slice_start + j).min(u32::MAX as usize) as u32,
        reference_end: (slice_start + target_end).min(u32::MAX as usize) as u32,
        matches: matches.min(u16::MAX as usize) as u16,
        mismatches: mismatches.min(u16::MAX as usize) as u16,
        gap_bases: gaps.min(u16::MAX as usize) as u16,
    })
}

pub fn align_read(
    index: &UceIndex,
    query: &[u8],
    locus: LocusId,
    band: usize,
) -> Option<AlignmentEvidence> {
    let seed = index.best_exact(query, locus)?;
    let reference = index.references.get(seed.sequence as usize)?;
    seeded_local_affine(query, &reference.bases, seed, reference.strand, band)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_alignment_reports_identity_and_terminal_overhang() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let query = b"TTTTACGTTGCAACGATTCGGTACC";
        let seed = ExactSeed {
            sequence: 0,
            read_start: 4,
            read_end: 25,
            reference_start: 0,
            reference_end: 21,
        };
        let alignment = seeded_local_affine(query, reference, seed, 1, 16).unwrap();
        assert_eq!(alignment.matches, 21);
        assert_eq!(alignment.mismatches, 0);
        assert_eq!(alignment.identity(), 1.0);
        let terminal = terminal_evidence(alignment, query.len(), reference.len(), 4);
        assert_eq!(terminal.left_overhang, 4);
        assert!(terminal.extends_left_terminal);
        assert!(!terminal.extends_right_terminal);
    }

    #[test]
    fn reverse_alignment_swaps_terminal_sides() {
        let alignment = AlignmentEvidence {
            sequence: 0,
            strand: 2,
            score: 20,
            exact_seed_length: 10,
            query_start: 3,
            query_end: 13,
            reference_start: 0,
            reference_end: 10,
            matches: 10,
            mismatches: 0,
            gap_bases: 0,
        };
        let terminal = terminal_evidence(alignment, 13, 100, 5);
        assert!(terminal.near_right_terminal);
        assert_eq!(terminal.right_overhang, 3);
        assert!(terminal.extends_right_terminal);
    }
}
