use crate::model::Candidate;

const AUTO_TARGET_DEPTH: f64 = 80.0;
const AUTO_ACTIVATION_DEPTH: f64 = 160.0;
const AUTO_MIN_CANDIDATES: usize = 512;
const AUTO_MIN_COVERED_BINS: u32 = 48;
const AUTO_TERMINAL_LIMIT: usize = 768;
const AUTO_TERMINAL_PER_EXTENSION: usize = 4;
const AUTO_MIN_RETAIN_NUMERATOR: usize = 3;
const AUTO_MIN_RETAIN_DENOMINATOR: usize = 5;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LegacyDecision {
    pub minimum_exact: Option<usize>,
    pub thinning_interval: Option<usize>,
    pub total_bases: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AutoMode {
    PassThrough,
    Core,
    Rescue,
    LegacyFallback,
}

impl AutoMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::PassThrough => "pass-through",
            Self::Core => "core",
            Self::Rescue => "rescue",
            Self::LegacyFallback => "legacy-fallback",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AutoDecision {
    pub selected_ids: Vec<u32>,
    pub mode: AutoMode,
    pub eligible_fragments: usize,
    pub covered_bins: u32,
    pub target_core_fragments: usize,
    pub left_terminal_candidates: usize,
    pub right_terminal_candidates: usize,
    pub selected_left_terminal: usize,
    pub selected_right_terminal: usize,
}

pub fn choose_legacy(
    candidates: &[Candidate],
    effective_length: f64,
    initial_k: usize,
    min_depth: i64,
    max_depth: i64,
    max_size_mb: i64,
) -> LegacyDecision {
    if candidates.is_empty() || effective_length <= 0.0 {
        return LegacyDecision {
            minimum_exact: None,
            thinning_interval: None,
            total_bases: 0,
        };
    }
    let maximum = candidates
        .iter()
        .map(|v| v.max_exact as usize)
        .max()
        .unwrap_or(0)
        .max(initial_k + 8);
    let mut histogram = vec![0_u64; maximum + 2];
    for candidate in candidates {
        histogram[candidate.max_exact as usize] += candidate.fragment_bases as u64;
    }
    for i in (0..histogram.len() - 1).rev() {
        histogram[i] += histogram[i + 1];
    }
    let all_bases: u64 = candidates.iter().map(|v| v.fragment_bases as u64).sum();
    let mut total = all_bases;
    let mut coverage = total as f64 / effective_length;
    let mut too_deep = coverage > max_depth as f64;
    let mut too_large = total / 1_000_000 > max_size_mb as u64;
    if !too_deep && !too_large {
        return LegacyDecision {
            minimum_exact: None,
            thinning_interval: None,
            total_bases: total,
        };
    }
    let min_depth = (min_depth as f64).min(max_depth as f64 / 4.0);
    let mut current_k = initial_k;
    while current_k < 64 && (too_deep || too_large) {
        let last_k = current_k;
        let last_total = total;
        current_k +=
            if coverage > 8.0 * max_depth as f64 || total / 1_000_000 > (6 * max_size_mb) as u64 {
                6
            } else {
                2
            };
        total = histogram.get(current_k).copied().unwrap_or(0);
        coverage = total as f64 / effective_length;
        too_deep = coverage > max_depth as f64;
        too_large = total / 1_000_000 > max_size_mb as u64;
        if coverage < min_depth {
            current_k = last_k;
            total = last_total;
            too_large = total / 1_000_000 > max_size_mb as u64;
            break;
        }
    }
    if current_k == initial_k && !too_large {
        return LegacyDecision {
            minimum_exact: None,
            thinning_interval: None,
            total_bases: all_bases,
        };
    }
    let interval =
        too_large.then(|| ((total as f64 / 1e6 / max_size_mb as f64).trunc() as usize).max(2));
    LegacyDecision {
        minimum_exact: Some(current_k),
        thinning_interval: interval,
        total_bases: total,
    }
}

pub fn selected(
    candidate: &Candidate,
    decision: &LegacyDecision,
    passing_ordinal: &mut usize,
) -> bool {
    if decision
        .minimum_exact
        .is_some_and(|k| (candidate.max_exact as usize) < k)
    {
        return false;
    }
    *passing_ordinal += 1;
    decision
        .thinning_interval
        .is_none_or(|interval| passing_ordinal.is_multiple_of(interval))
}

fn legacy_ids(candidates: &[Candidate], decision: &LegacyDecision) -> Vec<u32> {
    let mut passing = 0_usize;
    candidates
        .iter()
        .filter(|candidate| selected(candidate, decision, &mut passing))
        .map(|candidate| candidate.fragment_id)
        .collect()
}

fn quality_key(candidate: &Candidate) -> (u16, std::cmp::Reverse<u16>, std::cmp::Reverse<u8>, u32) {
    (
        candidate.locus_count,
        std::cmp::Reverse(candidate.max_exact),
        std::cmp::Reverse(candidate.aligned_mates),
        candidate.fragment_id,
    )
}

fn terminal_quality_key(
    candidate: &Candidate,
    terminal_bit: u8,
) -> (
    std::cmp::Reverse<u16>,
    u16,
    std::cmp::Reverse<u16>,
    std::cmp::Reverse<u8>,
    u32,
) {
    let extension = if terminal_bit == 1 {
        candidate.left_extension
    } else {
        candidate.right_extension
    };
    (
        std::cmp::Reverse(extension),
        candidate.locus_count,
        std::cmp::Reverse(candidate.max_exact),
        std::cmp::Reverse(candidate.aligned_mates),
        candidate.fragment_id,
    )
}

fn midpoint_bin(mask: u64) -> Option<usize> {
    if mask == 0 {
        return None;
    }
    let first = mask.trailing_zeros() as usize;
    let last = 63 - mask.leading_zeros() as usize;
    Some((first + last) / 2)
}

/// Selects a bounded, reference-spanning subset only for saturated loci.
///
/// Low-depth loci pass through unchanged. A high-depth locus whose exact seeds
/// do not span most of the reference uses the legacy selector instead: this
/// prevents a local repeat pile-up from being mistaken for useful coverage.
pub fn choose_auto(
    candidates: &[Candidate],
    legacy: &LegacyDecision,
    effective_length: f64,
    reference_is_contig: bool,
) -> AutoDecision {
    let eligible: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, candidate)| {
            legacy
                .minimum_exact
                .is_none_or(|k| candidate.max_exact as usize >= k)
        })
        .map(|(index, _)| index)
        .collect();
    let total_bases: u64 = eligible
        .iter()
        .map(|&index| candidates[index].fragment_bases as u64)
        .sum();
    let covered = eligible
        .iter()
        .fold(0_u64, |mask, &index| mask | candidates[index].covered_bins);
    let covered_bins = covered.count_ones();
    let left_terminal_candidates = eligible
        .iter()
        .filter(|&&index| candidates[index].terminal_mask & 1 != 0)
        .count();
    let right_terminal_candidates = eligible
        .iter()
        .filter(|&&index| candidates[index].terminal_mask & 2 != 0)
        .count();
    let coverage = if effective_length > 0.0 {
        total_bases as f64 / effective_length
    } else {
        0.0
    };
    let pass_through = eligible.len() < AUTO_MIN_CANDIDATES
        || coverage <= AUTO_ACTIVATION_DEPTH
        || covered_bins < AUTO_MIN_COVERED_BINS;
    if pass_through {
        let (selected_ids, mode) = if legacy.thinning_interval.is_some() {
            (legacy_ids(candidates, legacy), AutoMode::LegacyFallback)
        } else {
            (
                eligible
                    .iter()
                    .map(|&index| candidates[index].fragment_id)
                    .collect(),
                AutoMode::PassThrough,
            )
        };
        let selected_left_terminal = candidates
            .iter()
            .filter(|candidate| {
                candidate.terminal_mask & 1 != 0
                    && selected_ids.binary_search(&candidate.fragment_id).is_ok()
            })
            .count();
        let selected_right_terminal = candidates
            .iter()
            .filter(|candidate| {
                candidate.terminal_mask & 2 != 0
                    && selected_ids.binary_search(&candidate.fragment_id).is_ok()
            })
            .count();
        return AutoDecision {
            selected_ids,
            mode,
            eligible_fragments: eligible.len(),
            covered_bins,
            target_core_fragments: eligible.len(),
            left_terminal_candidates,
            right_terminal_candidates,
            selected_left_terminal,
            selected_right_terminal,
        };
    }

    let average_bases = total_bases as f64 / eligible.len() as f64;
    let minimum_retained =
        (eligible.len() * AUTO_MIN_RETAIN_NUMERATOR).div_ceil(AUTO_MIN_RETAIN_DENOMINATOR);
    let target_core = ((AUTO_TARGET_DEPTH * effective_length / average_bases).ceil() as usize)
        .max(AUTO_MIN_CANDIDATES)
        .max(minimum_retained)
        .min(eligible.len());
    let mut ranked = eligible.clone();
    ranked.sort_unstable_by_key(|&index| quality_key(&candidates[index]));
    let mut chosen = vec![false; candidates.len()];

    // Bait edges are not biological contig ends, but reads crossing them carry
    // the UCE flanks. Protect them in both initial and rescue rounds.
    for terminal_bit in [1_u8, 2_u8] {
        let mut terminal_ranked: Vec<_> = eligible
            .iter()
            .copied()
            .filter(|&index| candidates[index].terminal_mask & terminal_bit != 0)
            .collect();
        terminal_ranked
            .sort_unstable_by_key(|&index| terminal_quality_key(&candidates[index], terminal_bit));
        let maximum_extension = terminal_ranked.first().map_or(0, |&index| {
            if terminal_bit == 1 {
                candidates[index].left_extension
            } else {
                candidates[index].right_extension
            }
        }) as usize;
        let mut per_extension = vec![0_u8; maximum_extension + 1];
        let mut terminal_count = 0_usize;
        for &index in &terminal_ranked {
            let extension = if terminal_bit == 1 {
                candidates[index].left_extension
            } else {
                candidates[index].right_extension
            } as usize;
            if per_extension[extension] as usize >= AUTO_TERMINAL_PER_EXTENSION {
                continue;
            }
            chosen[index] = true;
            per_extension[extension] += 1;
            terminal_count += 1;
            if terminal_count == AUTO_TERMINAL_LIMIT {
                break;
            }
        }
    }

    // First distribute core evidence across the reference, then fill by quality.
    let bin_quota = target_core.div_ceil(64).max(1);
    let mut bin_counts = [0_usize; 64];
    let mut core_count = 0_usize;
    for &index in &ranked {
        if core_count == target_core {
            break;
        }
        let Some(bin) = midpoint_bin(candidates[index].covered_bins) else {
            continue;
        };
        if bin_counts[bin] < bin_quota {
            chosen[index] = true;
            bin_counts[bin] += 1;
            core_count += 1;
        }
    }
    for &index in &ranked {
        if core_count == target_core {
            break;
        }
        if !chosen[index] {
            chosen[index] = true;
            core_count += 1;
        }
    }
    let selected_ids = candidates
        .iter()
        .enumerate()
        .filter(|(index, _)| chosen[*index])
        .map(|(_, candidate)| candidate.fragment_id)
        .collect();
    let selected_left_terminal = candidates
        .iter()
        .enumerate()
        .filter(|(index, candidate)| chosen[*index] && candidate.terminal_mask & 1 != 0)
        .count();
    let selected_right_terminal = candidates
        .iter()
        .enumerate()
        .filter(|(index, candidate)| chosen[*index] && candidate.terminal_mask & 2 != 0)
        .count();
    AutoDecision {
        selected_ids,
        mode: if reference_is_contig {
            AutoMode::Rescue
        } else {
            AutoMode::Core
        },
        eligible_fragments: eligible.len(),
        covered_bins,
        target_core_fragments: target_core,
        left_terminal_candidates,
        right_terminal_candidates,
        selected_left_terminal,
        selected_right_terminal,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(ordinal: u64, bases: u32, maximum: u16) -> Candidate {
        Candidate {
            fragment_id: ordinal as u32,
            fragment_bases: bases,
            max_exact: maximum,
            covered_bins: u64::MAX,
            terminal_mask: 0,
            left_extension: 0,
            right_extension: 0,
            aligned_mates: 1,
            locus_count: 1,
        }
    }

    #[test]
    fn low_depth_keeps_every_run_accepted_fragment() {
        let candidates = vec![candidate(0, 300, 19), candidate(1, 300, 40)];
        let decision = choose_legacy(&candidates, 1_000.0, 31, 50, 768, 6);
        assert_eq!(decision.minimum_exact, None);
        assert_eq!(decision.thinning_interval, None);
    }

    #[test]
    fn high_depth_raises_exact_threshold() {
        let mut candidates = Vec::new();
        for i in 0..10_000 {
            candidates.push(candidate(i, 300, if i % 2 == 0 { 33 } else { 45 }));
        }
        let decision = choose_legacy(&candidates, 1_000.0, 31, 50, 100, 100);
        assert!(decision.minimum_exact.is_some_and(|value| value > 31));
    }

    #[test]
    fn auto_passes_low_depth_loci_through() {
        let candidates = vec![candidate(0, 300, 31), candidate(1, 300, 40)];
        let legacy = choose_legacy(&candidates, 1_000.0, 31, 50, 768, 6);
        let decision = choose_auto(&candidates, &legacy, 1_000.0, false);
        assert_eq!(decision.mode, AutoMode::PassThrough);
        assert_eq!(decision.selected_ids, vec![0, 1]);
    }

    #[test]
    fn auto_compresses_only_saturated_reference_spanning_loci() {
        let candidates: Vec<_> = (0..2_000)
            .map(|i| {
                let mut value = candidate(i, 300, 50);
                value.covered_bins = 1_u64 << (i % 64);
                value
            })
            .collect();
        let legacy = LegacyDecision {
            minimum_exact: None,
            thinning_interval: None,
            total_bases: 600_000,
        };
        let decision = choose_auto(&candidates, &legacy, 1_000.0, false);
        assert_eq!(decision.mode, AutoMode::Core);
        assert_eq!(decision.target_core_fragments, 1_200);
        assert_eq!(decision.selected_ids.len(), 1_200);
    }

    #[test]
    fn auto_rescue_keeps_separate_left_and_right_terminal_evidence() {
        let mut candidates: Vec<_> = (0..2_000)
            .map(|i| {
                let mut value = candidate(i, 300, 50);
                value.covered_bins = 1_u64 << (i % 64);
                value
            })
            .collect();
        candidates[1_900].terminal_mask = 1;
        candidates[1_900].left_extension = 40;
        candidates[1_901].terminal_mask = 2;
        candidates[1_901].right_extension = 40;
        let legacy = LegacyDecision {
            minimum_exact: None,
            thinning_interval: None,
            total_bases: 600_000,
        };
        let decision = choose_auto(&candidates, &legacy, 1_000.0, true);
        assert_eq!(decision.mode, AutoMode::Rescue);
        assert!(decision.selected_ids.contains(&1_900));
        assert!(decision.selected_ids.contains(&1_901));
    }

    #[test]
    fn auto_uses_legacy_fallback_for_deep_but_localized_pileups() {
        let mut candidates: Vec<_> = (0..1_000).map(|i| candidate(i, 300, 50)).collect();
        for candidate in &mut candidates {
            candidate.covered_bins = 1;
        }
        let legacy = LegacyDecision {
            minimum_exact: Some(40),
            thinning_interval: Some(5),
            total_bases: 300_000,
        };
        let decision = choose_auto(&candidates, &legacy, 1_000.0, false);
        assert_eq!(decision.mode, AutoMode::LegacyFallback);
        assert_eq!(decision.selected_ids.len(), 200);
    }
}
