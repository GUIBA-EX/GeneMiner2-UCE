use crate::model::{
    Args, AssemblyMode, ContigRecord, KmerInfo, Node, PathContig, PathStrategy, ReadSupport,
    RefKmer, SideContig,
};
use crate::seq::{
    bits_base, decode_kmer, for_each_kmer, kmer_mask, median, quartiles, reverse_complement,
    reverse_complement_kmer, valid_runs,
};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

pub fn build_reads_dictionary(sequences: &[Vec<u8>]) -> (HashMap<Vec<u8>, u64>, usize) {
    let Some(minimum) = sequences
        .iter()
        .filter(|seq| !seq.is_empty())
        .map(Vec::len)
        .min()
    else {
        return (HashMap::new(), 0);
    };
    let slice_len = ((minimum as f64) * 0.9) as usize;
    let slice_len = slice_len.max(1);
    let mut reads = HashMap::new();

    for sequence in sequences {
        if sequence.len() < slice_len {
            continue;
        }
        let start = (sequence.len() - slice_len) / 2;
        let forward = sequence[start..start + slice_len].to_vec();
        let reverse_sequence = reverse_complement(sequence);
        let reverse = reverse_sequence[start..start + slice_len].to_vec();
        *reads.entry(forward.clone()).or_default() += 1;
        if reverse != forward {
            *reads.entry(reverse).or_default() += 1;
        }
    }
    (reads, slice_len)
}

pub fn build_assemble_dictionary(
    sequences: &[Vec<u8>],
    k: usize,
    reference: &HashMap<u128, RefKmer>,
) -> HashMap<u128, KmerInfo> {
    let mut graph: HashMap<u128, KmerInfo> = HashMap::new();
    for sequence in sequences {
        let mut physical_read_kmers = HashSet::new();
        for (_, run) in valid_runs(sequence) {
            for_each_kmer(run, k, |_, forward, reverse| {
                physical_read_kmers.insert(forward);
                physical_read_kmers.insert(reverse);
            });
        }
        for kmer in physical_read_kmers {
            if let Some(info) = graph.get_mut(&kmer) {
                info.depth += 1;
                continue;
            }
            if let Some(reference_info) = reference.get(&kmer) {
                let position = if reference_info.is_reverse {
                    1000 - reference_info.position
                } else {
                    reference_info.position
                };
                graph.insert(
                    kmer,
                    KmerInfo {
                        depth: 1,
                        position,
                        is_reverse: reference_info.is_reverse,
                        reference_weight: reference_info.depth as i64,
                    },
                );
            } else {
                graph.insert(
                    kmer,
                    KmerInfo {
                        depth: 1,
                        position: 1023,
                        is_reverse: true,
                        reference_weight: 0,
                    },
                );
            }
        }
    }
    graph
}

pub fn filter_and_weight_graph(
    graph: &mut HashMap<u128, KmerInfo>,
    error_limit: u32,
    reference_count: usize,
) {
    if error_limit > 0 {
        graph.retain(|_, value| value.depth > error_limit as i64 || value.reference_weight > 0);
    }
    if graph.is_empty() {
        return;
    }

    let depths: Vec<i64> = graph.values().map(|value| value.depth).collect();
    let (q1, _, q3, _) = quartiles(&depths);
    let depth_upper = (((q3 - q1) * 1.5) + q3) as i64;
    for value in graph.values_mut() {
        if value.reference_weight != 0 {
            value.reference_weight = if value.depth > error_limit as i64 {
                (reference_count as i64 * depth_upper
                    / ((value.reference_weight - reference_count as i64).abs() + 1))
                    + 1
            } else {
                1
            };
        }
        value.depth = value.depth.min(depth_upper);
    }
}

fn outgoing(
    graph: &HashMap<u128, KmerInfo>,
    current: u128,
    k: usize,
    blocked: &HashSet<u128>,
    discarded: Option<&HashSet<u128>>,
    reference_tie_break: bool,
) -> Vec<Node> {
    let suffix_mask = kmer_mask(k - 1);
    let prefix = (current & suffix_mask) << 2;
    let mut nodes = Vec::with_capacity(4);
    for base in 0..4_u128 {
        let candidate = prefix | base;
        if blocked.contains(&candidate) || discarded.is_some_and(|items| items.contains(&candidate))
        {
            continue;
        }
        if let Some(value) = graph.get(&candidate) {
            nodes.push(Node {
                kmer: candidate,
                position: value.position,
                weight: value.depth + value.reference_weight,
            });
        }
    }
    nodes.sort_by(|left, right| {
        right.weight.cmp(&left.weight).then_with(|| {
            if reference_tie_break {
                (right.position > 0).cmp(&(left.position > 0))
            } else {
                Ordering::Equal
            }
        })
    });
    nodes
}

pub fn walk_backbone(
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    iteration: usize,
    lookahead: usize,
) -> (Vec<PathContig>, HashSet<u128>, Vec<i32>, i64) {
    let lookahead = lookahead.max(1);
    let mut path = vec![seed];
    let mut visited = HashSet::from([seed]);
    let mut discarded = HashSet::new();
    let mut positions = Vec::new();
    let mut result = PathContig::default();

    for _ in 0..iteration {
        let nodes = outgoing(
            graph,
            *path.last().expect("seed"),
            k,
            &visited,
            Some(&discarded),
            true,
        );
        if nodes.is_empty() {
            break;
        }

        let chosen = if nodes.len() == 1 {
            nodes[0]
        } else {
            let mut best_trace: Option<Vec<Node>> = None;
            for first in &nodes {
                let mut trace = Vec::new();
                let mut trace_seen = visited.clone();
                let mut node = *first;
                for _ in 0..lookahead {
                    if !trace_seen.insert(node.kmer) {
                        break;
                    }
                    trace.push(node);
                    let following =
                        outgoing(graph, node.kmer, k, &trace_seen, Some(&discarded), true);
                    let Some(next) = following.first() else {
                        break;
                    };
                    node = *next;
                }

                let replace = best_trace.as_ref().is_none_or(|best| {
                    let trace_sum: i64 = trace.iter().map(|node| node.weight).sum();
                    let best_sum: i64 = best.iter().map(|node| node.weight).sum();
                    (trace.len(), trace_sum, trace[0].weight)
                        > (best.len(), best_sum, best[0].weight)
                });
                if replace {
                    best_trace = Some(trace);
                }
            }
            let winner = best_trace.expect("outgoing branch has a trace")[0];
            for node in &nodes {
                if node.kmer != winner.kmer {
                    discarded.insert(node.kmer);
                }
            }
            winner
        };

        path.push(chosen.kmer);
        visited.insert(chosen.kmer);
        positions.push(chosen.position);
        result.weights.push(chosen.weight);
        result.bases.push((chosen.kmer & 3) as u8);
    }

    let weight = result.weights.iter().sum();
    (vec![result], visited, positions, weight)
}

pub fn walk_search(
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    mut iteration: usize,
) -> (Vec<PathContig>, HashSet<u128>, Vec<i32>, i64) {
    let mut path = vec![seed];
    let mut path_set = HashSet::from([seed]);
    let mut all_visited = HashSet::from([seed]);
    let mut positions = Vec::new();
    let mut current = PathContig::default();
    let mut contigs = Vec::new();
    let mut stack: Vec<(Vec<Node>, usize)> = Vec::new();
    let mut node_distance = 0_usize;
    let mut best_weight = 0_i64;

    while iteration > 0 {
        let mut nodes = outgoing(
            graph,
            *path.last().expect("seed"),
            k,
            &path_set,
            None,
            false,
        );
        if nodes.is_empty() {
            iteration -= 1;
            let weight: i64 = current.weights.iter().sum();
            best_weight = best_weight.max(weight);
            contigs.push(current.clone());

            for _ in 0..node_distance {
                if let Some(removed) = path.pop() {
                    path_set.remove(&removed);
                }
                current.weights.pop();
                current.bases.pop();
            }
            let Some((alternatives, distance)) = stack.pop() else {
                break;
            };
            nodes = alternatives;
            node_distance = distance;
        }

        if nodes.len() >= 2 {
            stack.push((nodes[1..].to_vec(), node_distance));
            node_distance = 0;
        }
        let chosen = nodes[0];
        path.push(chosen.kmer);
        path_set.insert(chosen.kmer);
        all_visited.insert(chosen.kmer);
        positions.push(chosen.position);
        current.weights.push(chosen.weight);
        current.bases.push((chosen.kmer & 3) as u8);
        node_distance += 1;
    }

    (contigs, all_visited, positions, best_weight)
}

fn locate_read_slices<'a>(
    sequence: &'a [u8],
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
) -> HashMap<&'a [u8], Option<usize>> {
    let mut matches = HashMap::new();
    if slice_len == 0 || sequence.len() < slice_len {
        return matches;
    }
    for position in 0..=sequence.len() - slice_len {
        let slice = &sequence[position..position + slice_len];
        if !reads.contains_key(slice) {
            continue;
        }
        if let Some(existing) = matches.get_mut(slice) {
            *existing = None;
        } else {
            matches.insert(slice, Some(position));
        }
    }
    matches
}

fn process_sides(
    mut contigs: Vec<PathContig>,
    max_weight: i64,
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
    soft_boundary: usize,
    mode: AssemblyMode,
) -> Vec<SideContig> {
    if mode == AssemblyMode::Reference {
        for contig in &mut contigs {
            if contig.weights.len() > 2 {
                let (q1, _, _, _) = quartiles(&contig.weights);
                let cut_position = contig
                    .weights
                    .iter()
                    .rposition(|weight| *weight as f64 >= q1);
                if let Some(cut) = cut_position {
                    let retained = cut.saturating_add(soft_boundary).saturating_add(1);
                    if retained < contig.weights.len() {
                        contig.weights.truncate(retained);
                        contig.bases.truncate(retained);
                    }
                }
            }
        }
    }

    let minimum = max_weight >> if mode == AssemblyMode::Uce { 2 } else { 1 };
    let mut processed = Vec::new();
    for contig in contigs {
        let weight: i64 = contig.weights.iter().sum();
        if weight <= minimum {
            continue;
        }
        let sequence: Vec<u8> = contig.bases.iter().map(|base| bits_base(*base)).collect();
        let matches = locate_read_slices(&sequence, slice_len, reads);
        let read_count = matches.keys().filter_map(|slice| reads.get(*slice)).sum();
        processed.push(SideContig {
            sequence,
            weight,
            read_count,
        });
    }

    if mode == AssemblyMode::Uce {
        processed.sort_by(|left, right| {
            right
                .sequence
                .len()
                .cmp(&left.sequence.len())
                .then_with(|| right.read_count.cmp(&left.read_count))
                .then_with(|| right.weight.cmp(&left.weight))
        });
    } else {
        processed.sort_by(|left, right| right.read_count.cmp(&left.read_count));
    }
    processed
}

pub fn calculate_read_support(
    sequence: &[u8],
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
) -> ReadSupport {
    let contig_len = sequence.len();
    let matches = locate_read_slices(sequence, slice_len, reads);
    let total_read_count = matches.keys().filter_map(|slice| reads.get(*slice)).sum();
    let mut unique_read_count = 0_u64;
    let mut multi_mapping_read_count = 0_u64;
    let mut intervals = Vec::new();

    for (slice, position) in matches {
        let count = reads.get(slice).copied().unwrap_or(0);
        if let Some(start) = position {
            unique_read_count += count;
            intervals.push((start, (start + slice_len).min(contig_len)));
        } else {
            multi_mapping_read_count += count;
        }
    }

    if intervals.is_empty() {
        return ReadSupport {
            total_read_count,
            multi_mapping_read_count,
            max_gap: contig_len,
            left_coord: contig_len,
            ..ReadSupport::default()
        };
    }

    intervals.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in intervals {
        if merged
            .last()
            .is_none_or(|(_, previous_end)| start > *previous_end)
        {
            merged.push((start, end));
        } else if let Some((_, previous_end)) = merged.last_mut() {
            *previous_end = (*previous_end).max(end);
        }
    }

    let left_coord = merged[0].0;
    let right_coord = merged.last().expect("nonempty").1;
    let supported_extent = right_coord - left_coord;
    let supported_bases = merged.iter().map(|(start, end)| end - start).sum();
    let breadth = if contig_len > 0 {
        supported_bases as f64 / contig_len as f64
    } else {
        0.0
    };
    let mut max_gap = left_coord.max(contig_len - right_coord);
    for pair in merged.windows(2) {
        max_gap = max_gap.max(pair[1].0 - pair[0].1);
    }
    let left_extension = left_coord;
    let right_extension = contig_len - right_coord;
    let flank_balance = if supported_extent == 0 {
        0.0
    } else if left_extension == 0 && right_extension == 0 {
        1.0
    } else {
        left_extension.min(right_extension) as f64 / left_extension.max(right_extension) as f64
    };

    ReadSupport {
        total_read_count,
        unique_read_count,
        multi_mapping_read_count,
        supported_extent,
        supported_bases,
        breadth,
        max_gap,
        flank_balance,
        left_coord,
        right_coord,
    }
}

fn depth_stats(sequence: &[u8], k: usize, graph: &HashMap<u128, KmerInfo>) -> (f64, f64, f64) {
    if sequence.len() < k {
        return (0.0, 0.0, 0.0);
    }
    let mut counts = Vec::with_capacity(sequence.len() - k + 1);
    for_each_kmer(sequence, k, |_, forward, _| {
        counts.push(graph.get(&forward).map_or(0, |value| value.depth));
    });
    let median_depth = median(&counts);
    if counts.is_empty() || median_depth <= 0.0 {
        return (median_depth, 0.0, 0.0);
    }
    let mean = counts.iter().sum::<i64>() as f64 / counts.len() as f64;
    let variance = counts
        .iter()
        .map(|count| (*count as f64 - mean).powi(2))
        .sum::<f64>()
        / counts.len() as f64;
    let cv = if mean > 0.0 {
        variance.sqrt() / mean
    } else {
        0.0
    };
    let maximum = counts.iter().copied().max().unwrap_or(0) as f64;
    (median_depth, cv, maximum / median_depth)
}

fn rejection_reasons(
    args: &Args,
    length: usize,
    unique_read_count: u64,
    supported_bases: usize,
    unique_read_density: f64,
    depth_cv: f64,
    max_depth_ratio: f64,
) -> Vec<&'static str> {
    let mut reasons = Vec::new();
    if unique_read_count == 0 {
        reasons.push("no_unique_read_support");
    }
    if supported_bases == 0 {
        reasons.push("no_positional_support");
    }
    if args.max_contig_length > 0 && length > args.max_contig_length {
        reasons.push("contig_too_long");
    }
    if length >= args.density_check_min_length && unique_read_density < args.min_read_density {
        reasons.push("low_unique_read_density");
    }
    if args.max_depth_cv > 0.0 && depth_cv > args.max_depth_cv {
        reasons.push("high_depth_cv");
    }
    if args.max_depth_ratio > 0.0 && max_depth_ratio > args.max_depth_ratio {
        reasons.push("repeat_depth_peak");
    }
    reasons
}

pub fn assemble_seed(
    args: &Args,
    reads: &HashMap<Vec<u8>, u64>,
    slice_len: usize,
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    soft_boundary: usize,
) -> (Vec<ContigRecord>, HashSet<u128>, i32) {
    let reverse_seed = reverse_complement_kmer(seed, k);
    let (right_paths, right_kmers, right_positions, right_weight) = if args.assembly_mode
        == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        walk_backbone(graph, seed, k, args.iteration, args.backbone_lookahead)
    } else {
        walk_search(graph, seed, k, args.iteration)
    };
    let (left_paths, left_kmers, left_positions, left_weight) = if args.assembly_mode
        == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        walk_backbone(
            graph,
            reverse_seed,
            k,
            args.iteration,
            args.backbone_lookahead,
        )
    } else {
        walk_search(graph, reverse_seed, k, args.iteration)
    };

    let mut positions: Vec<i64> = right_positions
        .into_iter()
        .chain(left_positions)
        .filter(|position| *position > 0 && *position < 1000)
        .map(i64::from)
        .collect();
    let contig_position = if positions.len() > 1 {
        median(&positions) as i32
    } else {
        -1
    };
    positions.clear();

    let mut right = process_sides(
        right_paths,
        right_weight,
        slice_len,
        reads,
        soft_boundary,
        args.assembly_mode,
    );
    let mut left = process_sides(
        left_paths,
        left_weight,
        slice_len,
        reads,
        soft_boundary,
        args.assembly_mode,
    );
    if right.is_empty() {
        right.push(SideContig::default());
    }
    if left.is_empty() {
        left.push(SideContig::default());
    }
    let candidate_limit = if args.assembly_mode == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        1
    } else if args.assembly_mode == AssemblyMode::Uce {
        args.side_candidates
    } else {
        3
    };

    let seed_sequence = decode_kmer(seed, k);
    let mut candidates = Vec::new();
    for left_side in left.iter().take(candidate_limit) {
        for right_side in right.iter().take(candidate_limit) {
            let mut sequence = reverse_complement(&left_side.sequence);
            sequence.extend_from_slice(&seed_sequence);
            sequence.extend_from_slice(&right_side.sequence);
            let weight = left_side.weight + right_side.weight;
            let mut support = calculate_read_support(&sequence, slice_len, reads);
            let mut length = sequence.len();
            let mut depth = if args.assembly_mode == AssemblyMode::Uce {
                depth_stats(&sequence, k, graph)
            } else {
                (0.0, 0.0, 0.0)
            };

            if args.min_coverage > 0.0 {
                let positional_count = if args.assembly_mode == AssemblyMode::Uce {
                    support.unique_read_count
                } else {
                    support.total_read_count
                };
                let positional_span = if args.assembly_mode == AssemblyMode::Uce {
                    support.supported_bases
                } else {
                    support.supported_extent
                };
                let coverage_depth = positional_count as f64 * slice_len as f64 / 0.9;
                if positional_span == 0
                    || coverage_depth / (positional_span as f64) < args.min_coverage
                {
                    continue;
                }
                if coverage_depth / (length as f64) < args.min_coverage {
                    sequence = sequence[support.left_coord..support.right_coord].to_vec();
                    length = sequence.len();
                    support = calculate_read_support(&sequence, slice_len, reads);
                    if args.assembly_mode == AssemblyMode::Uce {
                        depth = depth_stats(&sequence, k, graph);
                    }
                }
            }

            let read_density = if length > 0 {
                support.total_read_count as f64 / length as f64
            } else {
                0.0
            };
            let support_fraction = if length > 0 {
                support.supported_extent as f64 / length as f64
            } else {
                0.0
            };
            let unique_density = if length > 0 {
                support.unique_read_count as f64 / length as f64
            } else {
                0.0
            };
            let reasons = if args.assembly_mode == AssemblyMode::Uce {
                rejection_reasons(
                    args,
                    length,
                    support.unique_read_count,
                    support.supported_bases,
                    unique_density,
                    depth.1,
                    depth.2,
                )
            } else {
                Vec::new()
            };

            candidates.push(ContigRecord {
                sequence,
                position: contig_position,
                weight,
                read_count: support.total_read_count,
                supported_span: support.supported_extent,
                flank_balance: support.flank_balance,
                read_density,
                support_fraction,
                kmer_median_depth: depth.0,
                kmer_depth_cv: depth.1,
                kmer_max_depth_ratio: depth.2,
                unique_read_count: support.unique_read_count,
                multi_mapping_read_count: support.multi_mapping_read_count,
                supported_bases: support.supported_bases,
                support_breadth: support.breadth,
                max_support_gap: support.max_gap,
                accepted: reasons.is_empty(),
                rejection_reason: reasons.join(";"),
                ..ContigRecord::default()
            });
        }
    }

    let all_kmers = right_kmers.union(&left_kmers).copied().collect();
    (candidates, all_kmers, contig_position)
}

pub fn compare_contigs(left: &ContigRecord, right: &ContigRecord, mode: AssemblyMode) -> Ordering {
    if mode == AssemblyMode::Reference {
        return left
            .read_count
            .cmp(&right.read_count)
            .then_with(|| left.weight.cmp(&right.weight));
    }

    fn score(contig: &ContigRecord) -> [f64; 9] {
        let length = contig.sequence.len();
        let unique_density = if length > 0 {
            contig.unique_read_count as f64 / length as f64
        } else {
            0.0
        };
        let density_factor = (unique_density / 0.01).min(1.0);
        let continuity_factor = 1.0 / (1.0 + contig.kmer_depth_cv);
        let repeat_factor = (10.0 / contig.kmer_max_depth_ratio.max(1.0)).min(1.0);
        let effective =
            contig.supported_bases as f64 * density_factor * continuity_factor * repeat_factor;
        let gap_fraction = if length > 0 {
            contig.max_support_gap as f64 / length as f64
        } else {
            1.0
        };
        [
            effective,
            contig.support_breadth,
            unique_density,
            -gap_fraction,
            length as f64,
            contig.unique_read_count as f64,
            contig.read_count as f64,
            contig.flank_balance,
            contig.weight as f64,
        ]
    }

    let left_score = score(left);
    let right_score = score(right);
    for (left_value, right_value) in left_score.iter().zip(right_score.iter()) {
        let ordering = left_value.total_cmp(right_value);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::seq::encode_kmer;
    fn info(depth: i64) -> KmerInfo {
        KmerInfo {
            depth,
            position: 100,
            is_reverse: false,
            reference_weight: 0,
        }
    }

    #[test]
    fn backbone_commits_long_arm_and_discards_sibling() {
        let k = 4;
        let entries = [
            ("AAAA", 5),
            ("AAAC", 2),
            ("AACG", 2),
            ("ACGT", 2),
            ("CGTT", 2),
            ("GTTT", 2),
            ("AAAG", 20),
            ("AAGA", 20),
            ("AGAC", 20),
        ];
        let graph: HashMap<_, _> = entries
            .into_iter()
            .map(|(sequence, depth)| (encode_kmer(sequence.as_bytes()).unwrap(), info(depth)))
            .collect();
        let seed = encode_kmer(b"AAAA").unwrap();
        let (paths, visited, _, _) = walk_backbone(&graph, seed, k, 100, 24);
        let extension: Vec<u8> = paths[0].bases.iter().map(|base| bits_base(*base)).collect();
        assert_eq!(extension, b"CGTTT");
        assert!(visited.contains(&encode_kmer(b"AAAC").unwrap()));
        assert!(!visited.contains(&encode_kmer(b"AAAG").unwrap()));
    }

    #[test]
    fn backbone_cycle_stops() {
        let graph: HashMap<_, _> = ["AAA", "AAC", "ACA", "CAA"]
            .into_iter()
            .map(|sequence| (encode_kmer(sequence.as_bytes()).unwrap(), info(5)))
            .collect();
        let seed = encode_kmer(b"AAA").unwrap();
        let (paths, visited, _, _) = walk_backbone(&graph, seed, 3, 100, 24);
        assert_eq!(paths[0].bases.len(), 3);
        assert_eq!(visited.len(), 4);
    }
}
