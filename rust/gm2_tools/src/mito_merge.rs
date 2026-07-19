use crate::mito_reads::{count_junction_support, read_interleaved_pairs, FastqRecord};
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct MitoContig {
    pub id: String,
    pub sequence: Vec<u8>,
}

#[derive(Clone, Debug)]
struct Component {
    members: Vec<String>,
    sequence: Vec<u8>,
}

#[derive(Clone, Copy, Debug)]
struct Overlap {
    length: usize,
    matches: usize,
}

/// Evidence that a linear contig is an assembler walk around the same cycle
/// more than once (for example, `C + C + prefix(C)`).  This deliberately uses
/// only the contig itself; the GenBank reference is not consulted.
#[derive(Clone, Copy, Debug)]
struct UnrolledCycle {
    period: usize,
    compared_bases: usize,
    matches: usize,
}

#[derive(Clone, Debug)]
pub struct LinkConfig {
    pub minimum_overlap: usize,
    pub minimum_identity: f64,
    pub terminal_window: usize,
    pub link_kmer: usize,
    pub minimum_link_hits: usize,
    pub minimum_pair_support: usize,
    pub bridge_kmer: usize,
    pub bridge_minimum_depth: usize,
    pub maximum_bridge: usize,
    pub minimum_junction_support: usize,
    pub expected_length: usize,
}

fn rc(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn best_overlap(left: &[u8], right: &[u8], minimum: usize, identity: f64) -> Option<Overlap> {
    let maximum = left.len().min(right.len());
    for length in (minimum..=maximum).rev() {
        let matches = left[left.len() - length..]
            .iter()
            .zip(&right[..length])
            .filter(|(a, b)| a == b)
            .count();
        if matches as f64 / length as f64 >= identity {
            return Some(Overlap { length, matches });
        }
    }
    None
}

fn closure_overlap(sequence: &[u8], minimum: usize, identity: f64) -> Option<Overlap> {
    if sequence.len() < minimum * 2 {
        return None;
    }
    for length in (minimum..=sequence.len() / 2).rev() {
        let matches = sequence[sequence.len() - length..]
            .iter()
            .zip(&sequence[..length])
            .filter(|(a, b)| a == b)
            .count();
        if matches as f64 / length as f64 >= identity {
            return Some(Overlap { length, matches });
        }
    }
    None
}

fn unrolled_cycle(sequence: &[u8], k: usize) -> Option<UnrolledCycle> {
    const MINIMUM_PERIOD: usize = 1_000;
    const MINIMUM_IDENTITY: f64 = 0.995;
    const MINIMUM_SECOND_LAP: f64 = 0.90;
    if k == 0 || sequence.len() < MINIMUM_PERIOD * 2 || sequence.len() < k * 2 {
        return None;
    }

    // Repeated k-mers nominate offsets cheaply. Every nominated offset is
    // then tested over the full available overlap, so a local repeat cannot
    // turn a contig into a circular genome.
    let mut anchors = HashMap::<u128, Vec<usize>>::new();
    let mut offsets = HashSet::<usize>::new();
    let lower = sequence.len() / 4;
    let upper = sequence.len() / 2 + 1;
    for (position, value) in kmers(sequence, k).into_iter().enumerate() {
        let positions = anchors.entry(value).or_default();
        for previous in positions.iter().copied() {
            let offset = position - previous;
            if offset >= MINIMUM_PERIOD && offset >= lower && offset <= upper {
                offsets.insert(offset);
            }
        }
        if positions.len() < 4 {
            positions.push(position);
        }
    }

    let mut best = None;
    for period in offsets {
        let compared_bases = sequence.len() - period;
        if compared_bases as f64 / (period as f64) < MINIMUM_SECOND_LAP {
            continue;
        }
        let matches = sequence[..compared_bases]
            .iter()
            .zip(&sequence[period..])
            .filter(|(left, right)| left == right)
            .count();
        if matches as f64 / (compared_bases as f64) < MINIMUM_IDENTITY {
            continue;
        }
        let candidate = UnrolledCycle {
            period,
            compared_bases,
            matches,
        };
        if best.as_ref().is_none_or(|old: &UnrolledCycle| {
            (
                candidate.matches,
                candidate.compared_bases,
                candidate.period,
            ) > (old.matches, old.compared_bases, old.period)
        }) {
            best = Some(candidate);
        }
    }
    best
}
fn merge_sequence(left: &[u8], right: &[u8], overlap: usize) -> Vec<u8> {
    let mut merged = left.to_vec();
    let offset = merged.len() - overlap;
    for (index, base) in right[..overlap].iter().enumerate() {
        if merged[offset + index] != *base {
            merged[offset + index] = b'N';
        }
    }
    merged.extend_from_slice(&right[overlap..]);
    merged
}

fn collapse_overlaps(mut components: Vec<Component>, config: &LinkConfig) -> Vec<Component> {
    loop {
        let mut best: Option<(usize, usize, bool, Overlap)> = None;
        for left in 0..components.len() {
            for right in 0..components.len() {
                if left == right {
                    continue;
                }
                for reverse in [false, true] {
                    let oriented = if reverse {
                        rc(&components[right].sequence)
                    } else {
                        components[right].sequence.clone()
                    };
                    let Some(overlap) = best_overlap(
                        &components[left].sequence,
                        &oriented,
                        config.minimum_overlap,
                        config.minimum_identity,
                    ) else {
                        continue;
                    };
                    if best.as_ref().is_none_or(|(_, _, _, old)| {
                        (overlap.length, overlap.matches) > (old.length, old.matches)
                    }) {
                        best = Some((left, right, reverse, overlap));
                    }
                }
            }
        }
        let Some((left, right, reverse, overlap)) = best else {
            break;
        };
        let oriented = if reverse {
            rc(&components[right].sequence)
        } else {
            components[right].sequence.clone()
        };
        let mut members = components[left].members.clone();
        members.extend(components[right].members.iter().cloned());
        let merged = Component {
            members,
            sequence: merge_sequence(&components[left].sequence, &oriented, overlap.length),
        };
        let high = left.max(right);
        let low = left.min(right);
        components.remove(high);
        components.remove(low);
        components.push(merged);
    }
    components
}

fn encode(sequence: &[u8]) -> Option<u128> {
    let mut value = 0_u128;
    for base in sequence {
        let bits = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        value = (value << 2) | bits;
    }
    Some(value)
}

fn mask(k: usize) -> u128 {
    if k == 64 {
        u128::MAX
    } else {
        (1_u128 << (2 * k)) - 1
    }
}

fn kmers(sequence: &[u8], k: usize) -> Vec<u128> {
    if k == 0 || k > 63 || sequence.len() < k {
        return Vec::new();
    }
    sequence.windows(k).filter_map(encode).collect()
}

fn oriented_sequence(component: &Component, reverse: bool) -> Vec<u8> {
    if reverse {
        rc(&component.sequence)
    } else {
        component.sequence.clone()
    }
}

fn terminal_indexes(
    components: &[Component],
    k: usize,
    window: usize,
) -> (HashMap<u128, Vec<usize>>, HashMap<u128, Vec<usize>>) {
    let mut suffix = HashMap::<u128, Vec<usize>>::new();
    let mut prefix = HashMap::<u128, Vec<usize>>::new();
    for (component, item) in components.iter().enumerate() {
        for reverse in [false, true] {
            let oriented = oriented_sequence(item, reverse);
            let oriented_id = component * 2 + usize::from(reverse);
            let left_end = oriented.len().min(window);
            let right_start = oriented.len().saturating_sub(window);
            for value in kmers(&oriented[..left_end], k) {
                prefix.entry(value).or_default().push(oriented_id);
            }
            for value in kmers(&oriented[right_start..], k) {
                suffix.entry(value).or_default().push(oriented_id);
            }
        }
    }
    (suffix, prefix)
}

fn unique_hit(
    read: &[u8],
    index: &HashMap<u128, Vec<usize>>,
    k: usize,
    minimum_hits: usize,
) -> Option<usize> {
    let mut counts = HashMap::<usize, usize>::new();
    let mut seen = HashSet::new();
    for value in kmers(read, k) {
        if !seen.insert(value) {
            continue;
        }
        if let Some(targets) = index.get(&value) {
            for target in targets {
                *counts.entry(*target).or_default() += 1;
            }
        }
    }
    let mut ranked: Vec<_> = counts.into_iter().collect();
    ranked.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    let (target, hits) = *ranked.first()?;
    if hits < minimum_hits || ranked.get(1).is_some_and(|second| second.1 == hits) {
        return None;
    }
    Some(target)
}

fn add_pair_edge(
    left: &FastqRecord,
    right: &FastqRecord,
    suffix: &HashMap<u128, Vec<usize>>,
    prefix: &HashMap<u128, Vec<usize>>,
    config: &LinkConfig,
    links: &mut HashMap<(usize, usize), usize>,
) {
    let from = unique_hit(
        &left.sequence,
        suffix,
        config.link_kmer,
        config.minimum_link_hits,
    );
    let to = unique_hit(
        &rc(&right.sequence),
        prefix,
        config.link_kmer,
        config.minimum_link_hits,
    );
    if let (Some(from), Some(to)) = (from, to) {
        *links.entry((from, to)).or_default() += 1;
    }
}

fn mate_links(
    components: &[Component],
    pairs: &[(FastqRecord, FastqRecord)],
    config: &LinkConfig,
) -> HashMap<(usize, usize), usize> {
    let (suffix, prefix) = terminal_indexes(components, config.link_kmer, config.terminal_window);
    let mut links = HashMap::new();
    for (first, second) in pairs {
        add_pair_edge(first, second, &suffix, &prefix, config, &mut links);
        add_pair_edge(second, first, &suffix, &prefix, config, &mut links);
    }
    links
}

fn read_graph(pairs: &[(FastqRecord, FastqRecord)], k: usize) -> HashMap<u128, usize> {
    let mut graph = HashMap::new();
    for (first, second) in pairs {
        for read in [&first.sequence, &second.sequence] {
            for sequence in [read.to_vec(), rc(read)] {
                let mut physical = HashSet::new();
                physical.extend(kmers(&sequence, k));
                for value in physical {
                    *graph.entry(value).or_default() += 1;
                }
            }
        }
    }
    graph
}

fn bridge_path(
    left: &[u8],
    right: &[u8],
    graph: &HashMap<u128, usize>,
    config: &LinkConfig,
) -> Option<Vec<u8>> {
    let k = config.bridge_kmer;
    if left.len() < k || right.len() < k {
        return None;
    }
    let start = encode(&left[left.len() - k..])?;
    let goal = encode(&right[..k])?;
    if start == goal {
        return Some(Vec::new());
    }
    let suffix_mask = mask(k - 1);
    let maximum_steps = config.maximum_bridge.saturating_add(k);
    let mut queue = VecDeque::from([start]);
    let mut distance = HashMap::from([(start, 0_usize)]);
    let mut ways = HashMap::from([(start, 1_u8)]);
    let mut parent = HashMap::<u128, (u128, u8)>::new();
    let mut visited_nodes = 0_usize;
    while let Some(node) = queue.pop_front() {
        visited_nodes += 1;
        if visited_nodes > 250_000 {
            return None;
        }
        let depth = distance[&node];
        if depth >= maximum_steps {
            continue;
        }
        let prefix = (node & suffix_mask) << 2;
        for base in 0..4_u8 {
            let next = prefix | base as u128;
            if next != goal && graph.get(&next).copied().unwrap_or(0) < config.bridge_minimum_depth
            {
                continue;
            }
            let candidate_distance = depth + 1;
            match distance.get(&next).copied() {
                None => {
                    distance.insert(next, candidate_distance);
                    ways.insert(next, ways[&node]);
                    parent.insert(next, (node, base));
                    queue.push_back(next);
                }
                Some(old) if old == candidate_distance => {
                    let count = ways.get(&next).copied().unwrap_or(0);
                    ways.insert(next, count.saturating_add(ways[&node]).min(2));
                }
                _ => {}
            }
        }
    }
    if ways.get(&goal).copied()? != 1 {
        return None;
    }
    let mut bases = Vec::new();
    let mut cursor = goal;
    while cursor != start {
        let (previous, base) = *parent.get(&cursor)?;
        bases.push(match base {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
        cursor = previous;
    }
    bases.reverse();
    Some(bases)
}

fn bridge_components(
    mut components: Vec<Component>,
    pairs: &[(FastqRecord, FastqRecord)],
    graph: &HashMap<u128, usize>,
    config: &LinkConfig,
    accepted_links: &mut Vec<(String, String, usize, usize)>,
) -> Vec<Component> {
    loop {
        let links = mate_links(&components, pairs, config);
        let mut ranked: Vec<_> = links.into_iter().collect();
        ranked.sort_by(|left, right| right.1.cmp(&left.1));
        let mut joined = None;
        for ((from, to), support) in ranked {
            if support < config.minimum_pair_support || from / 2 == to / 2 {
                continue;
            }
            let left = oriented_sequence(&components[from / 2], from % 2 == 1);
            let right = oriented_sequence(&components[to / 2], to % 2 == 1);
            let Some(path) = bridge_path(&left, &right, graph, config) else {
                continue;
            };
            let mut sequence = left;
            sequence.extend_from_slice(&path);
            sequence.extend_from_slice(&right[config.bridge_kmer.min(right.len())..]);
            let mut members = components[from / 2].members.clone();
            members.extend(components[to / 2].members.iter().cloned());
            accepted_links.push((
                components[from / 2].members.join(","),
                components[to / 2].members.join(","),
                support,
                path.len().saturating_sub(config.bridge_kmer),
            ));
            joined = Some((from / 2, to / 2, Component { members, sequence }));
            break;
        }
        let Some((left, right, merged)) = joined else {
            break;
        };
        components.remove(left.max(right));
        components.remove(left.min(right));
        components.push(merged);
        components = collapse_overlaps(components, config);
    }
    components
}

pub fn assemble_and_write(
    output: &Path,
    contigs: Vec<MitoContig>,
    paired_reads: &Path,
    config: &LinkConfig,
) -> Result<String, String> {
    fs::create_dir_all(output).map_err(|error| error.to_string())?;
    let pairs = read_interleaved_pairs(paired_reads)?;
    let mut seen = HashSet::new();
    let mut components: Vec<Component> = contigs
        .into_iter()
        .filter(|contig| !contig.sequence.is_empty())
        .filter(|contig| {
            let reverse = rc(&contig.sequence);
            let key = contig.sequence.clone().min(reverse);
            seen.insert(key)
        })
        .map(|contig| Component {
            members: vec![contig.id],
            sequence: contig.sequence,
        })
        .collect();
    let input_contigs = components.len();
    components = collapse_overlaps(components, config);
    let initial_links = mate_links(&components, &pairs, config);
    let candidate_links: Vec<_> = initial_links
        .iter()
        .filter(|(_, support)| **support >= config.minimum_pair_support)
        .map(|((from, to), support)| {
            (
                format!(
                    "{}{}",
                    components[*from / 2].members.join(","),
                    if *from % 2 == 1 { "(-)" } else { "(+)" }
                ),
                format!(
                    "{}{}",
                    components[*to / 2].members.join(","),
                    if *to % 2 == 1 { "(-)" } else { "(+)" }
                ),
                *support,
            )
        })
        .collect();
    let maximum_pair_support = initial_links.values().copied().max().unwrap_or(0);
    let graph = read_graph(&pairs, config.bridge_kmer);
    let mut accepted_links = Vec::new();
    components = bridge_components(components, &pairs, &graph, config, &mut accepted_links);
    components.sort_by(|left, right| right.sequence.len().cmp(&left.sequence.len()));
    let raw_components = components.clone();

    let mut closure = "none";
    let mut closure_overlap_length = 0_usize;
    let mut cycle_period = 0_usize;
    let mut cycle_compared_bases = 0_usize;
    let mut cycle_identity = 0.0_f64;
    if !components.is_empty()
        && components[0].sequence.len() >= config.expected_length.saturating_mul(3) / 5
    {
        if let Some(cycle) = unrolled_cycle(&components[0].sequence, config.link_kmer) {
            components[0].sequence.truncate(cycle.period);
            cycle_period = cycle.period;
            cycle_compared_bases = cycle.compared_bases;
            cycle_identity = cycle.matches as f64 / cycle.compared_bases as f64;
            closure = "unrolled_cycle";
        } else if let Some(overlap) = closure_overlap(
            &components[0].sequence,
            config.minimum_overlap,
            config.minimum_identity,
        ) {
            let new_length = components[0].sequence.len() - overlap.length;
            components[0].sequence.truncate(new_length);
            closure_overlap_length = overlap.length;
            closure = "overlap";
        } else if components.len() == 1 {
            let links = mate_links(&components, &pairs, config);
            let best_self = links
                .into_iter()
                .filter(|((from, to), support)| {
                    from / 2 == to / 2
                        && from % 2 == to % 2
                        && *support >= config.minimum_pair_support
                })
                .max_by_key(|(_, support)| *support);
            if let Some(((from, _), _support)) = best_self {
                let oriented = oriented_sequence(&components[0], from % 2 == 1);
                if let Some(path) = bridge_path(&oriented, &oriented, &graph, config) {
                    if path.len() >= config.bridge_kmer {
                        let addition = path.len() - config.bridge_kmer;
                        let mut closed = oriented;
                        closed.extend_from_slice(&path[..addition]);
                        components[0].sequence = closed;
                        closure = "mate_bridge";
                    }
                }
            }
        }
    }

    let junction_support = if components.len() == 1 && closure != "none" {
        count_junction_support(paired_reads, &components[0].sequence, config.link_kmer)?
    } else {
        0
    };
    let status = if components.len() == 1
        && closure != "none"
        && junction_support >= config.minimum_junction_support
        && !components[0].sequence.contains(&b'N')
    {
        "circular"
    } else if components.len() == 1 {
        "linear_single_contig"
    } else {
        "partial_multi_contig"
    };

    let mut fasta = File::create(output.join("mitochondrial_assembly.fasta"))
        .map_err(|error| error.to_string())?;
    for (index, component) in components.iter().enumerate() {
        writeln!(
            fasta,
            ">mito_contig_{} status={} length={} members={}",
            index + 1,
            if index == 0 { status } else { "alternative" },
            component.sequence.len(),
            component.members.len()
        )
        .map_err(|error| error.to_string())?;
        for line in component.sequence.chunks(80) {
            writeln!(fasta, "{}", String::from_utf8_lossy(line))
                .map_err(|error| error.to_string())?;
        }
    }
    let mut raw_fasta = File::create(output.join("mitochondrial_assembly_raw.fasta"))
        .map_err(|error| error.to_string())?;
    for (index, component) in raw_components.iter().enumerate() {
        writeln!(
            raw_fasta,
            ">mito_raw_contig_{} length={} members={}",
            index + 1,
            component.sequence.len(),
            component.members.len()
        )
        .map_err(|error| error.to_string())?;
        for line in component.sequence.chunks(80) {
            writeln!(raw_fasta, "{}", String::from_utf8_lossy(line))
                .map_err(|error| error.to_string())?;
        }
    }
    let mut links = File::create(output.join("mitochondrial_mate_links.tsv"))
        .map_err(|error| error.to_string())?;
    writeln!(
        links,
        "from\tto\tpair_support\tstatus\tresolved_bridge_bases"
    )
    .map_err(|error| error.to_string())?;
    for (from, to, support) in &candidate_links {
        writeln!(links, "{from}\t{to}\t{support}\tcandidate\t")
            .map_err(|error| error.to_string())?;
    }
    for (from, to, support, bases) in &accepted_links {
        writeln!(links, "{from}\t{to}\t{support}\tresolved\t{bases}")
            .map_err(|error| error.to_string())?;
    }
    let longest = components
        .first()
        .map_or(0, |component| component.sequence.len());
    let mut summary = File::create(output.join("mitochondrial_assembly_summary.tsv"))
        .map_err(|error| error.to_string())?;
    writeln!(summary, "metric\tvalue").map_err(|error| error.to_string())?;
    writeln!(summary, "status\t{status}").map_err(|error| error.to_string())?;
    writeln!(summary, "input_read_pairs\t{}", pairs.len()).map_err(|error| error.to_string())?;
    writeln!(summary, "input_contigs\t{input_contigs}").map_err(|error| error.to_string())?;
    writeln!(summary, "merged_components\t{}", components.len())
        .map_err(|error| error.to_string())?;
    writeln!(summary, "longest_contig\t{longest}").map_err(|error| error.to_string())?;
    writeln!(summary, "resolved_mate_links\t{}", accepted_links.len())
        .map_err(|error| error.to_string())?;
    writeln!(summary, "candidate_mate_links\t{}", candidate_links.len())
        .map_err(|error| error.to_string())?;
    writeln!(summary, "maximum_pair_support\t{maximum_pair_support}")
        .map_err(|error| error.to_string())?;
    writeln!(summary, "closure_method\t{closure}").map_err(|error| error.to_string())?;
    writeln!(summary, "closure_overlap\t{closure_overlap_length}")
        .map_err(|error| error.to_string())?;
    writeln!(summary, "cycle_period\t{cycle_period}").map_err(|error| error.to_string())?;
    writeln!(summary, "cycle_compared_bases\t{cycle_compared_bases}")
        .map_err(|error| error.to_string())?;
    writeln!(summary, "cycle_identity\t{cycle_identity:.6}").map_err(|error| error.to_string())?;
    writeln!(summary, "junction_read_support\t{junction_support}")
        .map_err(|error| error.to_string())?;
    Ok(status.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> LinkConfig {
        LinkConfig {
            minimum_overlap: 4,
            minimum_identity: 1.0,
            terminal_window: 20,
            link_kmer: 4,
            minimum_link_hits: 1,
            minimum_pair_support: 1,
            bridge_kmer: 4,
            bridge_minimum_depth: 1,
            maximum_bridge: 20,
            minimum_junction_support: 1,
            expected_length: 16,
        }
    }

    #[test]
    fn overlap_merge_uses_reverse_complement_candidate() {
        let components = vec![
            Component {
                members: vec!["a".into()],
                sequence: b"AAAACCCC".to_vec(),
            },
            Component {
                members: vec!["b".into()],
                sequence: b"CCCCAAAA".to_vec(),
            },
        ];
        let merged = collapse_overlaps(components, &config());
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].sequence, b"AAAACCCCAAAA");
    }

    #[test]
    fn unique_read_graph_path_resolves_bridge() {
        let pairs = vec![(
            FastqRecord {
                header: "@x/1".into(),
                sequence: b"AAAACCCCGGGG".to_vec(),
                plus: "+".into(),
                quality: "F".repeat(12),
            },
            FastqRecord {
                header: "@x/2".into(),
                sequence: b"CCCCGGGGTTTT".to_vec(),
                plus: "+".into(),
                quality: "F".repeat(12),
            },
        )];
        let graph = read_graph(&pairs, 4);
        let path = bridge_path(b"AAAACCCC", b"GGGGTTTT", &graph, &config()).unwrap();
        assert!(!path.is_empty());
    }

    #[test]
    fn detects_two_laps_and_a_partial_third_lap() {
        let mut state = 17_u64;
        let cycle: Vec<u8> = (0..1_200)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                match (state >> 32) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            })
            .collect();
        let mut unrolled = cycle.clone();
        unrolled.extend_from_slice(&cycle);
        unrolled.extend_from_slice(&cycle[..29]);
        let detected = unrolled_cycle(&unrolled, 21).unwrap();
        assert_eq!(detected.period, cycle.len());
        assert_eq!(detected.compared_bases, cycle.len() + 29);
        assert_eq!(detected.matches, cycle.len() + 29);
    }
}
