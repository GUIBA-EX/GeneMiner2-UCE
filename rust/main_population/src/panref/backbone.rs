use super::dbg::{unitig_edges, KmerCounter, Unitig};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Evidence used only by PanRefV2 to rank a local graph path.
pub(crate) struct PathEvidence<'a> {
    pub(crate) sample_support: &'a [u32],
    pub(crate) pe_support: &'a [u64],
    pub(crate) depth_stability: &'a [u64],
    pub(crate) edge_support: &'a HashMap<(usize, usize), u64>,
    pub(crate) edges: &'a [(usize, usize)],
    pub(crate) k: usize,
}

pub(crate) struct ResolvedBackbone {
    pub(crate) sequence: Vec<u8>,
    /// Nodes are stored in the graph-forward orientation. `reversed` states
    /// whether `sequence` is their reverse complement.
    pub(crate) nodes: Vec<usize>,
    pub(crate) reversed: bool,
}

pub(crate) fn assemble_backbone(
    reads_path: &Path,
    baits: &[Vec<u8>],
) -> Result<Option<(Vec<u8>, u64)>, String> {
    let mut counter = KmerCounter::new(31).expect("valid fixed k-mer size");
    stream_interleaved_sequences(reads_path, |sequence| counter.add_read(sequence))?;
    let unitigs = counter.into_unitigs(2);
    let pe_support = pair_support(reads_path, &unitigs)?;
    let edges = unitig_edges(&unitigs, 31);
    let Some((best_id, best)) = unitigs.iter().enumerate().max_by_key(|(id, unitig)| {
        (
            anchor_score(&unitig.sequence, baits),
            pe_support[*id],
            unitig.kmer_count,
            unitig.sequence.len(),
        )
    }) else {
        return Ok(None);
    };
    if anchor_score(&best.sequence, baits).0 == 0 {
        return Ok(None);
    }
    let sequence = extend_backbone(best_id, &unitigs, &edges, &pe_support, baits, 31);
    let reverse = reverse_complement(&sequence);
    let sequence = if anchor_score(&reverse, baits) > anchor_score(&sequence, baits) {
        reverse
    } else {
        sequence
    };
    Ok(Some((sequence, pe_support[best_id])))
}

/// Resolve a bait-anchored coordinate backbone from pre-built unitigs.
///
/// Ranking order is intentionally biological rather than depth-only:
/// bait agreement, number of supporting samples, PE/read-path links, stable
/// per-sample depth, then length as a deterministic final tie-breaker.
pub(crate) fn assemble_backbone_from_unitigs(
    unitigs: &[Unitig],
    baits: &[Vec<u8>],
    evidence: PathEvidence<'_>,
) -> Option<ResolvedBackbone> {
    let (sequence, nodes) = if let Some(path) = resolve_global_path(unitigs, baits, &evidence) {
        path
    } else {
        let best_id = unitigs
            .iter()
            .enumerate()
            .max_by_key(|(id, unitig)| path_score(*id, unitig, baits, &evidence))
            .map(|(id, _)| id)?;
        extend_backbone_with_evidence(
            best_id,
            unitigs,
            evidence.edges,
            baits,
            &evidence,
            evidence.k,
        )
    };
    if anchor_score(&sequence, baits).0 == 0 {
        return None;
    }
    let reverse = reverse_complement(&sequence);
    let (sequence, reversed) = if anchor_score(&reverse, baits) > anchor_score(&sequence, baits) {
        (reverse, true)
    } else {
        (sequence, false)
    };
    Some(ResolvedBackbone {
        sequence,
        nodes,
        reversed,
    })
}

type GlobalScore = (u64, u64, u32, u64, u64, u64);

fn extend_score(left: GlobalScore, right: GlobalScore, edge_support: u64) -> GlobalScore {
    (
        left.0 + right.0,
        left.1 + right.1,
        left.2.min(right.2),
        left.3 + right.3 + edge_support,
        left.4.min(right.4),
        left.5 + right.5,
    )
}

fn resolve_global_path(
    unitigs: &[Unitig],
    baits: &[Vec<u8>],
    evidence: &PathEvidence<'_>,
) -> Option<(Vec<u8>, Vec<usize>)> {
    if unitigs.is_empty() {
        return None;
    }
    let mut outgoing = vec![Vec::new(); unitigs.len()];
    let mut indegree = vec![0_usize; unitigs.len()];
    for &(from, to) in evidence.edges {
        if from < unitigs.len() && to < unitigs.len() {
            outgoing[from].push(to);
            indegree[to] += 1;
        }
    }
    let mut queue = std::collections::VecDeque::new();
    for (id, &degree) in indegree.iter().enumerate() {
        if degree == 0 {
            queue.push_back(id);
        }
    }
    let mut order = Vec::with_capacity(unitigs.len());
    while let Some(node) = queue.pop_front() {
        order.push(node);
        for &next in &outgoing[node] {
            indegree[next] -= 1;
            if indegree[next] == 0 {
                queue.push_back(next);
            }
        }
    }
    if order.len() != unitigs.len() {
        return None;
    }
    let anchors = unitigs
        .iter()
        .map(|unitig| anchor_score(&unitig.sequence, baits))
        .collect::<Vec<_>>();
    let node_score = |id: usize| -> GlobalScore {
        (
            anchors[id].0 as u64,
            anchors[id].1 as u64,
            evidence.sample_support[id],
            evidence.pe_support[id],
            evidence.depth_stability[id],
            unitigs[id].sequence.len() as u64,
        )
    };
    let mut best = (0..unitigs.len())
        .map(|id| (node_score(id), vec![id]))
        .collect::<Vec<_>>();
    for from in order {
        let current = best[from].clone();
        for &to in &outgoing[from] {
            let edge = evidence
                .edge_support
                .get(&(from, to))
                .copied()
                .unwrap_or_default();
            if edge == 0 {
                continue;
            }
            let candidate = extend_score(current.0, node_score(to), edge);
            if candidate > best[to].0 {
                let mut nodes = current.1.clone();
                nodes.push(to);
                best[to] = (candidate, nodes);
            }
        }
    }
    let (_, nodes) = best
        .into_iter()
        .filter(|(score, _)| score.0 > 0)
        .max_by_key(|entry| entry.0)?;
    let mut sequence = unitigs[*nodes.first()?].sequence.clone();
    for &id in nodes.iter().skip(1) {
        sequence.extend_from_slice(&unitigs[id].sequence[evidence.k - 1..]);
    }
    Some((sequence, nodes))
}

fn path_score(
    id: usize,
    unitig: &Unitig,
    baits: &[Vec<u8>],
    evidence: &PathEvidence<'_>,
) -> ((usize, usize), u32, u64, u64, usize) {
    (
        anchor_score(&unitig.sequence, baits),
        evidence.sample_support[id],
        evidence.pe_support[id],
        evidence.depth_stability[id],
        unitig.sequence.len(),
    )
}

fn extend_backbone_with_evidence(
    root: usize,
    unitigs: &[Unitig],
    edges: &[(usize, usize)],
    baits: &[Vec<u8>],
    evidence: &PathEvidence<'_>,
    k: usize,
) -> (Vec<u8>, Vec<usize>) {
    let node_score = |id: usize| path_score(id, &unitigs[id], baits, evidence);
    let transition_score = |from: usize, to: usize| {
        let (_, breadth, _, stability, length) = node_score(to);
        (
            anchor_score(&unitigs[to].sequence, baits),
            breadth,
            evidence
                .edge_support
                .get(&(from, to))
                .copied()
                .unwrap_or_default(),
            stability,
            length,
        )
    };
    let mut seen = std::collections::HashSet::from([root]);
    let mut left = Vec::new();
    let mut current = root;
    while let Some(previous) = edges
        .iter()
        .filter_map(|(from, to)| (*to == current && !seen.contains(from)).then_some(*from))
        .max_by_key(|id| transition_score(*id, current))
    {
        seen.insert(previous);
        left.push(previous);
        current = previous;
    }
    left.reverse();
    let mut right = Vec::new();
    current = root;
    while let Some(next) = edges
        .iter()
        .filter_map(|(from, to)| (*from == current && !seen.contains(to)).then_some(*to))
        .max_by_key(|id| transition_score(current, *id))
    {
        seen.insert(next);
        right.push(next);
        current = next;
    }
    let mut path = left;
    path.push(root);
    path.extend(right);
    let mut sequence = unitigs[path[0]].sequence.clone();
    for &id in path.iter().skip(1) {
        sequence.extend_from_slice(&unitigs[id].sequence[k - 1..]);
    }
    (sequence, path)
}

fn extend_backbone(
    root: usize,
    unitigs: &[Unitig],
    edges: &[(usize, usize)],
    pe: &[u64],
    baits: &[Vec<u8>],
    k: usize,
) -> Vec<u8> {
    let score = |id: usize| {
        (
            anchor_score(&unitigs[id].sequence, baits),
            pe[id],
            unitigs[id].sequence.len(),
        )
    };
    let mut seen = std::collections::HashSet::from([root]);
    let mut left = Vec::new();
    let mut current = root;
    while let Some(previous) = edges
        .iter()
        .filter_map(|(from, to)| (*to == current && !seen.contains(from)).then_some(*from))
        .max_by_key(|id| score(*id))
    {
        seen.insert(previous);
        left.push(previous);
        current = previous;
    }
    left.reverse();
    let mut right = Vec::new();
    current = root;
    while let Some(next) = edges
        .iter()
        .filter_map(|(from, to)| (*from == current && !seen.contains(to)).then_some(*to))
        .max_by_key(|id| score(*id))
    {
        seen.insert(next);
        right.push(next);
        current = next;
    }
    let mut path = left;
    path.push(root);
    path.extend(right);
    let mut sequence = unitigs[path[0]].sequence.clone();
    for id in path.into_iter().skip(1) {
        sequence.extend_from_slice(&unitigs[id].sequence[k - 1..]);
    }
    sequence
}

fn pair_support(path: &Path, unitigs: &[super::dbg::Unitig]) -> Result<Vec<u64>, String> {
    let mut support = vec![0; unitigs.len()];
    stream_interleaved_pairs(path, |first, second| {
        let left = unique_unitig(first, unitigs);
        let right = unique_unitig(second, unitigs);
        if let (Some(left), Some(right)) = (left, right) {
            if left != right {
                support[left] += 1;
                support[right] += 1;
            }
        }
    })?;
    Ok(support)
}

fn unique_unitig(read: &[u8], unitigs: &[super::dbg::Unitig]) -> Option<usize> {
    if read.len() < 21 {
        return None;
    }
    let mut found = None;
    for (id, unitig) in unitigs.iter().enumerate() {
        if read
            .windows(21)
            .any(|seed| unitig.sequence.windows(21).any(|target| target == seed))
        {
            if found.is_some() {
                return None;
            }
            found = Some(id);
        }
    }
    found
}

fn stream_interleaved_pairs(
    path: &Path,
    mut consume: impl FnMut(&[u8], &[u8]),
) -> Result<(), String> {
    let file = File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    while let Some((header, first)) = read_fastq_sequence(&mut reader)? {
        if header.ends_with("/0") {
            continue;
        }
        let Some((_, second)) = read_fastq_sequence(&mut reader)? else {
            return Err("interleaved recruited FASTQ has an unpaired mate".into());
        };
        consume(&first, &second);
    }
    Ok(())
}

fn read_fastq_sequence(reader: &mut dyn BufRead) -> Result<Option<(String, Vec<u8>)>, String> {
    let mut header = String::new();
    if reader.read_line(&mut header).map_err(|e| e.to_string())? == 0 {
        return Ok(None);
    }
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    for line in [&mut sequence, &mut plus, &mut quality] {
        if reader.read_line(line).map_err(|e| e.to_string())? == 0 {
            return Err("truncated recruited FASTQ".into());
        }
    }
    if !header.starts_with('@') || !plus.starts_with('+') {
        return Err("invalid recruited FASTQ".into());
    }
    Ok(Some((
        header.trim().to_string(),
        sequence.trim().as_bytes().to_ascii_uppercase(),
    )))
}

fn stream_interleaved_sequences(path: &Path, mut consume: impl FnMut(&[u8])) -> Result<(), String> {
    let file = File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    loop {
        let mut header = String::new();
        if reader.read_line(&mut header).map_err(|e| e.to_string())? == 0 {
            break;
        }
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut quality = String::new();
        for line in [&mut sequence, &mut plus, &mut quality] {
            if reader.read_line(line).map_err(|e| e.to_string())? == 0 {
                return Err("truncated recruited FASTQ".into());
            }
        }
        if !header.starts_with('@') || !plus.starts_with('+') {
            return Err("invalid recruited FASTQ".into());
        }
        consume(sequence.trim().as_bytes());
    }
    Ok(())
}

fn anchor_score(sequence: &[u8], baits: &[Vec<u8>]) -> (usize, usize) {
    let mut best = 0;
    let mut anchored = 0;
    for bait in baits {
        if sequence.len() < 21 || bait.len() < 21 {
            continue;
        }
        let bait_kmers: HashSet<&[u8]> = bait.windows(21).collect();
        let score = sequence
            .windows(21)
            .filter(|kmer| bait_kmers.contains(kmer))
            .count();
        best = best.max(score);
        anchored += usize::from(score > 0);
    }
    (anchored, best)
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

#[cfg(test)]
mod tests {
    use super::{assemble_backbone, pair_support, path_score, resolve_global_path, PathEvidence};
    use crate::panref::dbg::Unitig;
    use std::collections::HashMap;
    #[test]
    fn chooses_and_orients_a_bait_supported_unitig() {
        let path =
            std::env::temp_dir().join(format!("gm2-panref-backbone-{}.fq", std::process::id()));
        std::fs::write(&path, "@x/1\nAACCGGTTAACCGGTTAACCGGTTAACCGGTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n@x/2\nAACCGGTTAACCGGTTAACCGGTTAACCGGTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n").unwrap();
        let bait = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTT";
        assert_eq!(
            assemble_backbone(&path, &[bait.to_vec()]).unwrap(),
            Some((bait.to_vec(), 0))
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn counts_unique_links_between_distinct_unitigs() {
        let path = std::env::temp_dir().join(format!("gm2-panref-links-{}.fq", std::process::id()));
        let left = b"AAAAAAAAAAAAAAAAAAAAACCCCCCCCCCC";
        let right = b"GGGGGGGGGGGGGGGGGGGGTTTTTTTTTTT";
        std::fs::write(
            &path,
            format!(
                "@x/1\n{}\n+\n{}\n@x/2\n{}\n+\n{}\n",
                String::from_utf8_lossy(left),
                "F".repeat(left.len()),
                String::from_utf8_lossy(right),
                "F".repeat(right.len())
            ),
        )
        .unwrap();
        let unitigs = vec![
            Unitig {
                sequence: left.to_vec(),
                kmer_count: 1,
            },
            Unitig {
                sequence: right.to_vec(),
                kmer_count: 1,
            },
        ];
        assert_eq!(pair_support(&path, &unitigs).unwrap(), vec![1, 1]);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn preserves_reverse_orientation_for_graph_path_writers() {
        let sequence = b"AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCC".to_vec();
        let reverse = super::reverse_complement(&sequence);
        let mut bait = sequence[..21].to_vec();
        bait.extend_from_slice(&reverse);
        let unitigs = [Unitig {
            sequence,
            kmer_count: 1,
        }];
        let evidence = PathEvidence {
            sample_support: &[1],
            pe_support: &[0],
            depth_stability: &[0],
            edge_support: &HashMap::new(),
            edges: &[],
            k: 31,
        };
        let resolved = super::assemble_backbone_from_unitigs(&unitigs, &[bait], evidence).unwrap();
        assert!(resolved.reversed);
        assert_eq!(resolved.nodes, vec![0]);
        assert_eq!(resolved.sequence, reverse);
    }

    #[test]
    fn global_path_does_not_cross_an_unsupported_graph_edge() {
        let unitigs = [
            Unitig {
                sequence: vec![b'A'; 31],
                kmer_count: 1,
            },
            Unitig {
                sequence: vec![b'C'; 31],
                kmer_count: 1,
            },
        ];
        let edges = [(0, 1)];
        let evidence = PathEvidence {
            sample_support: &[2, 2],
            pe_support: &[0, 0],
            depth_stability: &[10, 10],
            edge_support: &HashMap::new(),
            edges: &edges,
            k: 31,
        };
        let (_, nodes) = resolve_global_path(&unitigs, &[vec![b'A'; 31]], &evidence).unwrap();
        assert_eq!(nodes, vec![0]);
    }

    #[test]
    fn global_path_uses_sample_bottleneck_instead_of_node_count() {
        let unitigs = [
            Unitig {
                sequence: vec![b'A'; 31],
                kmer_count: 1,
            },
            Unitig {
                sequence: vec![b'C'; 31],
                kmer_count: 1,
            },
            Unitig {
                sequence: vec![b'G'; 31],
                kmer_count: 1,
            },
            Unitig {
                sequence: vec![b'T'; 31],
                kmer_count: 1,
            },
        ];
        let edges = [(0, 1), (0, 2), (2, 3)];
        let edge_support = HashMap::from([((0, 1), 1), ((0, 2), 1), ((2, 3), 1)]);
        let evidence = PathEvidence {
            sample_support: &[10, 10, 6, 6],
            pe_support: &[0; 4],
            depth_stability: &[10; 4],
            edge_support: &edge_support,
            edges: &edges,
            k: 31,
        };
        let (_, nodes) = resolve_global_path(&unitigs, &[vec![b'A'; 31]], &evidence).unwrap();
        assert_eq!(nodes, vec![0, 1]);
    }

    #[test]
    fn sample_breadth_outranks_pe_depth_and_length() {
        let unitigs = [
            Unitig {
                sequence: b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
                kmer_count: 1,
            },
            Unitig {
                sequence: b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),
                kmer_count: 1,
            },
        ];
        let evidence = PathEvidence {
            sample_support: &[2, 8],
            pe_support: &[100, 1],
            depth_stability: &[10_000, 1],
            edge_support: &HashMap::new(),
            edges: &[],
            k: 31,
        };
        let baits: Vec<Vec<u8>> = Vec::new();
        assert!(
            path_score(1, &unitigs[1], &baits, &evidence)
                > path_score(0, &unitigs[0], &baits, &evidence)
        );
    }
}
