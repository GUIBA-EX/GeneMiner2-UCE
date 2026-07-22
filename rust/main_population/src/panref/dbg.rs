use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Eq, PartialEq)]
pub(crate) struct Unitig {
    pub(crate) sequence: Vec<u8>,
    pub(crate) kmer_count: usize,
}

/// Directed local DBG on bait-oriented reads.  `k <= 31` keeps each k-mer in
/// a u64, avoiding one allocation per observed k-mer.
#[derive(Clone)]
pub(crate) struct KmerCounter {
    k: usize,
    counts: HashMap<u64, u32>,
}

impl KmerCounter {
    pub(crate) fn new(k: usize) -> Option<Self> {
        (3..=31).contains(&k).then_some(Self {
            k,
            counts: HashMap::new(),
        })
    }
    pub(crate) fn add_read(&mut self, read: &[u8]) {
        count_read(read, self.k, &mut self.counts);
    }
    /// True when a read touches an already supported local graph k-mer.
    pub(crate) fn has_supported_kmer(&self, read: &[u8], min_count: u32) -> bool {
        let mut value = 0_u64;
        let mut valid = 0_usize;
        let mask = (1_u64 << (2 * self.k)) - 1;
        for &base in read {
            let Some(bits) = encode_base(base) else {
                value = 0;
                valid = 0;
                continue;
            };
            value = ((value << 2) | bits as u64) & mask;
            valid += 1;
            if valid >= self.k && self.counts.get(&value).copied().unwrap_or_default() >= min_count
            {
                return true;
            }
        }
        false
    }
    pub(crate) fn into_unitigs(self, min_count: u32) -> Vec<Unitig> {
        build_from_counts(self.counts, self.k, min_count)
    }
}

#[cfg(test)]
pub(crate) fn build_unitigs(reads: &[Vec<u8>], k: usize, min_count: u32) -> Vec<Unitig> {
    let Some(mut counter) = KmerCounter::new(k) else {
        return Vec::new();
    };
    for read in reads {
        counter.add_read(read);
    }
    counter.into_unitigs(min_count)
}

fn build_from_counts(counts: HashMap<u64, u32>, k: usize, min_count: u32) -> Vec<Unitig> {
    let nodes: HashSet<u64> = counts
        .into_iter()
        .filter_map(|(kmer, count)| (count >= min_count).then_some(kmer))
        .collect();
    if nodes.is_empty() {
        return Vec::new();
    }
    let mask = (1_u64 << (2 * k)) - 1;
    let suffix_mask = (1_u64 << (2 * (k - 1))) - 1;
    let successors = |node: u64| -> Vec<u64> {
        let prefix = (node & suffix_mask) << 2;
        (0..4)
            .map(|base| prefix | base)
            .filter(|next| nodes.contains(next))
            .collect()
    };
    let mut indegree = nodes
        .iter()
        .map(|node| (*node, 0_usize))
        .collect::<HashMap<_, _>>();
    for node in &nodes {
        for next in successors(*node) {
            *indegree.get_mut(&next).expect("node") += 1;
        }
    }
    let mut starts = nodes
        .iter()
        .filter(|node| indegree[node] != 1)
        .copied()
        .collect::<Vec<_>>();
    starts.sort_unstable();
    let mut visited = HashSet::new();
    let mut unitigs = Vec::new();
    let mut remaining = nodes.iter().copied().collect::<Vec<_>>();
    remaining.sort_unstable();
    for start in starts.into_iter().chain(remaining) {
        if !visited.insert(start) {
            continue;
        }
        let mut path = vec![start];
        loop {
            let next = successors(*path.last().expect("path"));
            if next.len() != 1 || indegree[&next[0]] != 1 || !visited.insert(next[0]) {
                break;
            }
            path.push(next[0]);
        }
        let mut sequence = decode_kmer(path[0], k);
        sequence.extend(path.iter().skip(1).map(|node| decode_base(*node as u8 & 3)));
        unitigs.push(Unitig {
            sequence,
            kmer_count: path.len(),
        });
    }
    unitigs.sort_by(|left, right| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    let _ = mask; // documents the full-kmer width used by the rolling encoder.
    unitigs
}

fn count_read(read: &[u8], k: usize, counts: &mut HashMap<u64, u32>) {
    let mut value = 0_u64;
    let mut valid = 0_usize;
    let mask = (1_u64 << (2 * k)) - 1;
    for &base in read {
        let Some(bits) = encode_base(base) else {
            value = 0;
            valid = 0;
            continue;
        };
        value = ((value << 2) | bits as u64) & mask;
        valid += 1;
        if valid >= k {
            *counts.entry(value).or_default() += 1;
        }
    }
}

fn encode_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}
fn decode_base(bits: u8) -> u8 {
    b"ACGT"[bits as usize]
}
fn decode_kmer(mut value: u64, k: usize) -> Vec<u8> {
    let mut out = vec![b'A'; k];
    for base in out.iter_mut().rev() {
        *base = decode_base((value & 3) as u8);
        value >>= 2;
    }
    out
}

/// Recover directed unitig overlap edges from a compact DBG.  An edge means
/// a real k-1 base overlap, never an inferred gap.
pub(crate) fn unitig_edges(unitigs: &[Unitig], k: usize) -> Vec<(usize, usize)> {
    if k < 2 {
        return Vec::new();
    }
    let overlap = k - 1;
    let mut prefixes: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (to, unitig) in unitigs.iter().enumerate() {
        if unitig.sequence.len() >= overlap {
            prefixes
                .entry(&unitig.sequence[..overlap])
                .or_default()
                .push(to);
        }
    }
    let mut edges = Vec::new();
    for (from, left) in unitigs.iter().enumerate() {
        if left.sequence.len() < overlap {
            continue;
        }
        let suffix = &left.sequence[left.sequence.len() - overlap..];
        if let Some(targets) = prefixes.get(suffix) {
            edges.extend(targets.iter().map(|&to| (from, to)));
        }
    }
    edges.sort_unstable();
    edges
}

#[cfg(test)]
mod tests {
    use super::{build_unitigs, unitig_edges};
    #[test]
    fn compacts_a_supported_linear_path() {
        let reads = vec![
            b"AACCGGTT".to_vec(),
            b"AACCGGTT".to_vec(),
            b"ACCGGTTA".to_vec(),
            b"ACCGGTTA".to_vec(),
        ];
        assert!(build_unitigs(&reads, 4, 2)
            .iter()
            .any(|unitig| unitig.sequence == b"AACCGGTTA"));
    }
    #[test]
    fn canonicalizes_cycle_unitig_start_across_rebuilds() {
        let reads = vec![b"ACGACG".to_vec(), b"ACGACG".to_vec()];
        let expected = build_unitigs(&reads, 3, 2);
        assert_eq!(expected.len(), 1);
        assert_eq!(expected[0].sequence, b"ACGAC");
        for _ in 0..32 {
            assert_eq!(build_unitigs(&reads, 3, 2), expected);
        }
    }

    #[test]
    fn keeps_a_linear_path_maximal_when_an_internal_kmer_sorts_first() {
        let reads = vec![b"CAACC".to_vec(), b"CAACC".to_vec()];
        let unitigs = build_unitigs(&reads, 3, 2);
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].sequence, b"CAACC");
    }

    #[test]
    fn preserves_a_compacted_cycle_as_a_self_edge() {
        let reads = vec![b"ACGACG".to_vec(), b"ACGACG".to_vec()];
        let unitigs = build_unitigs(&reads, 3, 2);
        assert_eq!(unitig_edges(&unitigs, 3), vec![(0, 0)]);
    }

    #[test]
    fn removes_singleton_error_kmers() {
        let reads = vec![
            b"AACCGGTT".to_vec(),
            b"AACCGGTT".to_vec(),
            b"AACCGATT".to_vec(),
        ];
        assert!(build_unitigs(&reads, 4, 2)
            .iter()
            .all(|unitig| !unitig.sequence.windows(4).any(|kmer| kmer == b"CGAT")));
    }
}
