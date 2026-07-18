use crate::fasta::{read_fasta, write_fasta, FastaRecord};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter};
use std::path::Path;

fn pair_is_similar(left: &[u8], right: &[u8], identity_threshold: f64) -> bool {
    let mut overlap = 0_usize;
    let mut identical = 0_usize;

    for (&left, &right) in left.iter().zip(right) {
        let left = left.to_ascii_uppercase();
        let right = right.to_ascii_uppercase();
        if matches!(left, b'-' | b'?') || matches!(right, b'-' | b'?') {
            continue;
        }
        overlap += 1;
        identical += usize::from(left == right);
    }

    // 原实现要求至少 7 个可比较碱基，以降低偶然连边的概率。
    overlap > 6 && identical as f64 / overlap as f64 >= identity_threshold
}

fn find_bridges(adjacency: &[Vec<usize>]) -> HashSet<(usize, usize)> {
    fn visit(
        node: usize,
        parent: usize,
        adjacency: &[Vec<usize>],
        discovery: &mut [usize],
        low: &mut [usize],
        clock: &mut usize,
        bridges: &mut HashSet<(usize, usize)>,
    ) {
        *clock += 1;
        discovery[node] = *clock;
        low[node] = *clock;

        for &next in &adjacency[node] {
            if discovery[next] == 0 {
                visit(next, node, adjacency, discovery, low, clock, bridges);
                low[node] = low[node].min(low[next]);
                if low[next] > discovery[node] {
                    bridges.insert((node.min(next), node.max(next)));
                }
            } else if next != parent {
                low[node] = low[node].min(discovery[next]);
            }
        }
    }

    let count = adjacency.len();
    let mut discovery = vec![0; count];
    let mut low = vec![0; count];
    let mut clock = 0;
    let mut bridges = HashSet::new();
    for node in 0..count {
        if discovery[node] == 0 {
            visit(
                node,
                usize::MAX,
                adjacency,
                &mut discovery,
                &mut low,
                &mut clock,
                &mut bridges,
            );
        }
    }
    bridges
}

#[derive(Debug)]
struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(count: usize) -> Self {
        Self {
            parent: (0..count).collect(),
            rank: vec![0; count],
        }
    }

    fn find(&mut self, value: usize) -> usize {
        if self.parent[value] != value {
            self.parent[value] = self.find(self.parent[value]);
        }
        self.parent[value]
    }

    fn merge(&mut self, left: usize, right: usize) {
        let left = self.find(left);
        let right = self.find(right);
        if left == right {
            return;
        }
        if self.rank[left] < self.rank[right] {
            self.parent[left] = right;
        } else {
            self.parent[right] = left;
            if self.rank[left] == self.rank[right] {
                self.rank[left] += 1;
            }
        }
    }
}

fn bridge_connected_components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let bridges = find_bridges(adjacency);
    let mut sets = DisjointSet::new(adjacency.len());
    for (left, neighbours) in adjacency.iter().enumerate() {
        for &right in neighbours {
            if left < right && !bridges.contains(&(left, right)) {
                sets.merge(left, right);
            }
        }
    }

    let mut components: Vec<Vec<usize>> = Vec::new();
    let mut roots = Vec::new();
    for node in 0..adjacency.len() {
        let root = sets.find(node);
        if let Some(position) = roots.iter().position(|&known| known == root) {
            components[position].push(node);
        } else {
            roots.push(root);
            components.push(vec![node]);
        }
    }
    components.sort_by_key(|component| std::cmp::Reverse(component.len()));
    components
}

/// 清理一个 alignment。返回 None 表示该 locus 应被删除。
pub fn clean_records(
    records: &[FastaRecord],
    minimum_sequences: usize,
    maximum_difference: f64,
) -> Option<Vec<FastaRecord>> {
    if records.len() <= 1 {
        return None;
    }

    let mut adjacency = vec![Vec::new(); records.len()];
    let identity_threshold = 1.0 - maximum_difference;
    for left in 0..records.len() - 1 {
        for right in left + 1..records.len() {
            if pair_is_similar(
                &records[left].sequence,
                &records[right].sequence,
                identity_threshold,
            ) {
                adjacency[left].push(right);
                adjacency[right].push(left);
            }
        }
    }

    let minimum_sequences = minimum_sequences.max(2);
    let components: Vec<Vec<usize>> = bridge_connected_components(&adjacency)
        .into_iter()
        .filter(|component| component.len() >= minimum_sequences)
        .collect();
    if components.is_empty() {
        return None;
    }

    let mut retained = vec![false; records.len()];
    for component in &components {
        for &index in component {
            retained[index] = true;
        }
    }

    let mut output = vec![Vec::new(); records.len()];
    for component in components {
        let shortest = component
            .iter()
            .map(|&index| records[index].sequence.len())
            .min()
            .unwrap_or(0);
        let keep_columns: Vec<bool> = (0..shortest)
            .map(|column| {
                component
                    .iter()
                    .any(|&index| records[index].sequence[column] != b'-')
            })
            .collect();
        let block_length = keep_columns.iter().filter(|&&keep| keep).count();
        let in_component: HashSet<usize> = component.into_iter().collect();

        for index in 0..records.len() {
            if !retained[index] {
                continue;
            }
            if in_component.contains(&index) {
                output[index].extend(
                    records[index]
                        .sequence
                        .iter()
                        .zip(&keep_columns)
                        .filter_map(|(&base, &keep)| keep.then_some(base)),
                );
            } else {
                output[index].extend(std::iter::repeat_n(b'-', block_length));
            }
        }
    }

    Some(
        records
            .iter()
            .enumerate()
            .filter(|(index, _)| retained[*index])
            .map(|(index, record)| FastaRecord {
                name: record.name.clone(),
                sequence: std::mem::take(&mut output[index]),
            })
            .collect(),
    )
}

pub fn clean_file(
    path: &Path,
    minimum_sequences: usize,
    maximum_difference: f64,
) -> io::Result<()> {
    let records = read_fasta(BufReader::new(File::open(path)?))?;
    let Some(records) = clean_records(&records, minimum_sequences, maximum_difference) else {
        fs::remove_file(path)?;
        return Ok(());
    };
    write_fasta(BufWriter::new(File::create(path)?), &records)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(name: &str, sequence: &[u8]) -> FastaRecord {
        FastaRecord {
            name: name.to_string(),
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn removes_singletons_and_all_gap_columns() {
        let records = vec![
            record("a", b"AAA--AAAA"),
            record("b", b"AAA--AAAA"),
            record("c", b"AAA--AAAA"),
            record("singleton", b"TTT--TTTT"),
        ];
        assert_eq!(
            clean_records(&records, 2, 0.0).unwrap(),
            vec![
                record("a", b"AAAAAAA"),
                record("b", b"AAAAAAA"),
                record("c", b"AAAAAAA"),
            ]
        );
    }

    #[test]
    fn concatenates_independent_supported_components() {
        let records = vec![
            record("a1", b"AAAAAAAA"),
            record("a2", b"AAAAAAAA"),
            record("a3", b"AAAAAAAA"),
            record("t1", b"TTTTTTTT"),
            record("t2", b"TTTTTTTT"),
            record("t3", b"TTTTTTTT"),
        ];
        assert_eq!(
            clean_records(&records, 2, 0.0).unwrap(),
            vec![
                record("a1", b"AAAAAAAA--------"),
                record("a2", b"AAAAAAAA--------"),
                record("a3", b"AAAAAAAA--------"),
                record("t1", b"--------TTTTTTTT"),
                record("t2", b"--------TTTTTTTT"),
                record("t3", b"--------TTTTTTTT"),
            ]
        );
    }

    #[test]
    fn requires_seven_comparable_bases() {
        let records = vec![record("a", b"AAAAAA-"), record("b", b"AAAAAA?")];
        assert!(clean_records(&records, 2, 0.0).is_none());
    }
}
