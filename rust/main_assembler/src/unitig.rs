use crate::model::KmerInfo;
use crate::seq::{bits_base, decode_kmer, kmer_mask};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[derive(Debug)]
pub struct Unitig {
    pub kmers: Vec<u128>,
    pub sequence: Vec<u8>,
}
#[derive(Debug)]
pub struct UnitigGraph {
    pub unitigs: Vec<Unitig>,
    pub edges: Vec<(usize, usize)>,
}

fn successors(graph: &HashMap<u128, KmerInfo>, node: u128, k: usize) -> Vec<u128> {
    let prefix = (node & kmer_mask(k - 1)) << 2;
    (0..4)
        .map(|b| prefix | b)
        .filter(|n| graph.contains_key(n))
        .collect()
}

pub fn build_unitig_graph(graph: &HashMap<u128, KmerInfo>, k: usize) -> UnitigGraph {
    let mut indegree: HashMap<u128, usize> = graph.keys().map(|n| (*n, 0)).collect();
    for node in graph.keys() {
        for next in successors(graph, *node, k) {
            *indegree.entry(next).or_default() += 1;
        }
    }
    let mut order: Vec<u128> = graph.keys().copied().collect();
    order.sort_unstable_by_key(|n| (indegree[n] == 1 && successors(graph, *n, k).len() == 1, *n));
    let mut assigned = HashSet::new();
    let mut unitigs = Vec::new();
    for start in order {
        if assigned.contains(&start) {
            continue;
        }
        let mut path = vec![start];
        assigned.insert(start);
        loop {
            let next = successors(graph, *path.last().unwrap(), k);
            if next.len() != 1 {
                break;
            }
            let candidate = next[0];
            if indegree[&candidate] != 1 || assigned.contains(&candidate) {
                break;
            }
            assigned.insert(candidate);
            path.push(candidate);
        }
        let mut sequence = decode_kmer(path[0], k);
        sequence.extend(path.iter().skip(1).map(|n| bits_base((n & 3) as u8)));
        unitigs.push(Unitig {
            kmers: path,
            sequence,
        });
    }
    let mut owner = HashMap::new();
    for (id, unitig) in unitigs.iter().enumerate() {
        for node in &unitig.kmers {
            owner.insert(*node, id);
        }
    }
    let mut edges = HashSet::new();
    for (id, unitig) in unitigs.iter().enumerate() {
        for next in successors(graph, *unitig.kmers.last().unwrap(), k) {
            if let Some(target) = owner.get(&next) {
                edges.insert((id, *target));
            }
        }
    }
    let mut edges: Vec<_> = edges.into_iter().collect();
    edges.sort_unstable();
    UnitigGraph { unitigs, edges }
}

pub fn write_graphs(
    output: &Path,
    key: &str,
    graph: &HashMap<u128, KmerInfo>,
    k: usize,
    gfa: bool,
    dot: bool,
) -> io::Result<()> {
    let compact = build_unitig_graph(graph, k);
    fs::create_dir_all(output)?;
    if gfa {
        let mut out = BufWriter::new(File::create(output.join(format!("{key}.gfa")))?);
        writeln!(out, "H\tVN:Z:1.0")?;
        for (id, unitig) in compact.unitigs.iter().enumerate() {
            writeln!(
                out,
                "S\tu{id}\t{}\tLN:i:{}",
                String::from_utf8_lossy(&unitig.sequence),
                unitig.sequence.len()
            )?;
        }
        for (from, to) in &compact.edges {
            writeln!(out, "L\tu{from}\t+\tu{to}\t+\t{}M", k.saturating_sub(1))?;
        }
    }
    if dot {
        let mut out = BufWriter::new(File::create(output.join(format!("{key}.dot")))?);
        writeln!(out, "digraph assembly {{")?;
        for (id, unitig) in compact.unitigs.iter().enumerate() {
            writeln!(
                out,
                "  u{id} [label=\"u{id} len={}\"];",
                unitig.sequence.len()
            )?;
        }
        for (from, to) in &compact.edges {
            writeln!(out, "  u{from} -> u{to};")?;
        }
        writeln!(out, "}}")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_kmer;
    fn info() -> KmerInfo {
        KmerInfo {
            depth: 1,
            position: 0,
            is_reverse: false,
            reference_weight: 0,
        }
    }
    #[test]
    fn compresses_linear_path() {
        let graph = ["ATG", "TGC", "GCC", "CCA"]
            .into_iter()
            .map(|s| (encode_kmer(s.as_bytes()).unwrap(), info()))
            .collect();
        let compact = build_unitig_graph(&graph, 3);
        assert_eq!(compact.unitigs.len(), 1);
        assert_eq!(compact.unitigs[0].sequence, b"ATGCCA");
    }
    #[test]
    fn writes_gfa_and_dot() {
        let graph = ["ATG", "TGC", "GCC", "CCA"]
            .into_iter()
            .map(|s| (encode_kmer(s.as_bytes()).unwrap(), info()))
            .collect();
        let dir = std::env::temp_dir().join(format!("gm2-unitig-{}", std::process::id()));
        write_graphs(&dir, "locus", &graph, 3, true, true).unwrap();
        let gfa = std::fs::read_to_string(dir.join("locus.gfa")).unwrap();
        let dot = std::fs::read_to_string(dir.join("locus.dot")).unwrap();
        assert!(gfa.contains("S\tu0\tATGCCA"));
        assert!(dot.contains("digraph assembly"));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
