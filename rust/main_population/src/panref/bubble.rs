//! Conservative topology QC for local PanRef graphs.
//!
//! A `simple_bubble` is only reported when every off-backbone successor of a
//! backbone node follows a non-branching, acyclic path and rejoins the same
//! later backbone node. Everything else remains `complex_branch`: this module
//! describes graph ambiguity but never infers alleles or sample haplotypes.

use super::backbone::ResolvedBackbone;
use super::dbg::Unitig;
use std::collections::{HashMap, HashSet};

const MAX_ALTERNATIVE_NODES: usize = 256;

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct BubbleRecord {
    pub(crate) entry: usize,
    pub(crate) exit: Option<usize>,
    pub(crate) status: &'static str,
    pub(crate) alternative_branches: usize,
    pub(crate) alternative_nodes: usize,
    pub(crate) alternative_min_samples: u32,
    pub(crate) alternative_min_edge_support: u64,
    pub(crate) canonical_min_samples: u32,
    pub(crate) canonical_min_edge_support: u64,
}

enum Trace {
    Rejoins {
        exit_position: usize,
        nodes: Vec<usize>,
    },
    Complex {
        nodes: Vec<usize>,
    },
}

/// Find only bounded, topology-unambiguous bubbles adjacent to a resolved
/// backbone. `sample_support` and `edge_support` are aggregate accepted-ledger
/// evidence and are reported for QC; they are not a phasing assertion.
pub(crate) fn inspect_backbone_branches(
    unitigs: &[Unitig],
    edges: &[(usize, usize)],
    backbone: &ResolvedBackbone,
    sample_support: &[u32],
    edge_support: &HashMap<(usize, usize), u64>,
) -> Vec<BubbleRecord> {
    if backbone.nodes.is_empty() {
        return Vec::new();
    }
    let mut outgoing = vec![Vec::new(); unitigs.len()];
    let mut incoming = vec![Vec::new(); unitigs.len()];
    for &(from, to) in edges {
        if from < unitigs.len() && to < unitigs.len() {
            outgoing[from].push(to);
            incoming[to].push(from);
        }
    }
    for neighbours in outgoing.iter_mut().chain(&mut incoming) {
        neighbours.sort_unstable();
        neighbours.dedup();
    }
    let indegree = incoming.iter().map(Vec::len).collect::<Vec<_>>();
    let positions = backbone
        .nodes
        .iter()
        .enumerate()
        .map(|(position, &node)| (node, position))
        .collect::<HashMap<_, _>>();
    let mut records = Vec::new();
    for (entry_position, &entry) in backbone.nodes.iter().enumerate() {
        let canonical_next = backbone.nodes.get(entry_position + 1).copied();
        let alternatives = outgoing[entry]
            .iter()
            .copied()
            .filter(|&node| Some(node) != canonical_next)
            .collect::<Vec<_>>();
        if alternatives.is_empty() {
            continue;
        }
        if canonical_next.is_none() {
            let traces = alternatives
                .iter()
                .map(|&start| {
                    trace_alternative(
                        start,
                        entry,
                        entry_position,
                        &outgoing,
                        &indegree,
                        &positions,
                    )
                })
                .collect::<Vec<_>>();
            let (alternative_min_samples, alternative_min_edge_support) = alternative_metrics(
                entry,
                &alternatives,
                &traces,
                None,
                sample_support,
                edge_support,
            );
            records.push(BubbleRecord {
                entry,
                exit: None,
                status: "terminal_branch",
                alternative_branches: alternatives.len(),
                alternative_nodes: traces.iter().map(trace_nodes).sum(),
                alternative_min_samples,
                alternative_min_edge_support,
                canonical_min_samples: 0,
                canonical_min_edge_support: 0,
            });
            continue;
        }
        let traces = alternatives
            .iter()
            .map(|&start| {
                trace_alternative(
                    start,
                    entry,
                    entry_position,
                    &outgoing,
                    &indegree,
                    &positions,
                )
            })
            .collect::<Vec<_>>();
        let rejoin_position = traces.iter().find_map(|trace| match trace {
            Trace::Rejoins { exit_position, .. } => Some(*exit_position),
            Trace::Complex { .. } => None,
        });
        let simple = rejoin_position.is_some()
            && traces.iter().all(|trace| {
                matches!(trace, Trace::Rejoins { exit_position, .. } if Some(*exit_position) == rejoin_position)
            })
            && canonical_arm_is_isolated(
                backbone,
                entry_position,
                rejoin_position.expect("rejoining alternatives"),
                &outgoing,
                &incoming,
                &traces,
            );
        let alternative_nodes = traces.iter().map(trace_nodes).sum::<usize>();
        let exit_node =
            simple.then(|| backbone.nodes[rejoin_position.expect("simple bubble has exit")]);
        let (alternative_min_samples, alternative_min_edge_support) = alternative_metrics(
            entry,
            &alternatives,
            &traces,
            exit_node,
            sample_support,
            edge_support,
        );
        let (exit, canonical_min_samples, canonical_min_edge_support) = if simple {
            let exit_position = rejoin_position.expect("simple bubble has exit");
            let canonical_nodes = &backbone.nodes[entry_position..=exit_position];
            (
                Some(backbone.nodes[exit_position]),
                canonical_nodes
                    .iter()
                    .map(|&node| sample_support.get(node).copied().unwrap_or_default())
                    .min()
                    .unwrap_or_default(),
                canonical_nodes
                    .windows(2)
                    .map(|edge| {
                        edge_support
                            .get(&(edge[0], edge[1]))
                            .copied()
                            .unwrap_or_default()
                    })
                    .min()
                    .unwrap_or_default(),
            )
        } else {
            (None, 0, 0)
        };
        records.push(BubbleRecord {
            entry,
            exit,
            status: if simple {
                "simple_bubble"
            } else {
                "complex_branch"
            },
            alternative_branches: alternatives.len(),
            alternative_nodes,
            alternative_min_samples,
            alternative_min_edge_support,
            canonical_min_samples,
            canonical_min_edge_support,
        });
    }
    records
}

fn canonical_arm_is_isolated(
    backbone: &ResolvedBackbone,
    entry_position: usize,
    exit_position: usize,
    outgoing: &[Vec<usize>],
    incoming: &[Vec<usize>],
    traces: &[Trace],
) -> bool {
    if exit_position <= entry_position {
        return false;
    }
    for position in entry_position + 1..exit_position {
        let node = backbone.nodes[position];
        let previous = backbone.nodes[position - 1];
        let next = backbone.nodes[position + 1];
        if incoming.get(node).map(Vec::as_slice) != Some(&[previous])
            || outgoing.get(node).map(Vec::as_slice) != Some(&[next])
        {
            return false;
        }
    }
    let exit = backbone.nodes[exit_position];
    let mut expected_incoming = HashSet::from([backbone.nodes[exit_position - 1]]);
    for trace in traces {
        let Trace::Rejoins { nodes, .. } = trace else {
            return false;
        };
        expected_incoming.insert(
            nodes
                .last()
                .copied()
                .unwrap_or(backbone.nodes[entry_position]),
        );
    }
    incoming
        .get(exit)
        .is_some_and(|nodes| nodes.iter().copied().collect::<HashSet<_>>() == expected_incoming)
}

fn trace_alternative(
    start: usize,
    entry: usize,
    entry_position: usize,
    outgoing: &[Vec<usize>],
    indegree: &[usize],
    positions: &HashMap<usize, usize>,
) -> Trace {
    let mut nodes = Vec::new();
    let mut seen = HashSet::from([entry]);
    let mut current = start;
    for _ in 0..MAX_ALTERNATIVE_NODES {
        if let Some(&position) = positions.get(&current) {
            return if position > entry_position {
                Trace::Rejoins {
                    exit_position: position,
                    nodes,
                }
            } else {
                Trace::Complex { nodes }
            };
        }
        if !seen.insert(current) || indegree.get(current).copied().unwrap_or_default() != 1 {
            return Trace::Complex { nodes };
        }
        nodes.push(current);
        let Some(children) = outgoing.get(current) else {
            return Trace::Complex { nodes };
        };
        if children.len() != 1 {
            return Trace::Complex { nodes };
        }
        current = children[0];
    }
    Trace::Complex { nodes }
}

fn trace_nodes(trace: &Trace) -> usize {
    match trace {
        Trace::Rejoins { nodes, .. } | Trace::Complex { nodes } => nodes.len(),
    }
}

fn alternative_metrics(
    entry: usize,
    starts: &[usize],
    traces: &[Trace],
    exit: Option<usize>,
    sample_support: &[u32],
    edge_support: &HashMap<(usize, usize), u64>,
) -> (u32, u64) {
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    for (&start, trace) in starts.iter().zip(traces) {
        let interior = match trace {
            Trace::Rejoins { nodes, .. } | Trace::Complex { nodes } => nodes,
        };
        let mut path = Vec::with_capacity(interior.len() + 2);
        path.push(entry);
        path.extend(interior.iter().copied());
        if let Some(exit) = exit {
            path.push(exit);
        }
        if path.len() == 1 {
            path.push(start);
        }
        nodes.extend(path.iter().copied());
        edges.extend(path.windows(2).map(|edge| (edge[0], edge[1])));
    }
    let min_samples = nodes
        .into_iter()
        .map(|node| sample_support.get(node).copied().unwrap_or_default())
        .min()
        .unwrap_or_default();
    let min_edge_support = edges
        .into_iter()
        .map(|edge| edge_support.get(&edge).copied().unwrap_or_default())
        .min()
        .unwrap_or_default();
    (min_samples, min_edge_support)
}

#[cfg(test)]
mod tests {
    use super::{inspect_backbone_branches, BubbleRecord};
    use crate::panref::backbone::ResolvedBackbone;
    use crate::panref::dbg::Unitig;
    use std::collections::HashMap;

    fn unitigs(count: usize) -> Vec<Unitig> {
        (0..count)
            .map(|_| Unitig {
                sequence: b"ACGTACGTACGTACGTACGTA".to_vec(),
                kmer_count: 1,
            })
            .collect()
    }

    #[test]
    fn reports_a_rejoining_non_branching_alternative_as_simple_bubble() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 2)];
        let evidence = HashMap::from([((0, 1), 9), ((1, 2), 9), ((0, 3), 3), ((3, 2), 3)]);
        let records = inspect_backbone_branches(
            &unitigs(4),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[8, 8, 8, 3],
            &evidence,
        );
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0],
            BubbleRecord {
                entry: 0,
                exit: Some(2),
                status: "simple_bubble",
                alternative_branches: 1,
                alternative_nodes: 1,
                alternative_min_samples: 3,
                alternative_min_edge_support: 3,
                canonical_min_samples: 8,
                canonical_min_edge_support: 9,
            }
        );
    }

    #[test]
    fn keeps_complex_alternative_evidence_local_to_each_branch() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 2), (0, 4)];
        let evidence = HashMap::from([((0, 3), 3), ((3, 2), 3), ((0, 4), 4)]);
        let records = inspect_backbone_branches(
            &unitigs(5),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[5; 5],
            &evidence,
        );
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].status, "complex_branch");
        assert_eq!(records[0].alternative_min_edge_support, 3);
    }

    #[test]
    fn rejects_an_alternative_with_an_external_incoming_edge() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 2), (4, 3)];
        let records = inspect_backbone_branches(
            &unitigs(5),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[1; 5],
            &HashMap::new(),
        );
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].status, "complex_branch");
        assert_eq!(records[0].exit, None);
    }

    #[test]
    fn rejects_a_bubble_when_the_canonical_arm_has_a_side_branch() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 2), (1, 4)];
        let records = inspect_backbone_branches(
            &unitigs(5),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[1; 5],
            &HashMap::new(),
        );
        assert_eq!(
            records
                .iter()
                .find(|record| record.entry == 0)
                .unwrap()
                .status,
            "complex_branch"
        );
    }

    #[test]
    fn rejects_a_bubble_with_an_external_edge_into_the_exit() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 2), (4, 2)];
        let records = inspect_backbone_branches(
            &unitigs(5),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[1; 5],
            &HashMap::new(),
        );
        assert_eq!(records[0].status, "complex_branch");
        assert_eq!(records[0].exit, None);
    }

    #[test]
    fn records_an_unselected_backbone_terminal_extension() {
        let edges = vec![(0, 1), (1, 2)];
        let records = inspect_backbone_branches(
            &unitigs(3),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1],
                reversed: false,
            },
            &[4, 4, 2],
            &HashMap::from([((1, 2), 2)]),
        );
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].status, "terminal_branch");
        assert_eq!(records[0].entry, 1);
        assert_eq!(records[0].exit, None);
    }

    #[test]
    fn never_calls_a_secondary_branch_an_allelic_bubble() {
        let edges = vec![(0, 1), (1, 2), (0, 3), (3, 4), (3, 5)];
        let records = inspect_backbone_branches(
            &unitigs(6),
            &edges,
            &ResolvedBackbone {
                sequence: Vec::new(),
                nodes: vec![0, 1, 2],
                reversed: false,
            },
            &[1; 6],
            &HashMap::new(),
        );
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].status, "complex_branch");
        assert_eq!(records[0].exit, None);
    }
}
