use std::collections::{BTreeSet, HashSet};

#[derive(Clone, Debug)]
pub struct Node {
    pub name: Option<String>,
    pub children: Vec<usize>,
}
#[derive(Clone, Debug)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: usize,
}
#[derive(Clone, Debug)]
pub struct Clade {
    pub node: usize,
    pub parent: Option<usize>,
    pub leaves: Vec<String>,
    pub samples: BTreeSet<String>,
}

fn skip(xs: &[u8], i: &mut usize) {
    while *i < xs.len() && xs[*i].is_ascii_whitespace() {
        *i += 1
    }
}
fn label(xs: &[u8], i: &mut usize) -> String {
    skip(xs, i);
    if *i < xs.len() && xs[*i] == b'\'' {
        *i += 1;
        let mut value = String::new();
        while *i < xs.len() {
            if xs[*i] == b'\'' {
                if *i + 1 < xs.len() && xs[*i + 1] == b'\'' {
                    value.push('\'');
                    *i += 2;
                } else {
                    *i += 1;
                    return value;
                }
            } else {
                value.push(xs[*i] as char);
                *i += 1;
            }
        }
        return value;
    }
    let start = *i;
    while *i < xs.len() && !b",();:".contains(&xs[*i]) {
        *i += 1
    }
    String::from_utf8_lossy(&xs[start..*i]).trim().to_owned()
}

fn newick_label(value: &str) -> String {
    if value.is_empty()
        || value.chars().any(|character| {
            character.is_whitespace() || matches!(character, ',' | ':' | '(' | ')' | ';' | '\'')
        })
    {
        format!("'{}'", value.replace('\'', "''"))
    } else {
        value.to_owned()
    }
}

fn branch(xs: &[u8], i: &mut usize) {
    skip(xs, i);
    if *i < xs.len() && xs[*i] == b':' {
        *i += 1;
        while *i < xs.len() && !b",();".contains(&xs[*i]) {
            *i += 1
        }
    }
}
fn node(xs: &[u8], i: &mut usize, nodes: &mut Vec<Node>) -> Result<usize, String> {
    skip(xs, i);
    let id = nodes.len();
    nodes.push(Node {
        name: None,
        children: vec![],
    });
    if *i >= xs.len() {
        return Err("unexpected end".into());
    }
    if xs[*i] == b'(' {
        *i += 1;
        loop {
            let child = node(xs, i, nodes)?;
            nodes[id].children.push(child);
            skip(xs, i);
            if *i >= xs.len() {
                return Err("unbalanced tree".into());
            }
            if xs[*i] == b',' {
                *i += 1
            } else if xs[*i] == b')' {
                *i += 1;
                break;
            } else {
                return Err("invalid tree separator".into());
            }
        }
        let x = label(xs, i);
        if !x.is_empty() {
            nodes[id].name = Some(x)
        }
    } else {
        let x = label(xs, i);
        if x.is_empty() {
            return Err("missing leaf".into());
        }
        nodes[id].name = Some(x)
    }
    branch(xs, i);
    Ok(id)
}
pub fn parse_newick(text: &str) -> Result<Tree, String> {
    let xs = text.as_bytes();
    let mut i = 0;
    let mut nodes = vec![];
    let root = node(xs, &mut i, &mut nodes)?;
    skip(xs, &mut i);
    if i < xs.len() && xs[i] == b';' {
        i += 1
    }
    skip(xs, &mut i);
    if i != xs.len() {
        return Err("trailing Newick data".into());
    }
    Ok(Tree { nodes, root })
}

fn adjacency(tree: &Tree) -> Vec<Vec<usize>> {
    let mut graph = vec![Vec::new(); tree.nodes.len()];
    for (parent, node) in tree.nodes.iter().enumerate() {
        for &child in &node.children {
            graph[parent].push(child);
            graph[child].push(parent);
        }
    }
    graph
}
fn leaves_from(
    tree: &Tree,
    graph: &[Vec<usize>],
    node: usize,
    parent: Option<usize>,
) -> Vec<String> {
    let children: Vec<_> = graph[node]
        .iter()
        .copied()
        .filter(|next| Some(*next) != parent)
        .collect();
    if children.is_empty() {
        return tree.nodes[node].name.clone().into_iter().collect();
    }
    children
        .into_iter()
        .flat_map(|next| leaves_from(tree, graph, next, Some(node)))
        .collect()
}
fn samples(leaves: &[String]) -> BTreeSet<String> {
    leaves
        .iter()
        .filter_map(|x| x.split('|').next().map(str::to_owned))
        .collect()
}

/// Select maximal, disjoint one-candidate-per-sample clades from an unrooted
/// tree. Every directed edge is considered, so a valid clade on the complement
/// of an arbitrary Newick root is not missed. With outgroups, their unique
/// monophyletic edge is first required and all selected clades must exclude it.
pub fn select_scogs(
    tree: &Tree,
    min_taxa: usize,
    outgroups: &BTreeSet<String>,
) -> Result<Vec<Clade>, String> {
    let graph = adjacency(tree);
    if !outgroups.is_empty() {
        let mut outgroup_edges = 0usize;
        for (left, neighbors) in graph.iter().enumerate() {
            for &right in neighbors {
                if left >= right {
                    continue;
                }
                let left_side = samples(&leaves_from(tree, &graph, left, Some(right)));
                let right_side = samples(&leaves_from(tree, &graph, right, Some(left)));
                if left_side == *outgroups || right_side == *outgroups {
                    outgroup_edges += 1;
                }
            }
        }
        if outgroup_edges != 1 {
            return Err("outgroup_missing_or_not_monophyletic".into());
        }
    }

    let mut candidates = Vec::new();
    if outgroups.is_empty() {
        let leaves = leaves_from(tree, &graph, tree.root, None);
        let taxa = samples(&leaves);
        if leaves.len() == taxa.len() && taxa.len() >= min_taxa {
            candidates.push(Clade {
                node: tree.root,
                parent: None,
                leaves,
                samples: taxa,
            });
        }
    }
    for (node, neighbors) in graph.iter().enumerate() {
        for &parent in neighbors {
            let leaves = leaves_from(tree, &graph, node, Some(parent));
            let taxa = samples(&leaves);
            if leaves.len() == taxa.len() && taxa.len() >= min_taxa && taxa.is_disjoint(outgroups) {
                candidates.push(Clade {
                    node,
                    parent: Some(parent),
                    leaves,
                    samples: taxa,
                });
            }
        }
    }
    candidates.sort_by(|left, right| {
        right
            .samples
            .len()
            .cmp(&left.samples.len())
            .then_with(|| left.leaves.cmp(&right.leaves))
    });
    let mut used = HashSet::new();
    Ok(candidates
        .into_iter()
        .filter_map(|clade| {
            if clade.leaves.iter().all(|leaf| !used.contains(leaf)) {
                used.extend(clade.leaves.iter().cloned());
                Some(clade)
            } else {
                None
            }
        })
        .collect())
}

pub fn all_leaf_names(tree: &Tree) -> Vec<String> {
    tree.nodes
        .iter()
        .filter(|node| node.children.is_empty())
        .filter_map(|node| node.name.clone())
        .collect()
}

pub fn render_clade_samples(tree: &Tree, node: usize, parent: Option<usize>) -> String {
    let graph = adjacency(tree);
    fn render(tree: &Tree, graph: &[Vec<usize>], node: usize, parent: Option<usize>) -> String {
        let children: Vec<_> = graph[node]
            .iter()
            .copied()
            .filter(|next| Some(*next) != parent)
            .collect();
        if children.is_empty() {
            return tree.nodes[node]
                .name
                .as_deref()
                .and_then(|name| name.split('|').next())
                .map(newick_label)
                .unwrap_or_else(|| newick_label(""));
        }
        format!(
            "({})",
            children
                .into_iter()
                .map(|next| render(tree, graph, next, Some(node)))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
    render(tree, &graph, node, parent)
}

pub fn render_clade(tree: &Tree, node: usize, parent: Option<usize>) -> String {
    let graph = adjacency(tree);
    fn render(tree: &Tree, graph: &[Vec<usize>], node: usize, parent: Option<usize>) -> String {
        let children: Vec<_> = graph[node]
            .iter()
            .copied()
            .filter(|next| Some(*next) != parent)
            .collect();
        if children.is_empty() {
            return tree.nodes[node]
                .name
                .as_deref()
                .map(newick_label)
                .unwrap_or_else(|| newick_label(""));
        }
        format!(
            "({})",
            children
                .into_iter()
                .map(|next| render(tree, graph, next, Some(node)))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
    render(tree, &graph, node, parent)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parses_iqtree_style_branches_and_selects_unique_clade() {
        let tree = parse_newick("((A|OG|candidate_1:0.1,B|OG|candidate_1:0.2)95:0.3,(C|OG|candidate_1:0.1,D|OG|candidate_1:0.2):0.3);").unwrap();
        let clades = select_scogs(&tree, 4, &BTreeSet::new()).unwrap();
        assert_eq!(clades.len(), 1);
        assert_eq!(clades[0].samples.len(), 4);
    }
    #[test]
    fn reroots_an_unrooted_iqtree_shape_at_monophyletic_outgroup_edge() {
        let tree = parse_newick("((A|OG|candidate_1,B|OG|candidate_1),(C|OG|candidate_1,(O|OG|candidate_1,O2|OG|candidate_1)),D|OG|candidate_1);").unwrap();
        let outgroups = ["O".to_string(), "O2".to_string()].into_iter().collect();
        let clades = select_scogs(&tree, 4, &outgroups).unwrap();
        assert_eq!(clades.len(), 1);
        assert_eq!(
            clades[0].samples,
            [
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string()
            ]
            .into_iter()
            .collect()
        );
    }
    #[test]
    fn escapes_newick_sensitive_sample_labels() {
        let tree = parse_newick("('A sample|OG1|candidate_1','B,C|OG1|candidate_1');").unwrap();
        assert_eq!(
            render_clade_samples(&tree, tree.root, None),
            "('A sample','B,C')"
        );
    }

    #[test]
    fn renders_strict_clade_with_sample_labels() {
        let tree = parse_newick("(A|OG1|candidate_1,B|OG1|candidate_2);").unwrap();
        assert_eq!(render_clade_samples(&tree, tree.root, None), "(A,B)");
    }

    #[test]
    fn finds_one_to_one_complement_across_arbitrary_newick_root() {
        let tree = parse_newick(
            "((A|OG|candidate_1,A|OG|candidate_2),B|OG|candidate_1,C|OG|candidate_1,D|OG|candidate_1);",
        )
        .unwrap();
        let clades = select_scogs(&tree, 3, &BTreeSet::new()).unwrap();
        assert!(clades.iter().any(|clade| {
            clade.samples
                == [
                    "A".to_string(),
                    "B".to_string(),
                    "C".to_string(),
                    "D".to_string(),
                ]
                .into_iter()
                .collect()
        }));
    }

    #[test]
    fn rejects_nonmonophyletic_outgroups() {
        let tree = parse_newick(
            "(A|OG|candidate_1,(O|OG|candidate_1,B|OG|candidate_1),O2|OG|candidate_1);",
        )
        .unwrap();
        let outgroups = ["O".to_string(), "O2".to_string()].into_iter().collect();
        assert_eq!(
            select_scogs(&tree, 2, &outgroups).unwrap_err(),
            "outgroup_missing_or_not_monophyletic"
        );
    }
}
