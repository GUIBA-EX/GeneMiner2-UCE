//! PanRefV2: a bounded, two-pass UCE local-graph reference builder.
//!
//! It deliberately never writes recruited FASTQ.  The first scan accepts only
//! unique paired bait hits; the second scan admits a single-mate rescue only
//! when it overlaps the already supported local graph.
use super::backbone::{assemble_backbone_from_unitigs, PathEvidence};
use super::bait::BaitCatalog;
use super::bait_index::BaitIndex;
use super::dbg::{unitig_edges, KmerCounter, Unitig};
use super::recruit::stream_recruited_pairs;
use crate::{
    safe_locus_name, write_reference_manifest, write_sample_manifest, AppResult, Args, Sample,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};

const K: usize = 31;
const MIN_KMER_COUNT: u32 = 2;
const MAX_CORE_PAIRS_PER_SAMPLE_LOCUS: usize = 4_000;
const MAX_RESCUE_PAIRS_PER_SAMPLE_LOCUS: usize = 1_000;
const UNITIG_SEED_LEN: usize = 21;
const AMBIGUOUS_UNITIG: usize = usize::MAX;

struct LocusAccumulator {
    counter: KmerCounter,
    core_pairs: u64,
    rescued_pairs: u64,
    sample_support: u32,
}

impl LocusAccumulator {
    fn new() -> Self {
        Self {
            counter: KmerCounter::new(K).expect("fixed valid k-mer size"),
            core_pairs: 0,
            rescued_pairs: 0,
            sample_support: 0,
        }
    }

    fn add_pair(&mut self, first: &[u8], second: Option<&[u8]>) {
        self.counter.add_read(first);
        if let Some(second) = second {
            self.counter.add_read(second);
        }
    }
}

#[derive(Default)]
struct SampleSummary {
    name: String,
    strong: u64,
    candidate_rescues: u64,
    accepted_rescues: u64,
    ambiguous: u64,
    loci_with_core: usize,
}

fn base_bits(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn build_unitig_seed_index(unitigs: &[Unitig]) -> HashMap<u64, usize> {
    let mut index = HashMap::new();
    let mask = (1_u64 << (2 * UNITIG_SEED_LEN)) - 1;
    for (unitig_id, unitig) in unitigs.iter().enumerate() {
        let mut code = 0_u64;
        let mut valid = 0_usize;
        for &base in &unitig.sequence {
            let Some(bits) = base_bits(base) else {
                code = 0;
                valid = 0;
                continue;
            };
            code = ((code << 2) | bits) & mask;
            valid += 1;
            if valid >= UNITIG_SEED_LEN {
                match index.get_mut(&code) {
                    None => {
                        index.insert(code, unitig_id);
                    }
                    Some(existing) if *existing != unitig_id => *existing = AMBIGUOUS_UNITIG,
                    Some(_) => {}
                }
            }
        }
    }
    index
}

fn read_unitig_path(read: &[u8], index: &HashMap<u64, usize>) -> Vec<usize> {
    if read.len() < UNITIG_SEED_LEN {
        return Vec::new();
    }
    let mask = (1_u64 << (2 * UNITIG_SEED_LEN)) - 1;
    let mut code = 0_u64;
    let mut valid = 0_usize;
    let mut path = Vec::new();
    for &base in read {
        let Some(bits) = base_bits(base) else {
            code = 0;
            valid = 0;
            continue;
        };
        code = ((code << 2) | bits) & mask;
        valid += 1;
        if valid >= UNITIG_SEED_LEN {
            let Some(&unitig) = index.get(&code) else {
                continue;
            };
            if unitig != AMBIGUOUS_UNITIG && path.last().copied() != Some(unitig) {
                path.push(unitig);
            }
        }
    }
    path
}

fn unique_unitig(path: &[usize]) -> Option<usize> {
    if path.len() == 1 {
        Some(path[0])
    } else {
        None
    }
}

fn stable_depth_score(depths: &[u32]) -> u64 {
    let mut nonzero = depths
        .iter()
        .copied()
        .filter(|&depth| depth > 0)
        .collect::<Vec<_>>();
    if nonzero.is_empty() {
        return 0;
    }
    nonzero.sort_unstable();
    let median = nonzero[nonzero.len() / 2];
    let mut deviations = nonzero
        .iter()
        .map(|depth| depth.abs_diff(median))
        .collect::<Vec<_>>();
    deviations.sort_unstable();
    let mad = deviations[deviations.len() / 2];
    u64::from(median) * 100 / (u64::from(mad) + 1)
}

fn write_gfa_record(
    out: &mut BufWriter<File>,
    locus_id: usize,
    unitigs: &[Unitig],
) -> AppResult<()> {
    for (id, unitig) in unitigs.iter().enumerate() {
        writeln!(
            out,
            "S\tL{locus_id}_U{id}\t{}\tKC:i:{}",
            String::from_utf8_lossy(&unitig.sequence),
            unitig.kmer_count
        )
        .map_err(|e| e.to_string())?;
    }
    for (from, to) in unitig_edges(unitigs, K) {
        writeln!(
            out,
            "L\tL{locus_id}_U{from}\t+\tL{locus_id}_U{to}\t+\t{}M",
            K - 1
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Build an independent PanRefV2 reference and its graph/QC artifacts.
pub(crate) fn build_reference(args: &Args, samples: &[Sample]) -> AppResult<std::path::PathBuf> {
    let baits = args
        .panref_baits
        .as_ref()
        .ok_or("--engine panrefv2 requires --panref-baits")?;
    let reference_dir = args.output.join("population").join("reference");
    let output_dir = reference_dir.join("panrefv2");
    fs::create_dir_all(&output_dir)
        .map_err(|e| format!("unable to create {}: {e}", output_dir.display()))?;
    let catalog = BaitCatalog::read(baits)?;
    let index = BaitIndex::build_catalog(&catalog)?;
    let loci = catalog
        .loci
        .iter()
        .map(|x| x.name.clone())
        .collect::<Vec<_>>();
    index.write_metadata(&output_dir.join("index_metadata.tsv"), &loci)?;
    write_sample_manifest(
        &args.output.join("population").join("sample_manifest.tsv"),
        samples,
    )?;

    let mut accumulators = (0..catalog.loci.len())
        .map(|_| LocusAccumulator::new())
        .collect::<Vec<_>>();
    let mut summaries = Vec::with_capacity(samples.len());

    // Pass 1: only a unique locus hit in both mates enters the core graph.
    for sample in samples {
        let mut per_locus = vec![0_usize; catalog.loci.len()];
        let mut summary = SampleSummary {
            name: sample.internal.clone(),
            ..SampleSummary::default()
        };
        let stats = stream_recruited_pairs(
            &index,
            &catalog,
            &sample.read1,
            &sample.read2,
            args.threads,
            false,
            |locus, strong, first, second| {
                if !strong {
                    return;
                }
                let id = locus as usize;
                if per_locus[id] == 0 {
                    accumulators[id].sample_support += 1;
                }
                if per_locus[id] < MAX_CORE_PAIRS_PER_SAMPLE_LOCUS {
                    accumulators[id].add_pair(first, second);
                    accumulators[id].core_pairs += 1;
                    per_locus[id] += 1;
                }
            },
        )?;
        summary.strong = stats.strong_pairs;
        summary.candidate_rescues = stats.rescued_pairs;
        summary.ambiguous = stats.ambiguous_pairs;
        summary.loci_with_core = per_locus.iter().filter(|&&n| n > 0).count();
        summaries.push(summary);
    }
    // Freeze the strict paired core. Rescue decisions must never become
    // self-reinforcing as new rescue reads are appended to the graph.
    let core_gates = accumulators
        .iter()
        .map(|entry| entry.counter.clone())
        .collect::<Vec<_>>();

    // Pass 2: a one-mate pair is rescued only when it extends the existing
    // supported graph. This is deliberately bounded per sample and locus.
    for (sample, summary) in samples.iter().zip(summaries.iter_mut()) {
        let mut per_locus = vec![0_usize; catalog.loci.len()];
        stream_recruited_pairs(
            &index,
            &catalog,
            &sample.read1,
            &sample.read2,
            args.threads,
            true,
            |locus, strong, first, second| {
                if strong {
                    return;
                }
                let id = locus as usize;
                if per_locus[id] >= MAX_RESCUE_PAIRS_PER_SAMPLE_LOCUS {
                    return;
                }
                let touches_core = core_gates[id].has_supported_kmer(first, MIN_KMER_COUNT)
                    || second.is_some_and(|read| {
                        core_gates[id].has_supported_kmer(read, MIN_KMER_COUNT)
                    });
                if touches_core {
                    accumulators[id].add_pair(first, second);
                    accumulators[id].rescued_pairs += 1;
                    per_locus[id] += 1;
                    summary.accepted_rescues += 1;
                }
            },
        )?;
    }

    let unitigs = accumulators
        .iter_mut()
        .map(|entry| {
            std::mem::replace(
                &mut entry.counter,
                KmerCounter::new(K).expect("fixed valid k-mer size"),
            )
            .into_unitigs(MIN_KMER_COUNT)
        })
        .collect::<Vec<_>>();

    // Pass 3: replay the exact bounded acceptance policy, then project
    // accepted pairs onto final unitigs for sample and edge evidence.
    let unitig_seed_indexes = unitigs
        .iter()
        .map(|nodes| build_unitig_seed_index(nodes))
        .collect::<Vec<_>>();
    let graph_edges = unitigs
        .iter()
        .map(|nodes| unitig_edges(nodes, K).into_iter().collect::<HashSet<_>>())
        .collect::<Vec<_>>();
    let mut pe_support = unitigs
        .iter()
        .map(|nodes| vec![0_u64; nodes.len()])
        .collect::<Vec<_>>();
    let mut edge_support = unitigs
        .iter()
        .map(|_| HashMap::new())
        .collect::<Vec<HashMap<(usize, usize), u64>>>();
    let mut unitig_sample_depth = unitigs
        .iter()
        .map(|nodes| {
            nodes
                .iter()
                .map(|_| vec![0_u32; samples.len()])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    for (sample_index, sample) in samples.iter().enumerate() {
        let mut local_depth = unitigs
            .iter()
            .map(|nodes| vec![0_u32; nodes.len()])
            .collect::<Vec<_>>();
        let mut core_pairs = vec![0_usize; catalog.loci.len()];
        let mut rescue_pairs = vec![0_usize; catalog.loci.len()];
        stream_recruited_pairs(
            &index,
            &catalog,
            &sample.read1,
            &sample.read2,
            args.threads,
            true,
            |locus, strong, first, second| {
                let id = locus as usize;
                if strong {
                    if core_pairs[id] >= MAX_CORE_PAIRS_PER_SAMPLE_LOCUS {
                        return;
                    }
                    core_pairs[id] += 1;
                } else {
                    if rescue_pairs[id] >= MAX_RESCUE_PAIRS_PER_SAMPLE_LOCUS {
                        return;
                    }
                    let touches_core = core_gates[id].has_supported_kmer(first, MIN_KMER_COUNT)
                        || second.is_some_and(|read| {
                            core_gates[id].has_supported_kmer(read, MIN_KMER_COUNT)
                        });
                    if !touches_core {
                        return;
                    }
                    rescue_pairs[id] += 1;
                }
                let left_path = read_unitig_path(first, &unitig_seed_indexes[id]);
                let right_path = second
                    .map(|read| read_unitig_path(read, &unitig_seed_indexes[id]))
                    .unwrap_or_default();
                for &node in left_path.iter().chain(right_path.iter()) {
                    local_depth[id][node] += 1;
                }
                for path in [&left_path, &right_path] {
                    for edge in path.windows(2).map(|nodes| (nodes[0], nodes[1])) {
                        if graph_edges[id].contains(&edge) {
                            *edge_support[id].entry(edge).or_default() += 1;
                        }
                    }
                }
                let left = unique_unitig(&left_path);
                let right = unique_unitig(&right_path);
                if let (Some(left), Some(right)) = (left, right) {
                    if left != right {
                        let forward = (left, right);
                        let reverse = (right, left);
                        if graph_edges[id].contains(&forward) {
                            pe_support[id][left] += 1;
                            pe_support[id][right] += 1;
                            *edge_support[id].entry(forward).or_default() += 1;
                        }
                        if graph_edges[id].contains(&reverse) {
                            pe_support[id][left] += 1;
                            pe_support[id][right] += 1;
                            *edge_support[id].entry(reverse).or_default() += 1;
                        }
                    }
                }
            },
        )?;
        for (locus, depths) in local_depth.into_iter().enumerate() {
            for (unitig, depth) in depths.into_iter().enumerate() {
                unitig_sample_depth[locus][unitig][sample_index] = depth;
            }
        }
    }
    let unitig_sample_support = unitig_sample_depth
        .iter()
        .map(|nodes| {
            nodes
                .iter()
                .map(|depths| depths.iter().filter(|&&depth| depth > 0).count() as u32)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let unitig_depth_stability = unitig_sample_depth
        .iter()
        .map(|nodes| {
            nodes
                .iter()
                .map(|depths| stable_depth_score(depths))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut recruitment = BufWriter::new(
        File::create(output_dir.join("recruitment_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(recruitment, "sample\tstrong_pairs\tcandidate_single_mate_pairs\taccepted_graph_rescues\tambiguous_pairs\tloci_with_core")
        .map_err(|e| e.to_string())?;
    for row in &summaries {
        writeln!(
            recruitment,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.name,
            row.strong,
            row.candidate_rescues,
            row.accepted_rescues,
            row.ambiguous,
            row.loci_with_core
        )
        .map_err(|e| e.to_string())?;
    }
    recruitment.flush().map_err(|e| e.to_string())?;

    let fasta_path = reference_dir.join("population_reference.fasta");
    let mut fasta = BufWriter::new(File::create(&fasta_path).map_err(|e| e.to_string())?);
    let mut graph = BufWriter::new(
        File::create(output_dir.join("population_graph.gfa")).map_err(|e| e.to_string())?,
    );
    writeln!(graph, "H\tVN:Z:1.0\tTS:Z:GeneMiner2-PanRefV2").map_err(|e| e.to_string())?;
    let mut report = BufWriter::new(
        File::create(output_dir.join("locus_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(report, "locus\tstatus\tsequence_length\tunitigs\tcore_pairs\taccepted_rescues\tsupporting_samples\tpe_links\tbackbone_min_samples\tbackbone_edge_support\tbackbone_min_depth_stability")
        .map_err(|e| e.to_string())?;
    let mut names = HashSet::new();
    let mut written = 0_usize;
    for (id, locus) in catalog.loci.iter().enumerate() {
        write_gfa_record(&mut graph, id, &unitigs[id])?;
        let pe_links: u64 = pe_support[id].iter().sum::<u64>() / 2;
        let evidence = PathEvidence {
            sample_support: &unitig_sample_support[id],
            pe_support: &pe_support[id],
            depth_stability: &unitig_depth_stability[id],
            edge_support: &edge_support[id],
        };
        let backbone = assemble_backbone_from_unitigs(&unitigs[id], &locus.records, evidence);
        let backbone_min_samples = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .iter()
                .map(|&node| unitig_sample_support[id][node])
                .min()
                .unwrap_or_default()
        });
        let backbone_edge_support = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .windows(2)
                .map(|edge| {
                    edge_support[id]
                        .get(&(edge[0], edge[1]))
                        .copied()
                        .unwrap_or_default()
                })
                .sum()
        });
        let backbone_min_depth_stability = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .iter()
                .map(|&node| unitig_depth_stability[id][node])
                .min()
                .unwrap_or_default()
        });
        let sequence = backbone.as_ref().map(|path| &path.sequence);
        let (status, length) = match sequence.as_ref() {
            None if accumulators[id].core_pairs == 0 => ("no_core", 0),
            None if unitigs[id].is_empty() => ("low_coverage", 0),
            None => ("complex", 0),
            Some(seq) if seq.len() < 100 => ("short", seq.len()),
            Some(seq) if accumulators[id].sample_support < 2 => ("low_sample_support", seq.len()),
            Some(seq) => ("pass", seq.len()),
        };
        if status == "pass" || args.panrefv2_include_low_confidence {
            let Some(sequence) = sequence else { continue };
            let name = safe_locus_name(&locus.name);
            if !names.insert(name.clone()) {
                return Err(format!(
                    "PanRefV2 locus names collide after FASTA sanitization: {}",
                    locus.name
                ));
            }
            writeln!(fasta, ">{name}").map_err(|e| e.to_string())?;
            for line in sequence.chunks(80) {
                writeln!(fasta, "{}", String::from_utf8_lossy(line)).map_err(|e| e.to_string())?;
            }
            written += 1;
        }
        writeln!(
            report,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            locus.name,
            status,
            length,
            unitigs[id].len(),
            accumulators[id].core_pairs,
            accumulators[id].rescued_pairs,
            accumulators[id].sample_support,
            pe_links,
            backbone_min_samples,
            backbone_edge_support,
            backbone_min_depth_stability
        )
        .map_err(|e| e.to_string())?;
    }
    fasta.flush().map_err(|e| e.to_string())?;
    graph.flush().map_err(|e| e.to_string())?;
    report.flush().map_err(|e| e.to_string())?;
    if written == 0 {
        return Err("PanRefV2 found no bait-anchored local backbones; inspect population/reference/panrefv2".into());
    }
    write_reference_manifest(&reference_dir, "panrefv2", &fasta_path)?;
    Ok(fasta_path)
}

#[cfg(test)]
mod tests {
    use super::{build_unitig_seed_index, read_unitig_path, unique_unitig};
    use crate::panref::dbg::Unitig;

    #[test]
    fn only_unique_unitig_seeds_are_used_for_pe_links() {
        let unitigs = vec![
            Unitig {
                sequence: b"AACCGGTTAACCGGTTAACCGGTT".to_vec(),
                kmer_count: 4,
            },
            Unitig {
                sequence: b"TTTTCCCCAAAAGGGGTTTTCCCC".to_vec(),
                kmer_count: 4,
            },
        ];
        let index = build_unitig_seed_index(&unitigs);
        assert_eq!(
            unique_unitig(&read_unitig_path(b"AACCGGTTAACCGGTTAACCGGTT", &index)),
            Some(0)
        );
        assert_eq!(
            unique_unitig(&read_unitig_path(b"ACGTACGTACGTACGTACGTACGT", &index)),
            None
        );
        let edge_unitigs = vec![
            Unitig {
                sequence: b"AAAAAAAAAAAAAAAAAAAAA".to_vec(),
                kmer_count: 1,
            },
            Unitig {
                sequence: b"CCCCCCCCCCCCCCCCCCCCC".to_vec(),
                kmer_count: 1,
            },
        ];
        let edge_index = build_unitig_seed_index(&edge_unitigs);
        assert_eq!(
            read_unitig_path(b"AAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCC", &edge_index),
            vec![0, 1]
        );
    }
}
