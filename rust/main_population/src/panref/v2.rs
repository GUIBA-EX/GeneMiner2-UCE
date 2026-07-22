//! PanRefV2.2: a single-scan UCE local-graph reference builder.
//!
//! It deliberately never writes recruited FASTQ. Strict paired hits stay
//! memory-bounded; one-mate candidates are streamed to an ephemeral spool and
//! replayed only after the supported core graph has been frozen.
use super::backbone::{assemble_backbone_from_unitigs, PathEvidence, ResolvedBackbone};
use super::bait::BaitCatalog;
use super::bait_index::BaitIndex;
use super::bubble::inspect_backbone_branches;
use super::dbg::{unitig_edges, KmerCounter, Unitig};
use super::recruit::stream_recruited_pairs;
use crate::{
    safe_locus_name, write_reference_manifest, write_sample_manifest, AppResult, Args, Sample,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};

const K: usize = 31;
const ADAPTIVE_K: usize = 25;
const MIN_KMER_COUNT: u32 = 2;
const MAX_CORE_PAIRS_PER_SAMPLE_LOCUS: usize = 4_000;
const MAX_RESCUE_PAIRS_PER_SAMPLE_LOCUS: usize = 1_000;
const UNITIG_SEED_LEN: usize = 21;
const AMBIGUOUS_UNITIG: usize = usize::MAX;

#[derive(Clone)]
struct LedgerPair {
    first: Vec<u8>,
    second: Option<Vec<u8>>,
}

impl LedgerPair {
    fn from_records(
        first: &super::recruit::FastqRecord,
        second: Option<&super::recruit::FastqRecord>,
    ) -> Self {
        Self {
            first: quality_masked_sequence(first),
            second: second.map(quality_masked_sequence),
        }
    }

    fn fingerprint(&self) -> u64 {
        let mut hash = 0xcbf2_9ce4_8422_2325_u64;
        for base in &self.first {
            hash ^= u64::from(*base);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        hash ^= 0xff; // preserve the mate boundary in the pair key.
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        if let Some(second) = &self.second {
            for base in second {
                hash ^= u64::from(*base);
                hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        hash
    }
}

#[derive(Default)]
struct LedgerBucket {
    core: Vec<LedgerPair>,
    accepted_rescues: Vec<LedgerPair>,
    core_seen: HashSet<u64>,
    accepted_seen: HashSet<u64>,
}

struct SampleLedger {
    buckets: Vec<Option<LedgerBucket>>,
    candidate_spool: std::path::PathBuf,
}

struct CandidateSpool {
    path: std::path::PathBuf,
}

impl CandidateSpool {
    fn create(output_dir: &std::path::Path) -> AppResult<Self> {
        let path = output_dir.join(format!(".candidate-spool-{}", std::process::id()));
        fs::create_dir_all(&path)
            .map_err(|e| format!("unable to create {}: {e}", path.display()))?;
        Ok(Self { path })
    }
}

impl Drop for CandidateSpool {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

impl SampleLedger {
    fn new(locus_count: usize, candidate_spool: std::path::PathBuf) -> Self {
        Self {
            buckets: (0..locus_count).map(|_| None).collect(),
            candidate_spool,
        }
    }

    fn bucket_mut(&mut self, locus: usize) -> &mut LedgerBucket {
        self.buckets[locus].get_or_insert_with(LedgerBucket::default)
    }
}

fn quality_masked_sequence(record: &super::recruit::FastqRecord) -> Vec<u8> {
    record
        .sequence
        .iter()
        .zip(&record.quality)
        .map(|(base, quality)| {
            if quality.saturating_sub(33) >= 20 {
                base.to_ascii_uppercase()
            } else {
                b'N'
            }
        })
        .collect()
}

fn write_candidate_record(
    out: &mut BufWriter<File>,
    locus: usize,
    pair: &LedgerPair,
) -> AppResult<()> {
    let locus = u32::try_from(locus).map_err(|_| "too many PanRefV2 loci".to_string())?;
    let first_len = u32::try_from(pair.first.len())
        .map_err(|_| "candidate read exceeds u32 length".to_string())?;
    let second_len = pair
        .second
        .as_ref()
        .map(|read| {
            u32::try_from(read.len()).map_err(|_| "candidate mate exceeds u32 length".to_string())
        })
        .transpose()?
        .unwrap_or(u32::MAX);
    out.write_all(&locus.to_le_bytes())
        .map_err(|e| e.to_string())?;
    out.write_all(&first_len.to_le_bytes())
        .map_err(|e| e.to_string())?;
    out.write_all(&second_len.to_le_bytes())
        .map_err(|e| e.to_string())?;
    out.write_all(&pair.first).map_err(|e| e.to_string())?;
    if let Some(second) = &pair.second {
        out.write_all(second).map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn read_candidate_record(
    input: &mut BufReader<File>,
    locus_count: usize,
) -> AppResult<Option<(usize, LedgerPair)>> {
    let mut header = [0_u8; 12];
    match input.read(&mut header[..1]) {
        Ok(0) => return Ok(None),
        Ok(1) => input
            .read_exact(&mut header[1..])
            .map_err(|e| format!("truncated PanRefV2 candidate spool: {e}"))?,
        Ok(_) => unreachable!("single-byte read buffer"),
        Err(error) => return Err(error.to_string()),
    }
    let locus = usize::try_from(u32::from_le_bytes(
        header[0..4].try_into().expect("slice length"),
    ))
    .expect("u32 fits usize");
    let first_len = usize::try_from(u32::from_le_bytes(
        header[4..8].try_into().expect("slice length"),
    ))
    .expect("u32 fits usize");
    let second_raw = u32::from_le_bytes(header[8..12].try_into().expect("slice length"));
    if locus >= locus_count || first_len > 1_000_000 {
        return Err("invalid PanRefV2 candidate spool record".into());
    }
    let second_len =
        (second_raw != u32::MAX).then_some(usize::try_from(second_raw).expect("u32 fits usize"));
    if second_len.is_some_and(|length| length > 1_000_000) {
        return Err("invalid PanRefV2 candidate mate length".into());
    }
    let mut first = vec![0_u8; first_len];
    input
        .read_exact(&mut first)
        .map_err(|e| format!("truncated PanRefV2 candidate sequence: {e}"))?;
    let second = if let Some(length) = second_len {
        let mut read = vec![0_u8; length];
        input
            .read_exact(&mut read)
            .map_err(|e| format!("truncated PanRefV2 candidate mate: {e}"))?;
        Some(read)
    } else {
        None
    };
    Ok(Some((locus, LedgerPair { first, second })))
}

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
    candidate_spool_bytes: u64,
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

type ActiveEvidence<'a> = (&'a [u32], &'a [u64], &'a HashMap<(usize, usize), u64>);

struct LocusEvidence {
    sample_support: Vec<u32>,
    pe_support: Vec<u64>,
    depth_stability: Vec<u64>,
    edge_support: HashMap<(usize, usize), u64>,
    sample_depth: Vec<Vec<(usize, u32)>>,
}

fn project_locus_evidence(
    ledgers: &[SampleLedger],
    locus: usize,
    unitigs: &[Unitig],
    edges: &[(usize, usize)],
) -> LocusEvidence {
    let seed_index = build_unitig_seed_index(unitigs);
    let edge_set = edges.iter().copied().collect::<HashSet<_>>();
    let mut pe_support = vec![0_u64; unitigs.len()];
    let mut edge_support = HashMap::new();
    let mut sample_depth = unitigs
        .iter()
        .map(|_| Vec::new())
        .collect::<Vec<Vec<(usize, u32)>>>();
    for (sample_id, ledger) in ledgers.iter().enumerate() {
        let Some(bucket) = &ledger.buckets[locus] else {
            continue;
        };
        let mut local_depth = vec![0_u32; unitigs.len()];
        for pair in bucket.core.iter().chain(&bucket.accepted_rescues) {
            let left_path = read_unitig_path(&pair.first, &seed_index);
            let right_path = pair
                .second
                .as_deref()
                .map(|read| read_unitig_path(read, &seed_index))
                .unwrap_or_default();
            for &node in left_path.iter().chain(&right_path) {
                local_depth[node] += 1;
            }
            for path in [&left_path, &right_path] {
                for edge in path.windows(2).map(|nodes| (nodes[0], nodes[1])) {
                    if edge_set.contains(&edge) {
                        *edge_support.entry(edge).or_default() += 1;
                    }
                }
            }
            if let (Some(left), Some(right)) =
                (unique_unitig(&left_path), unique_unitig(&right_path))
            {
                if left != right {
                    for edge in [(left, right), (right, left)] {
                        if edge_set.contains(&edge) {
                            pe_support[left] += 1;
                            pe_support[right] += 1;
                            *edge_support.entry(edge).or_default() += 1;
                        }
                    }
                }
            }
        }
        for (unitig, depth) in local_depth.into_iter().enumerate() {
            if depth > 0 {
                sample_depth[unitig].push((sample_id, depth));
            }
        }
    }
    let sample_support = sample_depth
        .iter()
        .map(|depths| depths.len() as u32)
        .collect();
    let depth_stability = sample_depth
        .iter()
        .map(|depths| {
            stable_depth_score(&depths.iter().map(|(_, depth)| *depth).collect::<Vec<_>>())
        })
        .collect();
    LocusEvidence {
        sample_support,
        pe_support,
        depth_stability,
        edge_support,
        sample_depth,
    }
}

fn should_retry_adaptive(backbone: &Option<ResolvedBackbone>, core_pairs: u64) -> bool {
    backbone.is_none() && core_pairs > 0
}

fn rebuild_locus_unitigs(ledgers: &[SampleLedger], locus: usize, k: usize) -> Vec<Unitig> {
    let mut counter = KmerCounter::new(k).expect("fixed valid adaptive k-mer size");
    for ledger in ledgers {
        let Some(bucket) = &ledger.buckets[locus] else {
            continue;
        };
        for pair in bucket.core.iter().chain(&bucket.accepted_rescues) {
            counter.add_read(&pair.first);
            if let Some(second) = pair.second.as_deref() {
                counter.add_read(second);
            }
        }
    }
    counter.into_unitigs(MIN_KMER_COUNT)
}

fn gfa_path_segments(graph_name: &str, path: &ResolvedBackbone) -> String {
    if path.reversed {
        path.nodes
            .iter()
            .rev()
            .map(|node| format!("{graph_name}_U{node}-"))
            .collect::<Vec<_>>()
            .join(",")
    } else {
        path.nodes
            .iter()
            .map(|node| format!("{graph_name}_U{node}+"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn write_sample_backbone_paths(
    out: &mut BufWriter<File>,
    ledgers: &[SampleLedger],
    locus: usize,
    unitigs: &[Unitig],
    backbone: &ResolvedBackbone,
    graph_name: &str,
    summaries: &[SampleSummary],
) -> AppResult<()> {
    let index = build_unitig_seed_index(unitigs);
    let path_nodes = &backbone.nodes;
    let required_nodes = path_nodes.iter().copied().collect::<HashSet<_>>();
    let required_edges = path_nodes
        .windows(2)
        .map(|nodes| (nodes[0], nodes[1]))
        .collect::<HashSet<_>>();
    let emitted_path = gfa_path_segments(graph_name, backbone);
    for (sample_id, ledger) in ledgers.iter().enumerate() {
        let mut covered_nodes = HashSet::new();
        let mut supported_edges = HashSet::new();
        if let Some(bucket) = &ledger.buckets[locus] {
            for pair in bucket.core.iter().chain(&bucket.accepted_rescues) {
                let left_path = read_unitig_path(&pair.first, &index);
                let right_path = pair
                    .second
                    .as_deref()
                    .map(|read| read_unitig_path(read, &index))
                    .unwrap_or_default();
                for node in left_path.iter().chain(&right_path) {
                    if required_nodes.contains(node) {
                        covered_nodes.insert(*node);
                    }
                }
                for path in [&left_path, &right_path] {
                    for edge in path.windows(2).map(|nodes| (nodes[0], nodes[1])) {
                        if required_edges.contains(&edge) {
                            supported_edges.insert(edge);
                        }
                    }
                }
                if let (Some(left), Some(right)) =
                    (unique_unitig(&left_path), unique_unitig(&right_path))
                {
                    for edge in [(left, right), (right, left)] {
                        if required_edges.contains(&edge) {
                            supported_edges.insert(edge);
                        }
                    }
                }
            }
        }
        let complete = covered_nodes.len() == required_nodes.len()
            && supported_edges.len() == required_edges.len();
        let status = if complete {
            "backbone_path"
        } else if covered_nodes.is_empty() {
            "no_coverage"
        } else {
            "partial"
        };
        let path = if complete { emitted_path.as_str() } else { "" };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            summaries[sample_id].name,
            graph_name,
            status,
            covered_nodes.len(),
            required_nodes.len(),
            supported_edges.len(),
            path
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

struct BubbleQcView<'a> {
    locus: &'a str,
    graph_name: &'a str,
    unitigs: &'a [Unitig],
    edges: &'a [(usize, usize)],
    backbone: &'a ResolvedBackbone,
    sample_support: &'a [u32],
    edge_support: &'a HashMap<(usize, usize), u64>,
}

fn write_bubble_qc(out: &mut BufWriter<File>, view: BubbleQcView<'_>) -> AppResult<()> {
    let BubbleQcView {
        locus,
        graph_name,
        unitigs,
        edges,
        backbone,
        sample_support,
        edge_support,
    } = view;
    for record in inspect_backbone_branches(unitigs, edges, backbone, sample_support, edge_support)
    {
        let exit = record
            .exit
            .map(|node| format!("{graph_name}_U{node}"))
            .unwrap_or_default();
        writeln!(
            out,
            "{locus}\t{graph_name}\t{graph_name}_U{}\t{exit}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.entry,
            record.status,
            record.alternative_branches,
            record.alternative_nodes,
            record.alternative_min_samples,
            record.alternative_min_edge_support,
            record.canonical_min_samples,
            record.canonical_min_edge_support,
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(())
}

#[derive(Debug, Eq, PartialEq)]
struct BackboneCoordinate {
    node: usize,
    orientation: char,
    start: usize,
    end: usize,
}

const SHA256_INITIAL: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];
const SHA256_ROUND: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

fn sha256_bytes(bytes: &[u8]) -> [u8; 32] {
    let mut message = bytes.to_vec();
    let bit_length = (message.len() as u64).wrapping_mul(8);
    message.push(0x80);
    while !(message.len() + 8).is_multiple_of(64) {
        message.push(0);
    }
    message.extend_from_slice(&bit_length.to_be_bytes());
    let mut state = SHA256_INITIAL;
    for chunk in message.chunks_exact(64) {
        let mut words = [0_u32; 64];
        for (index, word) in words[..16].iter_mut().enumerate() {
            *word = u32::from_be_bytes(chunk[index * 4..index * 4 + 4].try_into().expect("word"));
        }
        for index in 16..64 {
            let sigma0 = words[index - 15].rotate_right(7)
                ^ words[index - 15].rotate_right(18)
                ^ (words[index - 15] >> 3);
            let sigma1 = words[index - 2].rotate_right(17)
                ^ words[index - 2].rotate_right(19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16]
                .wrapping_add(sigma0)
                .wrapping_add(words[index - 7])
                .wrapping_add(sigma1);
        }
        let mut working = state;
        for (index, &constant) in SHA256_ROUND.iter().enumerate() {
            let choice = (working[4] & working[5]) ^ ((!working[4]) & working[6]);
            let major =
                (working[0] & working[1]) ^ (working[0] & working[2]) ^ (working[1] & working[2]);
            let sum1 = working[4].rotate_right(6)
                ^ working[4].rotate_right(11)
                ^ working[4].rotate_right(25);
            let sum0 = working[0].rotate_right(2)
                ^ working[0].rotate_right(13)
                ^ working[0].rotate_right(22);
            let temporary1 = working[7]
                .wrapping_add(sum1)
                .wrapping_add(choice)
                .wrapping_add(constant)
                .wrapping_add(words[index]);
            let temporary2 = sum0.wrapping_add(major);
            working = [
                temporary1.wrapping_add(temporary2),
                working[0],
                working[1],
                working[2],
                working[3].wrapping_add(temporary1),
                working[4],
                working[5],
                working[6],
            ];
        }
        for (target, value) in state.iter_mut().zip(working) {
            *target = target.wrapping_add(value);
        }
    }
    let mut output = [0_u8; 32];
    for (index, value) in state.iter().enumerate() {
        output[index * 4..index * 4 + 4].copy_from_slice(&value.to_be_bytes());
    }
    output
}

fn sha256_hex(bytes: &[u8]) -> String {
    sha256_bytes(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn stable_backbone_id(locus: &str, sequence: &[u8]) -> String {
    let mut identity = b"GeneMiner2-PanRefV2/backbone/v1\0".to_vec();
    identity.extend_from_slice(locus.as_bytes());
    identity.push(0);
    identity.extend_from_slice(sequence);
    format!("panrefv2-sha256-{}", sha256_hex(&identity))
}

fn backbone_coordinates(
    unitigs: &[Unitig],
    k: usize,
    backbone: &ResolvedBackbone,
) -> AppResult<Vec<BackboneCoordinate>> {
    if k < 2 {
        return Err("PanRefV2 backbone coordinate k must be at least two".into());
    }
    let nodes = if backbone.reversed {
        backbone.nodes.iter().rev().copied().collect::<Vec<_>>()
    } else {
        backbone.nodes.clone()
    };
    let orientation = if backbone.reversed {
        char::from(45_u8)
    } else {
        char::from(43_u8)
    };
    let mut coordinates = Vec::with_capacity(nodes.len());
    let mut offset = 0_usize;
    for (rank, node) in nodes.into_iter().enumerate() {
        let unitig = unitigs
            .get(node)
            .ok_or_else(|| "PanRefV2 backbone references an absent unitig".to_string())?;
        let start = if rank == 0 {
            0
        } else {
            offset.checked_sub(k - 1).ok_or_else(|| {
                "PanRefV2 backbone coordinate is shorter than its graph overlap".to_string()
            })?
        };
        let end = start
            .checked_add(unitig.sequence.len())
            .ok_or_else(|| "PanRefV2 backbone coordinate overflow".to_string())?;
        coordinates.push(BackboneCoordinate {
            node,
            orientation,
            start,
            end,
        });
        offset = end;
    }
    if offset != backbone.sequence.len() {
        return Err("PanRefV2 backbone coordinates do not match sequence length".into());
    }
    Ok(coordinates)
}

fn write_backbone_coordinates(
    out: &mut BufWriter<File>,
    stable_id: &str,
    graph_name: &str,
    unitigs: &[Unitig],
    k: usize,
    backbone: &ResolvedBackbone,
) -> AppResult<()> {
    for (rank, coordinate) in backbone_coordinates(unitigs, k, backbone)?
        .into_iter()
        .enumerate()
    {
        writeln!(
            out,
            "{stable_id}\t{graph_name}\t{graph_name}_U{}\t{}\t{}\t{}\t{}",
            coordinate.node, coordinate.orientation, rank, coordinate.start, coordinate.end,
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(())
}

fn write_gfa_record(
    out: &mut BufWriter<File>,
    graph_name: &str,
    unitigs: &[Unitig],
    edges: &[(usize, usize)],
    k: usize,
) -> AppResult<()> {
    for (id, unitig) in unitigs.iter().enumerate() {
        writeln!(
            out,
            "S\t{graph_name}_U{id}\t{}\tKC:i:{}",
            String::from_utf8_lossy(&unitig.sequence),
            unitig.kmer_count
        )
        .map_err(|e| e.to_string())?;
    }
    for &(from, to) in edges {
        writeln!(
            out,
            "L\t{graph_name}_U{from}\t+\t{graph_name}_U{to}\t+\t{}M",
            k - 1
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
    let mut ledgers = Vec::with_capacity(samples.len());
    let spool_dir = CandidateSpool::create(&output_dir)?;

    // One raw FASTQ scan per sample. Core reads stay bounded in memory; all
    // one-mate candidates are append-only spool records until the core freezes.
    for (sample_index, sample) in samples.iter().enumerate() {
        let spool_path = spool_dir.path.join(format!("{sample_index:05}.bin"));
        let spool_file = File::create(&spool_path).map_err(|e| e.to_string())?;
        let mut candidate_out = BufWriter::new(spool_file);
        let mut spool_error = None;
        let mut ledger = SampleLedger::new(catalog.loci.len(), spool_path);
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
            true,
            |locus, strong, first, second| {
                if spool_error.is_some() {
                    return;
                }
                let id = locus as usize;
                let pair = LedgerPair::from_records(first, second);
                if strong {
                    let fingerprint = pair.fingerprint();
                    let bucket = ledger.bucket_mut(id);
                    if bucket.core.len() < MAX_CORE_PAIRS_PER_SAMPLE_LOCUS
                        && bucket.core_seen.insert(fingerprint)
                    {
                        bucket.core.push(pair);
                    }
                } else if let Err(error) = write_candidate_record(&mut candidate_out, id, &pair) {
                    spool_error = Some(error);
                }
            },
        )?;
        candidate_out.flush().map_err(|e| e.to_string())?;
        if let Some(error) = spool_error {
            return Err(error);
        }
        summary.candidate_spool_bytes = fs::metadata(&ledger.candidate_spool)
            .map_err(|e| e.to_string())?
            .len();
        summary.strong = stats.strong_pairs;
        summary.candidate_rescues = stats.rescued_pairs;
        summary.ambiguous = stats.ambiguous_pairs;
        summary.loci_with_core = ledger
            .buckets
            .iter()
            .filter(|bucket| bucket.as_ref().is_some_and(|entry| !entry.core.is_empty()))
            .count();
        ledgers.push(ledger);
        summaries.push(summary);
    }

    // Freeze paired core before evaluating rescues, preserving the previous
    // non-self-reinforcing acceptance rule without a second FASTQ decode.
    for ledger in &ledgers {
        for (id, bucket) in ledger.buckets.iter().enumerate() {
            let Some(bucket) = bucket else { continue };
            if !bucket.core.is_empty() {
                accumulators[id].sample_support += 1;
            }
            for pair in &bucket.core {
                accumulators[id].add_pair(&pair.first, pair.second.as_deref());
                accumulators[id].core_pairs += 1;
            }
        }
    }
    let core_gates = accumulators
        .iter()
        .map(|entry| entry.counter.clone())
        .collect::<Vec<_>>();
    for (ledger, summary) in ledgers.iter_mut().zip(summaries.iter_mut()) {
        let file = File::open(&ledger.candidate_spool).map_err(|e| e.to_string())?;
        let mut input = BufReader::new(file);
        while let Some((id, pair)) = read_candidate_record(&mut input, catalog.loci.len())? {
            let touches_core = core_gates[id].has_supported_kmer(&pair.first, MIN_KMER_COUNT)
                || pair
                    .second
                    .as_deref()
                    .is_some_and(|read| core_gates[id].has_supported_kmer(read, MIN_KMER_COUNT));
            if !touches_core {
                continue;
            }
            let fingerprint = pair.fingerprint();
            let bucket = ledger.bucket_mut(id);
            if bucket.accepted_rescues.len() >= MAX_RESCUE_PAIRS_PER_SAMPLE_LOCUS
                || !bucket.accepted_seen.insert(fingerprint)
            {
                continue;
            }
            accumulators[id].add_pair(&pair.first, pair.second.as_deref());
            accumulators[id].rescued_pairs += 1;
            summary.accepted_rescues += 1;
            bucket.accepted_rescues.push(pair);
        }
        drop(input);
        fs::remove_file(&ledger.candidate_spool).map_err(|e| e.to_string())?;
    }
    drop(spool_dir);

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
    let graph_edge_lists = unitigs
        .iter()
        .map(|nodes| unitig_edges(nodes, K))
        .collect::<Vec<_>>();
    let graph_edges = graph_edge_lists
        .iter()
        .map(|edges| edges.iter().copied().collect::<HashSet<_>>())
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
                .map(|_| Vec::<(usize, u32)>::new())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    // Replay only the accepted target ledger. This replaces the old third
    // decompression/recruitment pass and keeps sample colour evidence sparse.
    for (sample_index, ledger) in ledgers.iter().enumerate() {
        let mut local_depth = unitigs
            .iter()
            .map(|nodes| vec![0_u32; nodes.len()])
            .collect::<Vec<_>>();
        for (id, bucket) in ledger.buckets.iter().enumerate() {
            let Some(bucket) = bucket else { continue };
            for pair in bucket.core.iter().chain(&bucket.accepted_rescues) {
                let left_path = read_unitig_path(&pair.first, &unitig_seed_indexes[id]);
                let right_path = pair
                    .second
                    .as_deref()
                    .map(|read| read_unitig_path(read, &unitig_seed_indexes[id]))
                    .unwrap_or_default();
                for &node in left_path.iter().chain(&right_path) {
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
                        for edge in [(left, right), (right, left)] {
                            if graph_edges[id].contains(&edge) {
                                pe_support[id][left] += 1;
                                pe_support[id][right] += 1;
                                *edge_support[id].entry(edge).or_default() += 1;
                            }
                        }
                    }
                }
            }
        }
        for (locus, depths) in local_depth.into_iter().enumerate() {
            for (unitig, depth) in depths.into_iter().enumerate() {
                if depth > 0 {
                    unitig_sample_depth[locus][unitig].push((sample_index, depth));
                }
            }
        }
    }
    let unitig_sample_support = unitig_sample_depth
        .iter()
        .map(|nodes| {
            nodes
                .iter()
                .map(|depths| depths.len() as u32)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let unitig_depth_stability = unitig_sample_depth
        .iter()
        .map(|nodes| {
            nodes
                .iter()
                .map(|depths| {
                    stable_depth_score(&depths.iter().map(|(_, depth)| *depth).collect::<Vec<_>>())
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut colors = BufWriter::new(
        File::create(output_dir.join("unitig_color_evidence.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(colors, "locus\tunitig\tsample\tread_depth").map_err(|e| e.to_string())?;
    for (locus_id, nodes) in unitig_sample_depth.iter().enumerate() {
        for (unitig_id, sample_depths) in nodes.iter().enumerate() {
            for &(sample_id, depth) in sample_depths {
                writeln!(
                    colors,
                    "{}\tL{}_U{}\t{}\t{}",
                    catalog.loci[locus_id].name,
                    locus_id,
                    unitig_id,
                    summaries[sample_id].name,
                    depth
                )
                .map_err(|e| e.to_string())?;
            }
        }
    }
    let mut sample_paths = BufWriter::new(
        File::create(output_dir.join("sample_backbone_paths.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        sample_paths,
        "sample\tgraph\tstatus\tcovered_nodes\trequired_nodes\tsupported_edges\tpath"
    )
    .map_err(|e| e.to_string())?;
    let mut bubble_qc =
        BufWriter::new(File::create(output_dir.join("bubble_qc.tsv")).map_err(|e| e.to_string())?);
    writeln!(bubble_qc, "locus\tgraph\tentry_unitig\texit_unitig\tstatus\talternative_branches\talternative_nodes\talternative_min_samples\talternative_min_edge_support\tcanonical_min_samples\tcanonical_min_edge_support")
        .map_err(|e| e.to_string())?;
    let mut backbone_manifest = BufWriter::new(
        File::create(output_dir.join("backbone_manifest.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(backbone_manifest, "stable_id\tlocus\tfasta_record\tgraph\tgfa_path\torientation\tassembly_k\tsequence_length\tsequence_sha256\tnode_count")
        .map_err(|e| e.to_string())?;
    let mut backbone_coordinates_out = BufWriter::new(
        File::create(output_dir.join("backbone_coordinates.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        backbone_coordinates_out,
        "stable_id\tgraph\tunitig\torientation\tpath_rank\tstart_0based\tend_0based"
    )
    .map_err(|e| e.to_string())?;

    let mut recruitment = BufWriter::new(
        File::create(output_dir.join("recruitment_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(recruitment, "sample\tstrong_pairs\tcandidate_single_mate_pairs\tcandidate_spool_bytes\taccepted_graph_rescues\tambiguous_pairs\tloci_with_core")
        .map_err(|e| e.to_string())?;
    for row in &summaries {
        writeln!(
            recruitment,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.name,
            row.strong,
            row.candidate_rescues,
            row.candidate_spool_bytes,
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
    writeln!(graph, "H\tVN:Z:1.0\tTS:Z:GeneMiner2-PanRefV2.2").map_err(|e| e.to_string())?;
    let mut report = BufWriter::new(
        File::create(output_dir.join("locus_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(report, "locus\tstatus\tassembly_k\tsequence_length\tunitigs\tcore_pairs\taccepted_rescues\tsupporting_samples\tpe_links\tbackbone_min_samples\tbackbone_edge_support\tbackbone_min_depth_stability")
        .map_err(|e| e.to_string())?;
    let mut names = HashSet::new();
    let mut written = 0_usize;
    for (id, locus) in catalog.loci.iter().enumerate() {
        write_gfa_record(
            &mut graph,
            &format!("L{id}"),
            &unitigs[id],
            &graph_edge_lists[id],
            K,
        )?;
        let mut active_unitigs = unitigs[id].len();
        let mut pe_links: u64 = pe_support[id].iter().sum::<u64>() / 2;
        let evidence = PathEvidence {
            sample_support: &unitig_sample_support[id],
            pe_support: &pe_support[id],
            depth_stability: &unitig_depth_stability[id],
            edge_support: &edge_support[id],
            edges: &graph_edge_lists[id],
            k: K,
        };
        let mut backbone = assemble_backbone_from_unitigs(&unitigs[id], &locus.records, evidence);
        let mut assembly_k = K;
        let mut adaptive_evidence = None;
        if should_retry_adaptive(&backbone, accumulators[id].core_pairs) {
            let adaptive_unitigs = rebuild_locus_unitigs(&ledgers, id, ADAPTIVE_K);
            let adaptive_edges = unitig_edges(&adaptive_unitigs, ADAPTIVE_K);
            let projected =
                project_locus_evidence(&ledgers, id, &adaptive_unitigs, &adaptive_edges);
            backbone = assemble_backbone_from_unitigs(
                &adaptive_unitigs,
                &locus.records,
                PathEvidence {
                    sample_support: &projected.sample_support,
                    pe_support: &projected.pe_support,
                    depth_stability: &projected.depth_stability,
                    edge_support: &projected.edge_support,
                    edges: &adaptive_edges,
                    k: ADAPTIVE_K,
                },
            );
            if let Some(adaptive_backbone) = &backbone {
                assembly_k = ADAPTIVE_K;
                active_unitigs = adaptive_unitigs.len();
                pe_links = projected.pe_support.iter().sum::<u64>() / 2;
                let graph_name = format!("L{id}_K{ADAPTIVE_K}");
                write_gfa_record(
                    &mut graph,
                    &graph_name,
                    &adaptive_unitigs,
                    &adaptive_edges,
                    ADAPTIVE_K,
                )?;
                for (unitig_id, depths) in projected.sample_depth.iter().enumerate() {
                    for &(sample_id, depth) in depths {
                        writeln!(
                            colors,
                            "{}\t{}_U{}\t{}\t{}",
                            locus.name, graph_name, unitig_id, summaries[sample_id].name, depth
                        )
                        .map_err(|e| e.to_string())?;
                    }
                }
                write_sample_backbone_paths(
                    &mut sample_paths,
                    &ledgers,
                    id,
                    &adaptive_unitigs,
                    adaptive_backbone,
                    &graph_name,
                    &summaries,
                )?;
                write_bubble_qc(
                    &mut bubble_qc,
                    BubbleQcView {
                        locus: &locus.name,
                        graph_name: &graph_name,
                        unitigs: &adaptive_unitigs,
                        edges: &adaptive_edges,
                        backbone: adaptive_backbone,
                        sample_support: &projected.sample_support,
                        edge_support: &projected.edge_support,
                    },
                )?;
                let stable_id = stable_backbone_id(&locus.name, &adaptive_backbone.sequence);
                write_backbone_coordinates(
                    &mut backbone_coordinates_out,
                    &stable_id,
                    &graph_name,
                    &adaptive_unitigs,
                    ADAPTIVE_K,
                    adaptive_backbone,
                )?;
                adaptive_evidence = Some(projected);
            }
        }
        let graph_name = if assembly_k == K {
            format!("L{id}")
        } else {
            format!("L{id}_K{ADAPTIVE_K}")
        };
        if assembly_k == K {
            if let Some(path) = &backbone {
                write_sample_backbone_paths(
                    &mut sample_paths,
                    &ledgers,
                    id,
                    &unitigs[id],
                    path,
                    &graph_name,
                    &summaries,
                )?;
                write_bubble_qc(
                    &mut bubble_qc,
                    BubbleQcView {
                        locus: &locus.name,
                        graph_name: &graph_name,
                        unitigs: &unitigs[id],
                        edges: &graph_edge_lists[id],
                        backbone: path,
                        sample_support: &unitig_sample_support[id],
                        edge_support: &edge_support[id],
                    },
                )?;
                let stable_id = stable_backbone_id(&locus.name, &path.sequence);
                write_backbone_coordinates(
                    &mut backbone_coordinates_out,
                    &stable_id,
                    &graph_name,
                    &unitigs[id],
                    K,
                    path,
                )?;
            }
        }
        if let Some(path) = &backbone {
            let segments = gfa_path_segments(&graph_name, path);
            writeln!(graph, "P\t{graph_name}_backbone\t{segments}\t*")
                .map_err(|e| e.to_string())?;
        }
        let (active_sample_support, active_depth_stability, active_edge_support): ActiveEvidence<
            '_,
        > = if assembly_k == K {
            (
                &unitig_sample_support[id],
                &unitig_depth_stability[id],
                &edge_support[id],
            )
        } else {
            let projected = adaptive_evidence
                .as_ref()
                .expect("adaptive evidence for adaptive backbone");
            (
                &projected.sample_support,
                &projected.depth_stability,
                &projected.edge_support,
            )
        };
        let backbone_min_samples = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .iter()
                .map(|&node| active_sample_support[node])
                .min()
                .unwrap_or_default()
        });
        let backbone_edge_support = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .windows(2)
                .map(|edge| {
                    active_edge_support
                        .get(&(edge[0], edge[1]))
                        .copied()
                        .unwrap_or_default()
                })
                .sum()
        });
        let backbone_min_depth_stability = backbone.as_ref().map_or(0, |path| {
            path.nodes
                .iter()
                .map(|&node| active_depth_stability[node])
                .min()
                .unwrap_or_default()
        });
        let sequence = backbone.as_ref().map(|path| path.sequence.as_slice());
        let (status, length) = match sequence {
            None if accumulators[id].core_pairs == 0 => ("no_core", 0),
            None if unitigs[id].is_empty() => ("low_coverage", 0),
            None => ("complex", 0),
            Some(seq) if seq.len() < 100 => ("short", seq.len()),
            Some(seq) if accumulators[id].sample_support < 2 => ("low_sample_support", seq.len()),
            Some(seq) => ("pass", seq.len()),
        };
        let mut fasta_record = String::new();
        if let Some(sequence) =
            sequence.filter(|_| status == "pass" || args.panrefv2_include_low_confidence)
        {
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
            fasta_record = name;
            written += 1;
        }
        if let Some(path) = &backbone {
            let stable_id = stable_backbone_id(&locus.name, &path.sequence);
            let orientation = if path.reversed {
                char::from(45_u8)
            } else {
                char::from(43_u8)
            };
            writeln!(
                backbone_manifest,
                "{stable_id}\t{}\t{fasta_record}\t{graph_name}\t{graph_name}_backbone\t{orientation}\t{assembly_k}\t{}\t{}\t{}",
                locus.name,
                path.sequence.len(),
                sha256_hex(&path.sequence),
                path.nodes.len(),
            )
            .map_err(|error| error.to_string())?;
        }
        writeln!(
            report,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            locus.name,
            status,
            assembly_k,
            length,
            active_unitigs,
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
    colors.flush().map_err(|e| e.to_string())?;
    sample_paths.flush().map_err(|e| e.to_string())?;
    bubble_qc.flush().map_err(|e| e.to_string())?;
    backbone_manifest.flush().map_err(|e| e.to_string())?;
    backbone_coordinates_out
        .flush()
        .map_err(|e| e.to_string())?;
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
    use super::{
        backbone_coordinates, build_unitig_seed_index, gfa_path_segments, read_candidate_record,
        read_unitig_path, rebuild_locus_unitigs, sha256_hex, should_retry_adaptive,
        stable_backbone_id, unique_unitig, write_candidate_record, write_sample_backbone_paths,
        BackboneCoordinate, CandidateSpool, LedgerBucket, LedgerPair, ResolvedBackbone,
        SampleLedger, SampleSummary,
    };
    use crate::panref::dbg::Unitig;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write};

    #[test]
    fn stable_backbone_ids_are_content_and_locus_specific() {
        let first = stable_backbone_id("uce-1", b"ACGT");
        assert_eq!(first, stable_backbone_id("uce-1", b"ACGT"));
        assert_ne!(first, stable_backbone_id("uce-2", b"ACGT"));
        assert_ne!(first, stable_backbone_id("uce-1", b"ACGA"));
        assert!(first.starts_with("panrefv2-sha256-"));
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(sha256_hex(b"ACGT").len(), 64);
    }

    #[test]
    fn reverse_backbone_coordinates_follow_emitted_gfa_orientation() {
        let unitigs = vec![
            Unitig {
                sequence: b"ACGT".to_vec(),
                kmer_count: 1,
            },
            Unitig {
                sequence: b"GTAA".to_vec(),
                kmer_count: 1,
            },
        ];
        let coordinates = backbone_coordinates(
            &unitigs,
            3,
            &ResolvedBackbone {
                sequence: b"TTACGT".to_vec(),
                nodes: vec![0, 1],
                reversed: true,
            },
        )
        .unwrap();
        assert_eq!(
            coordinates,
            vec![
                BackboneCoordinate {
                    node: 1,
                    orientation: char::from(45_u8),
                    start: 0,
                    end: 4,
                },
                BackboneCoordinate {
                    node: 0,
                    orientation: char::from(45_u8),
                    start: 2,
                    end: 6,
                },
            ]
        );
    }

    #[test]
    fn gfa_path_uses_reverse_node_order_and_minus_orientation() {
        let path = ResolvedBackbone {
            sequence: b"ACGT".to_vec(),
            nodes: vec![2, 7],
            reversed: true,
        };
        assert_eq!(gfa_path_segments("L4", &path), "L4_U7-,L4_U2-");
    }

    #[test]
    fn candidate_spool_is_removed_on_drop() {
        let root =
            std::env::temp_dir().join(format!("gm2-panrefv21-spool-drop-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let spool = CandidateSpool::create(&root).unwrap();
        let path = spool.path.clone();
        assert!(path.is_dir());
        drop(spool);
        assert!(!path.exists());
        std::fs::remove_dir(root).unwrap();
    }

    #[test]
    fn candidate_spool_round_trips_locus_and_mates() {
        let path =
            std::env::temp_dir().join(format!("gm2-panrefv21-spool-{}.bin", std::process::id()));
        let pair = LedgerPair {
            first: b"ACGTN".to_vec(),
            second: Some(b"TTGCA".to_vec()),
        };
        let mut out = BufWriter::new(File::create(&path).unwrap());
        write_candidate_record(&mut out, 7, &pair).unwrap();
        out.flush().unwrap();
        drop(out);
        let mut input = BufReader::new(File::open(&path).unwrap());
        let decoded = read_candidate_record(&mut input, 8).unwrap().unwrap();
        assert_eq!(decoded.0, 7);
        assert_eq!(decoded.1.first, pair.first);
        assert_eq!(decoded.1.second, pair.second);
        assert!(read_candidate_record(&mut input, 8).unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn adaptive_k_can_recover_when_the_standard_graph_is_empty() {
        let sequence = b"ACGTTGCAAGTCCTAGGATTCGACCTGTAA".to_vec();
        let pair = LedgerPair {
            first: sequence.clone(),
            second: Some(sequence),
        };
        let ledgers = vec![SampleLedger {
            buckets: vec![Some(LedgerBucket {
                core: vec![pair],
                ..Default::default()
            })],
            candidate_spool: std::path::PathBuf::new(),
        }];
        assert!(rebuild_locus_unitigs(&ledgers, 0, 31).is_empty());
        assert!(!rebuild_locus_unitigs(&ledgers, 0, 25).is_empty());
        assert!(should_retry_adaptive(&None, 1));
    }

    #[test]
    fn sample_paths_require_complete_node_and_edge_evidence() {
        let unitigs = vec![
            Unitig {
                sequence: b"AAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
                kmer_count: 1,
            },
            Unitig {
                sequence: b"CCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),
                kmer_count: 1,
            },
        ];
        let ledgers = vec![
            SampleLedger {
                buckets: vec![Some(LedgerBucket {
                    core: vec![LedgerPair {
                        first: unitigs[0].sequence.clone(),
                        second: None,
                    }],
                    ..Default::default()
                })],
                candidate_spool: std::path::PathBuf::new(),
            },
            SampleLedger {
                buckets: vec![Some(LedgerBucket {
                    core: vec![LedgerPair {
                        first: unitigs[1].sequence.clone(),
                        second: Some(unitigs[0].sequence.clone()),
                    }],
                    ..Default::default()
                })],
                candidate_spool: std::path::PathBuf::new(),
            },
        ];
        let summaries = vec![
            SampleSummary {
                name: "partial".into(),
                ..Default::default()
            },
            SampleSummary {
                name: "complete".into(),
                ..Default::default()
            },
        ];
        let backbone = ResolvedBackbone {
            sequence: [unitigs[0].sequence.clone(), unitigs[1].sequence.clone()].concat(),
            nodes: vec![0, 1],
            reversed: false,
        };
        let path = std::env::temp_dir().join(format!(
            "gm2-panrefv23-sample-paths-{}.tsv",
            std::process::id()
        ));
        let mut out = BufWriter::new(File::create(&path).unwrap());
        write_sample_backbone_paths(&mut out, &ledgers, 0, &unitigs, &backbone, "L0", &summaries)
            .unwrap();
        out.flush().unwrap();
        drop(out);
        let text = std::fs::read_to_string(&path).unwrap();
        let rows = text.lines().collect::<Vec<_>>();
        assert_eq!(rows[0], "partial\tL0\tpartial\t1\t2\t0\t");
        assert_eq!(
            rows[1],
            "complete\tL0\tbackbone_path\t2\t2\t1\tL0_U0+,L0_U1+"
        );
        std::fs::remove_file(path).unwrap();
    }

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
