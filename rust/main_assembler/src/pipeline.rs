use crate::assembly::{
    add_pe_branch_support, add_read_slices, assemble_seed, branch_edges, build_graph_from_counts,
    calculate_read_support, compare_contigs, count_assemble_chunk_parallel,
    filter_and_weight_graph, SortedKmerCounts,
};
use crate::hash::{HashMap, HashSet};
use crate::io_utils::{
    find_filtered, for_each_sequence_chunk, has_interleaved_pairs, load_or_build_reference_kmers,
    minimum_sequence_length, read_fasta, read_linked_fragments,
};
use crate::model::{Args, AssemblyMode, ContigRecord, GraphFormat, LocusResult, LocusTask};
use crate::seq::{calculate_auto_k, reverse_complement, reverse_complement_kmer};
use crate::unitig::write_graphs;
use std::fs::{self, File, OpenOptions};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

pub const SUMMARY_HEADER: &str = "locus,status,accepted,rejection_reason,selected_contig_length,read_supported_span,slice_supported_bases,slice_support_breadth,max_slice_support_gap,read_count,unique_read_count,multi_mapping_read_count,read_density,unique_read_density,support_fraction,flank_balance,kmer_median_depth,kmer_depth_cv,kmer_max_depth_ratio,candidate_count,low_quality";

#[derive(Default)]
pub struct ProfileStats {
    processed_loci: AtomicU64,
    completed_loci: AtomicU64,
    failed_loci: AtomicU64,
    skipped_loci: AtomicU64,
    reads_decode_ns: AtomicU64,
    kmer_count_ns: AtomicU64,
    graph_filter_ns: AtomicU64,
    seed_path_ns: AtomicU64,
    qc_write_ns: AtomicU64,
}

impl ProfileStats {
    fn add(counter: &AtomicU64, started: Instant) {
        counter.fetch_add(
            started.elapsed().as_nanos().min(u64::MAX as u128) as u64,
            Ordering::Relaxed,
        );
    }
    fn add_duration(counter: &AtomicU64, elapsed: Duration) {
        counter.fetch_add(
            elapsed.as_nanos().min(u64::MAX as u128) as u64,
            Ordering::Relaxed,
        );
    }
    fn record_reads(&self, started: Instant) {
        Self::add(&self.reads_decode_ns, started);
    }
    fn record_reads_duration(&self, elapsed: Duration) {
        Self::add_duration(&self.reads_decode_ns, elapsed);
    }
    fn record_kmers(&self, started: Instant) {
        Self::add(&self.kmer_count_ns, started);
    }
    fn record_kmers_duration(&self, elapsed: Duration) {
        Self::add_duration(&self.kmer_count_ns, elapsed);
    }
    fn record_graph(&self, started: Instant) {
        Self::add(&self.graph_filter_ns, started);
    }
    fn record_path(&self, started: Instant) {
        Self::add(&self.seed_path_ns, started);
    }
    fn record_qc(&self, started: Instant) {
        Self::add(&self.qc_write_ns, started);
    }
    fn record_locus(&self, result: &LocusResult) {
        self.processed_loci.fetch_add(1, Ordering::Relaxed);
        if result.skipped {
            self.skipped_loci.fetch_add(1, Ordering::Relaxed);
        } else if matches!(result.status.as_str(), "success" | "low quality") {
            self.completed_loci.fetch_add(1, Ordering::Relaxed);
        } else {
            self.failed_loci.fetch_add(1, Ordering::Relaxed);
        }
    }
}

pub fn write_assembly_profile(output: &Path, profile: &ProfileStats) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(output.join("assembly_profile.tsv"))?);
    writeln!(
        out,
        "stage	seconds	processed_loci	completed_loci	failed_loci	skipped_loci"
    )?;
    let processed = profile.processed_loci.load(Ordering::Relaxed);
    let completed = profile.completed_loci.load(Ordering::Relaxed);
    let failed = profile.failed_loci.load(Ordering::Relaxed);
    let skipped = profile.skipped_loci.load(Ordering::Relaxed);
    for (name, value) in [
        ("reads_decode", &profile.reads_decode_ns),
        ("kmer_count", &profile.kmer_count_ns),
        ("graph_filter_weight", &profile.graph_filter_ns),
        ("seed_path_search", &profile.seed_path_ns),
        ("contig_qc_write", &profile.qc_write_ns),
    ] {
        writeln!(
            out,
            "{name}	{:.6}	{processed}	{completed}	{failed}	{skipped}",
            value.load(Ordering::Relaxed) as f64 / 1e9
        )?;
    }
    out.flush()
}

// 多线程写同一个日志得排队，不然字儿搅一块就没法看了。
pub fn log_line(output: &Path, lock: &Mutex<()>, message: &str) {
    let _guard = lock.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Ok(mut log) = OpenOptions::new()
        .create(true)
        .append(true)
        .open(output.join("log.txt"))
    {
        let _ = writeln!(log, "{message}");
    }
    println!("{message}");
}

fn remove_if_exists(path: &Path) {
    if path.is_file() {
        let _ = fs::remove_file(path);
    }
}

// 某个位点重跑或失败时，只扫它自己的旧产物，别误伤邻居。
fn clean_locus_outputs(args: &Args, key: &str) {
    remove_if_exists(&args.output.join("results").join(format!("{key}.fasta")));
    remove_if_exists(
        &args
            .output
            .join("results")
            .join(format!("{key}.its2_support.tsv")),
    );
    remove_if_exists(&args.output.join("contigs_all").join(format!("{key}.fasta")));
    remove_if_exists(
        &args
            .output
            .join("contigs_all_low")
            .join(format!("{key}.fasta")),
    );
}

fn format_header(contig: &ContigRecord, mode: AssemblyMode, prefix: &str) -> String {
    let mut header = format!(
        ">{prefix}_{}_{}_{}_{}_{}_span_{}",
        contig.sequence.len(),
        contig.seed_count,
        contig.position,
        contig.weight,
        contig.read_count,
        contig.supported_span
    );
    if mode == AssemblyMode::Uce {
        header.push_str(&format!("_supported_{}", contig.supported_bases));
    }
    if mode == AssemblyMode::Its2 {
        header.push_str(&format!(
            "_fragments_{}_paired_{}_diagnostic_{}_em_{:.3}_abundance_{:.6}",
            contig.fragment_support,
            contig.paired_fragment_support,
            contig.diagnostic_fragment_support,
            contig.em_fragment_support,
            contig.em_abundance
        ));
    }
    if mode == AssemblyMode::Its2 && !contig.label.is_empty() {
        header = format!(">{}_{}", contig.label, &header[1..]);
    }
    header.push_str(&format!("_balance_{:.3}", contig.flank_balance));
    header
}

// 把候选按等级落盘，标题里带够指标，后续人工复核不抓瞎。
fn write_contigs(
    path: &Path,
    contigs: &[ContigRecord],
    mode: AssemblyMode,
    prefix: &str,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for contig in contigs {
        writeln!(writer, "{}", format_header(contig, mode, prefix))?;
        writeln!(writer, "{}", String::from_utf8_lossy(&contig.sequence))?;
    }
    writer.flush()
}

// ITS2 的 mate 能不能落在同一候选上，是区分真支持和瞎碰的重要证据。
fn mate_matches_candidate(mate: &[u8], candidate: &[u8], slice_len: usize) -> bool {
    if mate.is_empty() || candidate.is_empty() {
        return false;
    }
    let width = slice_len.min(mate.len());
    let start = (mate.len() - width) / 2;
    let slice = &mate[start..start + width];
    let reverse = reverse_complement(slice);
    candidate.windows(width).any(|window| window == slice)
        || candidate
            .windows(width)
            .any(|window| window == reverse.as_slice())
}

fn append_equivalence_member(members: &mut String, label: &str) {
    if label.is_empty() || members.split(";").any(|member| member == label) {
        return;
    }
    if !members.is_empty() {
        members.push(';');
    }
    members.push_str(label);
}

fn prefer_its2_candidate(
    existing: &ContigRecord,
    candidate: &ContigRecord,
    mode: AssemblyMode,
) -> bool {
    if existing.label.is_empty() != candidate.label.is_empty() {
        return !candidate.label.is_empty();
    }
    compare_contigs(existing, candidate, mode).is_lt()
}

// ITS2 候选分组、配对支持和 EM 丰度都在这里补齐。
fn annotate_its2_candidates(
    candidates: &mut [ContigRecord],
    fragments: &[Vec<Vec<u8>>],
    slice_len: usize,
) {
    let mut compatibilities = Vec::new();
    for fragment in fragments {
        let compatible: Vec<(usize, usize)> = candidates
            .iter()
            .enumerate()
            .filter_map(|(index, candidate)| {
                let matching_mates = fragment
                    .iter()
                    .filter(|mate| mate_matches_candidate(mate, &candidate.sequence, slice_len))
                    .count();
                (matching_mates > 0).then_some((index, matching_mates))
            })
            .collect();
        for (index, matching_mates) in &compatible {
            candidates[*index].fragment_support += 1;
            if *matching_mates >= 2 {
                candidates[*index].paired_fragment_support += 1;
            }
        }
        let indices: Vec<usize> = compatible.iter().map(|(index, _)| *index).collect();
        if indices.len() == 1 {
            candidates[indices[0]].diagnostic_fragment_support += 1;
        }
        if !indices.is_empty() {
            compatibilities.push(indices);
        }
    }
    if compatibilities.is_empty() {
        return;
    }
    let mut abundance: Vec<f64> = candidates
        .iter()
        .map(|candidate| candidate.fragment_support.max(1) as f64)
        .collect();
    let normalizer: f64 = abundance.iter().sum();
    for value in &mut abundance {
        *value /= normalizer;
    }
    let mut expected = vec![0.0; candidates.len()];
    for _ in 0..64 {
        expected.fill(0.0);
        for compatible in &compatibilities {
            let denominator: f64 = compatible.iter().map(|index| abundance[*index]).sum();
            if denominator > 0.0 {
                for index in compatible {
                    expected[*index] += abundance[*index] / denominator;
                }
            }
        }
        let total: f64 = expected.iter().sum();
        if total > 0.0 {
            for (value, count) in abundance.iter_mut().zip(expected.iter()) {
                *value = *count / total;
            }
        }
    }
    for (candidate, (count, value)) in candidates
        .iter_mut()
        .zip(expected.into_iter().zip(abundance))
    {
        candidate.em_fragment_support = count;
        candidate.em_abundance = value;
        candidate.accepted = candidate.fragment_support >= 2;
        candidate.rejection_reason = if candidate.accepted {
            String::new()
        } else {
            "low_fragment_or_diagnostic_support".to_string()
        };
    }
}

fn its2_status(contig: &ContigRecord) -> &'static str {
    if !contig.accepted {
        "LOW_SUPPORT"
    } else if contig.diagnostic_fragment_support > 0 && contig.paired_fragment_support > 0 {
        "PASS"
    } else {
        "EQUIVALENCE_GROUP"
    }
}

// 把 ITS2 的支持证据另存表，主 FASTA 不塞得乱七八糟。
fn write_its2_support(path: &Path, contigs: &[ContigRecord]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "variant,equivalence_members,length,fragment_support,paired_fragment_support,diagnostic_fragment_support,em_fragment_support,em_abundance,status")?;
    for (index, contig) in contigs.iter().enumerate() {
        let label = if contig.label.is_empty() {
            format!("assembly_{}", index + 1)
        } else {
            contig.label.clone()
        };
        writeln!(
            writer,
            "{},{},{},{},{},{},{:.6},{:.8},{}",
            label,
            if contig.equivalence_members.is_empty() {
                label.clone()
            } else {
                contig.equivalence_members.clone()
            },
            contig.sequence.len(),
            contig.fragment_support,
            contig.paired_fragment_support,
            contig.diagnostic_fragment_support,
            contig.em_fragment_support,
            contig.em_abundance,
            its2_status(contig)
        )?;
    }
    writer.flush()
}

pub fn run_manifest(args: &Args, tasks: &[LocusTask]) -> io::Result<String> {
    let mut hasher = DefaultHasher::new();
    args.kmer_size.hash(&mut hasher);
    args.kmer_min.hash(&mut hasher);
    args.kmer_max.hash(&mut hasher);
    args.error_limit.hash(&mut hasher);
    args.iteration.hash(&mut hasher);
    args.min_coverage.to_bits().hash(&mut hasher);
    args.soft_boundary.hash(&mut hasher);
    args.assembly_mode.hash(&mut hasher);
    args.side_candidates.hash(&mut hasher);
    args.path_strategy.hash(&mut hasher);
    args.backbone_lookahead.hash(&mut hasher);
    args.max_contig_length.hash(&mut hasher);
    args.min_read_density.to_bits().hash(&mut hasher);
    args.density_check_min_length.hash(&mut hasher);
    args.max_depth_cv.to_bits().hash(&mut hasher);
    args.max_depth_ratio.to_bits().hash(&mut hasher);
    for task in tasks {
        task.key.hash(&mut hasher);
        let reference = fs::metadata(&task.reference_path)?;
        reference.len().hash(&mut hasher);
        reference
            .modified()
            .ok()
            .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|time| time.as_nanos())
            .unwrap_or(0)
            .hash(&mut hasher);
        match find_filtered(&args.output, &task.key) {
            Some((path, fasta)) => {
                path.hash(&mut hasher);
                fasta.hash(&mut hasher);
                let filtered = fs::metadata(path)?;
                filtered.len().hash(&mut hasher);
                filtered
                    .modified()
                    .ok()
                    .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
                    .map(|time| time.as_nanos())
                    .unwrap_or(0)
                    .hash(&mut hasher);
            }
            None => false.hash(&mut hasher),
        }
    }
    Ok(format!(
        "GM2_ASSEMBLER_RUN_MANIFEST_V1\n{:016x}\n",
        hasher.finish()
    ))
}

pub fn process_locus(
    args: &Args,
    task: &LocusTask,
    completed: &HashSet<String>,
    log_lock: &Arc<Mutex<()>>,
    profile: Option<&ProfileStats>,
) -> LocusResult {
    // 单个位点出岔子就收拾它自己的输出，别把整批组装拖下水。
    let result = match process_locus_inner(args, task, completed, log_lock, profile) {
        Ok(result) => result,
        Err(error) => {
            clean_locus_outputs(args, &task.key);
            log_line(
                &args.output,
                log_lock,
                &format!("error assembling {}: {error}", task.key),
            );
            LocusResult::failure(&task.key, "error")
        }
    };
    if let Some(profile) = profile {
        profile.record_locus(&result);
    }
    result
}

fn process_locus_inner(
    args: &Args,
    task: &LocusTask,
    completed: &HashSet<String>,
    log_lock: &Arc<Mutex<()>>,
    profile: Option<&ProfileStats>,
) -> io::Result<LocusResult> {
    // 一个 locus 从过滤 reads 到结果文件，全流程在这旮沓串起来。
    let best_path = args
        .output
        .join("results")
        .join(format!("{}.fasta", task.key));
    let all_path = args
        .output
        .join("contigs_all")
        .join(format!("{}.fasta", task.key));
    let low_path = args
        .output
        .join("contigs_all_low")
        .join(format!("{}.fasta", task.key));

    if completed.contains(&task.key)
        && best_path.is_file()
        && fs::metadata(&best_path).is_ok_and(|metadata| metadata.len() > 0)
    {
        return Ok(LocusResult {
            key: task.key.clone(),
            status: "skipped".to_string(),
            skipped: true,
            ..LocusResult::default()
        });
    }

    let Some((filtered_path, fasta)) = find_filtered(&args.output, &task.key) else {
        clean_locus_outputs(args, &task.key);
        return Ok(LocusResult::failure(&task.key, "no filtered file"));
    };
    let paired_fragments =
        args.assembly_mode == AssemblyMode::Uce && has_interleaved_pairs(&filtered_path, fasta)?;
    let count_chunk_size = if paired_fragments {
        (args.read_chunk_size.max(2) / 2) * 2
    } else {
        args.read_chunk_size
    };
    let reads_started = Instant::now();
    let linked_fragments = if args.assembly_mode == AssemblyMode::Its2 {
        read_linked_fragments(&filtered_path, fasta)?
    } else {
        Vec::new()
    };
    let minimum = minimum_sequence_length(&filtered_path, fasta, args.read_chunk_size)?;
    let slice_len = minimum.map_or(0, |length| ((length as f64 * 0.9) as usize).max(1));
    if slice_len == 0 {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            &format!("No reads were obtained for gene {}", task.key),
        );
        return Ok(LocusResult::failure(&task.key, "no reads"));
    }

    let mut reads = HashMap::new();
    let fixed_k = args.kmer_size != 0;
    if !fixed_k {
        for_each_sequence_chunk(&filtered_path, fasta, args.read_chunk_size, |chunk| {
            add_read_slices(&mut reads, chunk, slice_len);
            Ok(())
        })?;
    }
    if let Some(profile) = profile {
        profile.record_reads(reads_started);
    }
    if !fixed_k && reads.is_empty() {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            &format!("No reads were obtained for gene {}", task.key),
        );
        return Ok(LocusResult::failure(&task.key, "no reads"));
    }

    let kmer_started = Instant::now();
    let reference_records = read_fasta(&task.reference_path)?;
    let reference_sequences: Vec<Vec<u8>> = reference_records
        .iter()
        .map(|(_, sequence)| sequence.clone())
        .collect();
    let current_k = if fixed_k {
        args.kmer_size
    } else {
        // 用户没指定 k 就按 reads 和参考自动挑，别死抱默认值。
        calculate_auto_k(
            &reference_sequences,
            &reads,
            slice_len,
            args.kmer_min,
            args.kmer_max,
            args.error_limit,
        )
    };
    if current_k == 0 || current_k > 63 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("assembly k-mer must be between 1 and 63, got {current_k}"),
        ));
    }

    let soft_boundary = if args.soft_boundary < 0 {
        slice_len / 2
    } else {
        args.soft_boundary as usize
    };
    log_line(
        &args.output,
        log_lock,
        &format!("Use k={current_k} for assembling gene {}.", task.key),
    );
    log_line(
        &args.output,
        log_lock,
        &format!("Assembling {} {} / {}", task.key, task.ordinal, task.total),
    );

    let reference = load_or_build_reference_kmers(
        &task.reference_path,
        &reference_records,
        current_k,
        args.reference_cache_dir.as_deref(),
    )?;
    let mut count_runs = SortedKmerCounts::default();
    if fixed_k {
        if let Some(profile) = profile {
            profile.record_kmers(kmer_started);
        }
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            let reads_chunk_started = Instant::now();
            add_read_slices(&mut reads, chunk, slice_len);
            if let Some(profile) = profile {
                profile.record_reads_duration(reads_chunk_started.elapsed());
            }
            let kmers_chunk_started = Instant::now();
            count_runs.push(count_assemble_chunk_parallel(
                chunk,
                current_k,
                args.kmer_count_threads,
                paired_fragments,
                paired_fragments.then_some(&reference),
            ));
            if let Some(profile) = profile {
                profile.record_kmers_duration(kmers_chunk_started.elapsed());
            }
            Ok(())
        })?;
        if reads.is_empty() {
            clean_locus_outputs(args, &task.key);
            log_line(
                &args.output,
                log_lock,
                &format!("No reads were obtained for gene {}", task.key),
            );
            return Ok(LocusResult::failure(&task.key, "no reads"));
        }
    } else {
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            count_runs.push(count_assemble_chunk_parallel(
                chunk,
                current_k,
                args.kmer_count_threads,
                paired_fragments,
                paired_fragments.then_some(&reference),
            ));
            Ok(())
        })?;
        if let Some(profile) = profile {
            profile.record_kmers(kmer_started);
        }
    }
    let graph_started = Instant::now();
    let min_fragment_support = u32::from(paired_fragments) * 2;
    let counts = count_runs.into_counts();
    // Both candidate graphs share one k-mer count. The core graph supplies a
    // safe no-PE path; the assisted graph only needs to beat it on the normal
    // UCE evidence ranking to become the reported contig.
    let mut core_graph = build_graph_from_counts(counts.clone(), &reference, args.error_limit, 0);
    filter_and_weight_graph(&mut core_graph, args.error_limit, task.reference_count, 0);
    let mut graph =
        build_graph_from_counts(counts, &reference, args.error_limit, min_fragment_support);
    filter_and_weight_graph(
        &mut graph,
        args.error_limit,
        task.reference_count,
        min_fragment_support,
    );
    let mut branch_support = branch_edges(&graph, current_k);
    if paired_fragments && !branch_support.is_empty() {
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            add_pe_branch_support(chunk, current_k, &reference, &mut branch_support);
            Ok(())
        })?;
    }
    let (gfa, dot) = match args.graph_format {
        GraphFormat::None => (false, false),
        GraphFormat::Gfa => (true, false),
        GraphFormat::Dot => (false, true),
        GraphFormat::Both => (true, true),
    };
    if gfa || dot {
        write_graphs(
            &args.output.join("assembly_graphs"),
            &task.key,
            &graph,
            current_k,
            gfa,
            dot,
        )?;
    }
    if let Some(profile) = profile {
        profile.record_graph(graph_started);
    }
    if graph.len() < 3 {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            "Could not get enough reads from filter.",
        );
        return Ok(LocusResult::failure(
            &task.key,
            "insufficient genomic kmers",
        ));
    }

    let path_started = Instant::now();
    let mut seeds: Vec<(u128, i64, i32, i64)> = graph
        .iter()
        .filter(|(_, value)| value.position > 1 && value.position < 1000 && !value.is_reverse)
        .map(|(kmer, value)| (*kmer, value.depth, value.position, value.reference_weight))
        .collect();
    seeds.sort_by(|left, right| {
        right
            .3
            .cmp(&left.3)
            .then_with(|| right.1.cmp(&left.1))
            .then_with(|| left.0.cmp(&right.0))
    });
    if seeds.is_empty() {
        clean_locus_outputs(args, &task.key);
        log_line(&args.output, log_lock, "Could not get enough seeds.");
        return Ok(LocusResult::failure(&task.key, "no seed"));
    }

    let original_seed_count = seeds.len();
    let seed_set: HashSet<u128> = seeds.iter().map(|seed| seed.0).collect();
    let mut accepted = Vec::new();
    let mut rejected = Vec::new();
    let mut core_accepted = Vec::new();

    while seeds.len() > original_seed_count / 2 {
        let seed = seeds[0].0;
        let (mut candidates, used_kmers, _) = assemble_seed(
            args,
            &reads,
            slice_len,
            &graph,
            seed,
            current_k,
            soft_boundary,
            paired_fragments.then_some(&branch_support),
        );
        let mut core_start = None;
        let core_used_kmers = if paired_fragments {
            let (mut core_candidates, core_used_kmers, _) = assemble_seed(
                args,
                &reads,
                slice_len,
                &core_graph,
                seed,
                current_k,
                soft_boundary,
                None,
            );
            core_start = Some(candidates.len());
            candidates.append(&mut core_candidates);
            Some(core_used_kmers)
        } else {
            None
        };
        let seed_coverage = core_used_kmers.as_ref().unwrap_or(&used_kmers);
        let used_seed_count = seed_set.intersection(seed_coverage).count();
        for candidate in &mut candidates {
            candidate.seed_count = used_seed_count;
        }

        seeds.retain(|candidate| {
            !seed_coverage.contains(&candidate.0)
                && !seed_coverage.contains(&reverse_complement_kmer(candidate.0, current_k))
        });
        for (candidate_index, candidate) in candidates.into_iter().enumerate() {
            let is_core_candidate = core_start.is_some_and(|start| candidate_index >= start);
            if args.assembly_mode == AssemblyMode::Uce {
                if candidate.accepted {
                    if is_core_candidate {
                        core_accepted.push(candidate.clone());
                    }
                    accepted.push(candidate);
                } else {
                    rejected.push(candidate);
                }
            } else if args.assembly_mode == AssemblyMode::Its2
                || candidate.read_count as usize * slice_len > candidate.sequence.len()
            {
                accepted.push(candidate);
            } else {
                rejected.push(candidate);
            }
        }
    }

    if let Some(profile) = profile {
        profile.record_path(path_started);
    }
    let qc_started = Instant::now();
    if args.assembly_mode == AssemblyMode::Its2 {
        accepted.retain(|candidate| candidate.sequence.len() <= 400);
        for (label, sequence) in &reference_records {
            let support = calculate_read_support(sequence, slice_len, &reads);
            if support.total_read_count > 0 {
                accepted.push(ContigRecord {
                    label: label.clone(),
                    equivalence_members: label.clone(),
                    sequence: sequence.clone(),
                    read_count: support.total_read_count,
                    supported_span: support.supported_extent,
                    supported_bases: support.supported_bases,
                    support_breadth: support.breadth,
                    max_support_gap: support.max_gap,
                    flank_balance: support.flank_balance,
                    unique_read_count: support.unique_read_count,
                    multi_mapping_read_count: support.multi_mapping_read_count,
                    read_density: if sequence.is_empty() {
                        0.0
                    } else {
                        support.total_read_count as f64 / sequence.len() as f64
                    },
                    support_fraction: if sequence.is_empty() {
                        0.0
                    } else {
                        support.supported_extent as f64 / sequence.len() as f64
                    },
                    ..ContigRecord::default()
                });
            }
        }
        let mut unique: HashMap<Vec<u8>, ContigRecord> = HashMap::new();
        for mut candidate in accepted.drain(..).chain(rejected.drain(..)) {
            if candidate.equivalence_members.is_empty() {
                candidate.equivalence_members = candidate.label.clone();
            }
            match unique.get_mut(&candidate.sequence) {
                Some(existing) => {
                    let mut members = existing.equivalence_members.clone();
                    append_equivalence_member(&mut members, &candidate.label);
                    if prefer_its2_candidate(existing, &candidate, args.assembly_mode) {
                        candidate.equivalence_members = members;
                        *existing = candidate;
                    } else {
                        existing.equivalence_members = members;
                    }
                }
                None => {
                    unique.insert(candidate.sequence.clone(), candidate);
                }
            }
        }
        let mut equivalence: HashMap<Vec<usize>, ContigRecord> = HashMap::new();
        for mut candidate in unique.into_values() {
            let signature: Vec<usize> = linked_fragments
                .iter()
                .enumerate()
                .filter_map(|(index, fragment)| {
                    fragment
                        .iter()
                        .any(|mate| mate_matches_candidate(mate, &candidate.sequence, slice_len))
                        .then_some(index)
                })
                .collect();
            match equivalence.get_mut(&signature) {
                Some(existing) => {
                    let mut members = existing.equivalence_members.clone();
                    for label in candidate.equivalence_members.split(";") {
                        append_equivalence_member(&mut members, label);
                    }
                    if prefer_its2_candidate(existing, &candidate, args.assembly_mode) {
                        candidate.equivalence_members = members;
                        *existing = candidate;
                    } else {
                        existing.equivalence_members = members;
                    }
                }
                None => {
                    equivalence.insert(signature, candidate);
                }
            }
        }
        let mut candidates: Vec<ContigRecord> = equivalence.into_values().collect();
        annotate_its2_candidates(&mut candidates, &linked_fragments, slice_len);
        accepted = candidates
            .iter()
            .filter(|candidate| candidate.accepted)
            .cloned()
            .collect();
        rejected = candidates
            .into_iter()
            .filter(|candidate| !candidate.accepted)
            .collect();
    }

    let mut low_quality = accepted.is_empty();
    if low_quality && args.assembly_mode == AssemblyMode::Reference {
        accepted = rejected.clone();
        low_quality = accepted.is_empty();
    }
    let pool = if accepted.is_empty() {
        &rejected
    } else {
        &accepted
    };
    let Some(mut best) = pool
        .iter()
        .max_by(|left, right| compare_contigs(left, right, args.assembly_mode))
        .cloned()
    else {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            "Insufficient reads coverage, unable to build contigs.",
        );
        return Ok(LocusResult::failure(&task.key, "no contigs"));
    };
    if paired_fragments {
        if let Some(core_best) = core_accepted
            .iter()
            .max_by(|left, right| compare_contigs(left, right, args.assembly_mode))
        {
            // PE may improve evidence-based ranking, but it must never shorten
            // the final locus relative to its independently assembled core path.
            if best.sequence.len() < core_best.sequence.len() {
                best = core_best.clone();
            }
        }
    }

    accepted.sort_by(|left, right| compare_contigs(right, left, args.assembly_mode));
    rejected.sort_by(|left, right| compare_contigs(right, left, args.assembly_mode));

    if args.assembly_mode == AssemblyMode::Its2 {
        write_contigs(&best_path, &accepted, args.assembly_mode, "its2_variant")?;
        write_contigs(&all_path, &accepted, args.assembly_mode, "its2_variant")?;
        let mut support = accepted.clone();
        support.extend(rejected.clone());
        write_its2_support(
            &args
                .output
                .join("results")
                .join(format!("{}.its2_support.tsv", task.key)),
            &support,
        )?;
    } else if args.assembly_mode == AssemblyMode::Reference || !low_quality {
        write_contigs(
            &best_path,
            std::slice::from_ref(&best),
            args.assembly_mode,
            "contig",
        )?;
        write_contigs(&all_path, &accepted, args.assembly_mode, "contig")?;
    } else {
        remove_if_exists(&best_path);
        remove_if_exists(&all_path);
    }
    if matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) && !rejected.is_empty()
    {
        write_contigs(
            &low_path,
            &rejected,
            args.assembly_mode,
            "low_support_contig",
        )?;
    } else {
        remove_if_exists(&low_path);
    }

    let accepted_locus = args.assembly_mode == AssemblyMode::Reference || !low_quality;
    let length = best.sequence.len();
    if let Some(profile) = profile {
        profile.record_qc(qc_started);
    }
    Ok(LocusResult {
        key: task.key.clone(),
        status: if low_quality {
            "low quality".to_string()
        } else {
            "success".to_string()
        },
        value: best.read_count,
        accepted: accepted_locus,
        rejection_reason: if low_quality {
            best.rejection_reason.clone()
        } else {
            String::new()
        },
        selected_contig_length: length,
        read_supported_span: best.supported_span,
        slice_supported_bases: best.supported_bases,
        slice_support_breadth: best.support_breadth,
        max_slice_support_gap: best.max_support_gap,
        read_count: best.read_count,
        unique_read_count: best.unique_read_count,
        multi_mapping_read_count: best.multi_mapping_read_count,
        read_density: best.read_density,
        unique_read_density: if length > 0 {
            best.unique_read_count as f64 / length as f64
        } else {
            0.0
        },
        support_fraction: best.support_fraction,
        flank_balance: best.flank_balance,
        kmer_median_depth: best.kmer_median_depth,
        kmer_depth_cv: best.kmer_depth_cv,
        kmer_max_depth_ratio: best.kmer_max_depth_ratio,
        candidate_count: accepted.len() + rejected.len(),
        low_quality,
        skipped: false,
    })
}

pub fn read_result_dict(path: &Path) -> io::Result<HashMap<String, (String, u64)>> {
    let mut results = HashMap::new();
    if !path.is_file() {
        return Ok(results);
    }
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 3 && !fields[0].is_empty() {
            results.insert(
                fields[0].to_string(),
                (fields[1].to_string(), fields[2].parse().unwrap_or(0)),
            );
        }
    }
    Ok(results)
}

pub fn write_result_dict(path: &Path, results: &HashMap<String, (String, u64)>) -> io::Result<()> {
    let mut keys: Vec<&String> = results.keys().collect();
    keys.sort();
    let mut writer = BufWriter::new(File::create(path)?);
    for key in keys {
        let (status, value) = &results[key];
        writeln!(writer, "{key},{status},{value},")?;
    }
    writer.flush()
}

pub fn read_summary_lines(path: &Path) -> io::Result<HashMap<String, String>> {
    let mut rows = HashMap::new();
    if !path.is_file() {
        return Ok(rows);
    }
    for (index, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if index == 0 {
            continue;
        }
        let key = line.split(',').next().unwrap_or("");
        if !key.is_empty() {
            rows.insert(key.to_string(), line);
        }
    }
    Ok(rows)
}

fn rounded(value: f64, decimals: usize) -> String {
    let mut text = format!("{value:.decimals$}");
    while text.contains('.') && text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    if text == "-0" {
        "0".to_string()
    } else {
        text
    }
}

fn rounded_float(value: f64, decimals: usize) -> String {
    let mut text = rounded(value, decimals);
    if !text.contains('.') {
        text.push_str(".0");
    }
    text
}

// 汇总行把状态和关键指标压成一条，统计脚本拿来就能用。
pub fn summary_line(result: &LocusResult) -> String {
    format!(
        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        result.key,
        result.status,
        u8::from(result.accepted),
        result.rejection_reason,
        result.selected_contig_length,
        result.read_supported_span,
        result.slice_supported_bases,
        rounded(result.slice_support_breadth, 6),
        result.max_slice_support_gap,
        result.read_count,
        result.unique_read_count,
        result.multi_mapping_read_count,
        rounded(result.read_density, 6),
        rounded(result.unique_read_density, 6),
        rounded(result.support_fraction, 3),
        rounded(result.flank_balance, 3),
        rounded_float(result.kmer_median_depth, 3),
        rounded_float(result.kmer_depth_cv, 3),
        rounded_float(result.kmer_max_depth_ratio, 3),
        result.candidate_count,
        u8::from(result.low_quality),
    )
}

// 汇总按 locus 名排序再写，重复运行也不至于顺序乱飘。
pub fn write_summary(path: &Path, rows: &HashMap<String, String>) -> io::Result<()> {
    let mut keys: Vec<&String> = rows.keys().collect();
    keys.sort();
    let mut writer = BufWriter::new(File::create(path)?);
    write!(writer, "{SUMMARY_HEADER}\r\n")?;
    for key in keys {
        write!(writer, "{}\r\n", rows[key])?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(label: &str, sequence: &[u8]) -> ContigRecord {
        ContigRecord {
            label: label.to_string(),
            equivalence_members: label.to_string(),
            sequence: sequence.to_vec(),
            ..ContigRecord::default()
        }
    }

    #[test]
    fn its2_pe_evidence_assigns_diagnostic_support_and_em() {
        let mut candidates = vec![candidate("A", b"AAAATATA"), candidate("B", b"AAAACGCG")];
        let fragments = vec![
            vec![b"TATA".to_vec(), b"TATA".to_vec()],
            vec![b"CGCG".to_vec()],
            vec![b"AAAA".to_vec()],
        ];
        annotate_its2_candidates(&mut candidates, &fragments, 4);
        assert_eq!(candidates[0].fragment_support, 2);
        assert_eq!(candidates[0].paired_fragment_support, 1);
        assert_eq!(candidates[0].diagnostic_fragment_support, 1);
        assert!(candidates[0].accepted);
        assert_eq!(its2_status(&candidates[0]), "PASS");
        assert_eq!(candidates[1].fragment_support, 2);
        assert_eq!(candidates[1].diagnostic_fragment_support, 1);
        assert!(candidates[1].accepted);
        assert_eq!(its2_status(&candidates[1]), "EQUIVALENCE_GROUP");
        assert!((candidates[0].em_abundance - 0.5).abs() < 1e-6);
        assert!((candidates[1].em_abundance - 0.5).abs() < 1e-6);
    }

    #[test]
    fn its2_labelled_reference_beats_anonymous_graph_candidate() {
        let anonymous = candidate("", b"AAAATATA");
        let labelled = candidate("C40", b"AAAATATA");
        assert!(prefer_its2_candidate(
            &anonymous,
            &labelled,
            AssemblyMode::Its2
        ));
        assert!(!prefer_its2_candidate(
            &labelled,
            &anonymous,
            AssemblyMode::Its2
        ));
        let mut members = "C40".to_string();
        append_equivalence_member(&mut members, "C40e");
        append_equivalence_member(&mut members, "C40");
        assert_eq!(members, "C40;C40e");
    }
}
