use crate::alignment::{align_read, terminal_evidence};
use crate::evidence::{collect_runs_stats, infer_orientation};
use crate::index::{ExactSeed, UceIndex};
use crate::model::{default_spill_path, Candidate, Fragment, FragmentBank, LocusId};
use crate::selection::{choose_auto, choose_legacy, selected};
use gm2_tools::fastx::{FastxFormat, FastxReader, FastxRecord};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Config {
    pub references: PathBuf,
    pub recruit_references: Option<PathBuf>,
    pub read1: PathBuf,
    pub read2: PathBuf,
    pub output: PathBuf,
    pub kmer_size: usize,
    pub step: usize,
    pub min_depth: i64,
    pub max_depth: i64,
    pub max_size_mb: i64,
    pub max_fragments: u64,
    pub memory_limit_mib: u64,
    pub selection_auto: bool,
    pub reference_is_contig: bool,
    pub alignment_shadow: bool,
    pub shadow_per_locus: usize,
    pub shadow_band: usize,
    pub terminal_window: usize,
}

#[derive(Clone, Debug, Default)]
pub struct RunSummary {
    pub fragments_read: u64,
    pub fragments_retained_once: usize,
    pub fragment_bases_retained_once: u64,
    pub assignments: usize,
    pub loci_written: usize,
    pub fragment_memory_bytes: u64,
    pub fragment_spill_bytes: u64,
    pub shadow_sampled_assignments: usize,
    pub shadow_aligned_mates: usize,
    pub shadow_seconds: f64,
    pub elapsed_seconds: f64,
}

const LOCUS_BUFFER_BYTES: usize = 64 * 1024;

struct LocusOutput {
    path: PathBuf,
    buffer: Vec<u8>,
}

struct OutputRouter {
    outputs: Vec<Option<LocusOutput>>,
}

struct ShadowWriter {
    writer: BufWriter<File>,
    band: usize,
    terminal_window: usize,
    sampled_assignments: usize,
    aligned_mates: usize,
    elapsed_seconds: f64,
    loci: Vec<ShadowLocusStats>,
}

#[derive(Clone, Debug, Default)]
struct ShadowLocusStats {
    sampled_assignments: u64,
    aligned_mates: u64,
    both_mates_aligned: u64,
    linked_mate_only: u64,
    multi_locus_assignments: u64,
    identity_sum: f64,
    identity_ge_90: u64,
    overlap_sum: u64,
    near_terminal_assignments: u64,
    extends_left_assignments: u64,
    extends_right_assignments: u64,
    covered_bins: u64,
}

impl ShadowWriter {
    fn new(
        path: &Path,
        band: usize,
        terminal_window: usize,
        locus_count: usize,
    ) -> Result<Self, String> {
        let mut writer = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        writeln!(
            writer,
            "locus\tfragment_ordinal\tmate\tselected_locus_count\tstatus\treference_id\tstrand\treference_length\tterminal_window\tquery_length\tquery_start\tquery_end\treference_start\treference_end\texact_seed_length\tscore\tmatches\tmismatches\tgap_bases\toverlap\tidentity\tleft_overhang\tright_overhang\tnear_left_terminal\tnear_right_terminal\textends_left_terminal\textends_right_terminal"
        )
        .map_err(|e| e.to_string())?;
        Ok(Self {
            writer,
            band,
            terminal_window,
            sampled_assignments: 0,
            aligned_mates: 0,
            elapsed_seconds: 0.0,
            loci: vec![ShadowLocusStats::default(); locus_count],
        })
    }

    fn write_fragment(
        &mut self,
        index: &UceIndex,
        locus: LocusId,
        locus_name: &str,
        fragment: &Fragment,
        selected_locus_count: usize,
    ) -> Result<(), String> {
        self.sampled_assignments += 1;
        let mut mapped_mates = 0_u64;
        let mut near_terminal = false;
        let mut extends_left = false;
        let mut extends_right = false;
        for (mate, record) in [(1, &fragment.r1), (2, &fragment.r2)] {
            let started = Instant::now();
            let alignment = align_read(index, &record.sequence, locus, self.band);
            self.elapsed_seconds += started.elapsed().as_secs_f64();
            let Some(alignment) = alignment else {
                writeln!(
                    self.writer,
                    "{locus_name}\t{}\t{mate}\t{selected_locus_count}\tunmapped\tNA\tNA\tNA\tNA\t{}\tNA\tNA\tNA\tNA\t0\tNA\t0\t0\t0\t0\t0.000000\t0\t0\t0\t0\t0\t0",
                    fragment.ordinal,
                    record.sequence.len(),
                )
                .map_err(|e| e.to_string())?;
                continue;
            };
            self.aligned_mates += 1;
            mapped_mates += 1;
            let reference = &index.references[alignment.sequence as usize];
            let terminal = terminal_evidence(
                alignment,
                record.sequence.len(),
                reference.bases.len(),
                self.terminal_window,
            );
            near_terminal |= terminal.near_left_terminal || terminal.near_right_terminal;
            extends_left |= terminal.extends_left_terminal;
            extends_right |= terminal.extends_right_terminal;
            let stats = &mut self.loci[locus as usize];
            stats.aligned_mates += 1;
            stats.identity_sum += alignment.identity();
            stats.identity_ge_90 += u64::from(alignment.identity() >= 0.90);
            stats.overlap_sum += alignment.reference_overlap() as u64;
            let reference_length = reference.bases.len().max(1);
            let first_bin =
                (terminal.original_reference_start as usize * 64 / reference_length).min(63);
            let last_position = (terminal.original_reference_end as usize).saturating_sub(1);
            let last_bin = (last_position * 64 / reference_length).min(63);
            for bin in first_bin..=last_bin {
                stats.covered_bins |= 1_u64 << bin;
            }
            writeln!(
                self.writer,
                "{locus_name}\t{}\t{mate}\t{selected_locus_count}\taligned\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                fragment.ordinal,
                alignment.sequence,
                alignment.strand,
                reference.bases.len(),
                terminal.effective_window,
                record.sequence.len(),
                alignment.query_start,
                alignment.query_end,
                terminal.original_reference_start,
                terminal.original_reference_end,
                alignment.exact_seed_length,
                alignment.score,
                alignment.matches,
                alignment.mismatches,
                alignment.gap_bases,
                alignment.reference_overlap(),
                alignment.identity(),
                terminal.left_overhang,
                terminal.right_overhang,
                terminal.near_left_terminal as u8,
                terminal.near_right_terminal as u8,
                terminal.extends_left_terminal as u8,
                terminal.extends_right_terminal as u8,
            )
            .map_err(|e| e.to_string())?;
        }
        let stats = &mut self.loci[locus as usize];
        stats.sampled_assignments += 1;
        stats.both_mates_aligned += u64::from(mapped_mates == 2);
        stats.linked_mate_only += u64::from(mapped_mates == 1);
        stats.multi_locus_assignments += u64::from(selected_locus_count > 1);
        stats.near_terminal_assignments += u64::from(near_terminal);
        stats.extends_left_assignments += u64::from(extends_left);
        stats.extends_right_assignments += u64::from(extends_right);
        Ok(())
    }

    fn flush(&mut self) -> Result<(), String> {
        self.writer.flush().map_err(|e| e.to_string())
    }

    fn write_summary(&self, path: &Path, index: &UceIndex) -> Result<(), String> {
        let mut writer = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        writeln!(
            writer,
            "locus\tsampled_assignments\taligned_mates\tboth_mates_aligned\tlinked_mate_only\tmulti_locus_assignments\tmean_identity\tidentity_ge_0.90\tmean_overlap\tnear_terminal_assignments\textends_left_assignments\textends_right_assignments\tcovered_bins_64"
        )
        .map_err(|e| e.to_string())?;
        for (locus, stats) in index.loci.iter().zip(&self.loci) {
            if stats.sampled_assignments == 0 {
                continue;
            }
            let mean_identity = if stats.aligned_mates == 0 {
                0.0
            } else {
                stats.identity_sum / stats.aligned_mates as f64
            };
            let mean_overlap = if stats.aligned_mates == 0 {
                0.0
            } else {
                stats.overlap_sum as f64 / stats.aligned_mates as f64
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
                locus.name,
                stats.sampled_assignments,
                stats.aligned_mates,
                stats.both_mates_aligned,
                stats.linked_mate_only,
                stats.multi_locus_assignments,
                mean_identity,
                stats.identity_ge_90,
                mean_overlap,
                stats.near_terminal_assignments,
                stats.extends_left_assignments,
                stats.extends_right_assignments,
                stats.covered_bins.count_ones(),
            )
            .map_err(|e| e.to_string())?;
        }
        writer.flush().map_err(|e| e.to_string())
    }
}

fn evenly_sample(ids: &[u32], limit: usize) -> Vec<u32> {
    if limit == 0 || ids.len() <= limit {
        return ids.to_vec();
    }
    if limit == 1 {
        return vec![ids[ids.len() / 2]];
    }
    (0..limit)
        .map(|i| ids[i * (ids.len() - 1) / (limit - 1)])
        .collect()
}

impl OutputRouter {
    fn new(locus_count: usize) -> Self {
        Self {
            outputs: (0..locus_count).map(|_| None).collect(),
        }
    }

    fn register(&mut self, locus: usize, path: PathBuf) -> Result<(), String> {
        File::create(&path).map_err(|e| e.to_string())?;
        self.outputs[locus] = Some(LocusOutput {
            path,
            buffer: Vec::with_capacity(LOCUS_BUFFER_BYTES + 1024),
        });
        Ok(())
    }

    fn write_fragment(&mut self, locus: usize, fragment: &Fragment) -> Result<(), String> {
        let should_flush = {
            let output = self.outputs[locus]
                .as_mut()
                .ok_or_else(|| format!("locus {locus} was not registered"))?;
            write_record(&mut output.buffer, &fragment.r1)?;
            write_record(&mut output.buffer, &fragment.r2)?;
            output.buffer.len() >= LOCUS_BUFFER_BYTES
        };
        if should_flush {
            self.flush_one(locus)?;
        }
        Ok(())
    }

    fn flush_one(&mut self, locus: usize) -> Result<(), String> {
        let Some(output) = self.outputs[locus].as_mut() else {
            return Ok(());
        };
        if output.buffer.is_empty() {
            return Ok(());
        }
        let mut file = OpenOptions::new()
            .append(true)
            .open(&output.path)
            .map_err(|e| e.to_string())?;
        file.write_all(&output.buffer).map_err(|e| e.to_string())?;
        output.buffer.clear();
        Ok(())
    }

    fn flush(&mut self) -> Result<(), String> {
        for locus in 0..self.outputs.len() {
            self.flush_one(locus)?;
        }
        Ok(())
    }
}

fn write_record(out: &mut impl Write, record: &FastxRecord) -> Result<(), String> {
    out.write_all(&record.header).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    out.write_all(&record.sequence).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    out.write_all(if record.plus.is_empty() {
        b"+\n"
    } else {
        &record.plus
    })
    .map_err(|e| e.to_string())?;
    if !record.plus.is_empty() {
        out.write_all(b"\n").map_err(|e| e.to_string())?;
    }
    out.write_all(&record.quality).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

fn next_pair(
    r1: &mut FastxReader,
    r2: &mut FastxReader,
) -> Result<Option<(FastxRecord, FastxRecord)>, String> {
    let first = r1.next_record().map_err(|e| e.to_string())?;
    let second = r2.next_record().map_err(|e| e.to_string())?;
    match (first, second) {
        (None, None) => Ok(None),
        (Some(a), Some(b)) => Ok(Some((a, b))),
        _ => Err("paired input files contain different numbers of records".to_string()),
    }
}

#[inline]
fn keep_linked_pair(orient1: u8, orient2: u8) -> bool {
    !((1..=2).contains(&orient1) && orient1 == orient2 || orient1 == 0 && orient2 == 0)
}

#[derive(Clone, Copy, Debug, Default)]
struct FastEvidence {
    max_exact: u16,
    covered_bins: u64,
    terminal_mask: u8,
    left_extension: u16,
    right_extension: u16,
    aligned_mates: u8,
}

fn add_exact_evidence(
    index: &UceIndex,
    seed: ExactSeed,
    query_length: usize,
    terminal_window: usize,
    evidence: &mut FastEvidence,
) {
    evidence.max_exact = evidence
        .max_exact
        .max(seed.len().min(u16::MAX as usize) as u16);
    evidence.aligned_mates = evidence.aligned_mates.saturating_add(1);
    let reference = &index.references[seed.sequence as usize];
    let reference_length = reference.bases.len().max(1);
    let (original_start, original_end, left_overhang, right_overhang) = if reference.strand == 1 {
        (
            seed.reference_start as usize,
            seed.reference_end as usize,
            seed.read_start as usize,
            query_length.saturating_sub(seed.read_end as usize),
        )
    } else {
        (
            reference_length.saturating_sub(seed.reference_end as usize),
            reference_length.saturating_sub(seed.reference_start as usize),
            query_length.saturating_sub(seed.read_end as usize),
            seed.read_start as usize,
        )
    };
    let first_bin = (original_start * 64 / reference_length).min(63);
    let last_position = original_end.saturating_sub(1);
    let last_bin = (last_position * 64 / reference_length).min(63);
    for bin in first_bin..=last_bin {
        evidence.covered_bins |= 1_u64 << bin;
    }
    let effective_window = terminal_window.min(reference_length / 5);
    if original_start <= effective_window && left_overhang > 0 {
        evidence.terminal_mask |= 1;
        evidence.left_extension = evidence
            .left_extension
            .max(left_overhang.min(u16::MAX as usize) as u16);
    }
    if reference_length.saturating_sub(original_end) <= effective_window && right_overhang > 0 {
        evidence.terminal_mask |= 2;
        evidence.right_extension = evidence
            .right_extension
            .max(right_overhang.min(u16::MAX as usize) as u16);
    }
}

pub fn run(config: &Config) -> Result<RunSummary, String> {
    let started = Instant::now();
    fs::create_dir_all(&config.output).map_err(|e| e.to_string())?;
    let filtered = config.output.join("filtered");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&filtered).map_err(|e| e.to_string())?;
    let index_started = Instant::now();
    let index = UceIndex::build_split(
        config
            .recruit_references
            .as_deref()
            .unwrap_or(&config.references),
        &config.references,
        config.kmer_size,
    )?;
    eprintln!(
        "UCEFilter index: {} loci, {} positional anchors, k={}, run-k={} ({:.3}s)",
        index.loci.len(),
        index.anchor_entries(),
        index.k,
        index.run_k,
        index_started.elapsed().as_secs_f64()
    );
    let mut reader1 =
        FastxReader::open(&config.read1, FastxFormat::Fastq).map_err(|e| e.to_string())?;
    let mut reader2 =
        FastxReader::open(&config.read2, FastxFormat::Fastq).map_err(|e| e.to_string())?;
    let memory_limit = config.memory_limit_mib.saturating_mul(1024 * 1024);
    let mut bank = FragmentBank::new(memory_limit, default_spill_path(&config.output));
    let mut locus_candidates: Vec<Vec<Candidate>> = vec![Vec::new(); index.loci.len()];
    let mut coarse_counts = vec![0_u64; index.loci.len()];
    let mut ordinal = 0_u64;
    while let Some((r1, r2)) = next_pair(&mut reader1, &mut reader2)? {
        if config.max_fragments > 0 && ordinal >= config.max_fragments {
            break;
        }
        let mut loci = Vec::<LocusId>::new();
        index.recruit(&r1.sequence, config.step, &mut loci);
        index.recruit(&r2.sequence, config.step, &mut loci);
        loci.sort_unstable();
        for &locus in &loci {
            coarse_counts[locus as usize] += 2;
        }
        if loci.is_empty() {
            ordinal += 1;
            continue;
        }
        let events1 = index.orientation_events(&r1.sequence, &loci);
        let events2 = index.orientation_events(&r2.sequence, &loci);
        let mut evidence = Vec::new();
        for (i, &locus) in loci.iter().enumerate() {
            let orient1 = infer_orientation(collect_runs_stats(&events1[i]));
            let orient2 = infer_orientation(collect_runs_stats(&events2[i]));
            if !keep_linked_pair(orient1, orient2) {
                continue;
            }
            let mut fast = FastEvidence::default();
            if let Some(seed) = index.best_exact(&r1.sequence, locus) {
                add_exact_evidence(
                    &index,
                    seed,
                    r1.sequence.len(),
                    config.terminal_window,
                    &mut fast,
                );
            }
            if let Some(seed) = index.best_exact(&r2.sequence, locus) {
                add_exact_evidence(
                    &index,
                    seed,
                    r2.sequence.len(),
                    config.terminal_window,
                    &mut fast,
                );
            }
            evidence.push((locus, fast));
        }
        if !evidence.is_empty() {
            let fragment_bases = (r1.sequence.len() + r2.sequence.len()) as u32;
            let fragment_id = bank.insert(Fragment { ordinal, r1, r2 })?;
            let locus_count = evidence.len().min(u16::MAX as usize) as u16;
            for (locus, fast) in evidence {
                locus_candidates[locus as usize].push(Candidate {
                    fragment_id,
                    ordinal,
                    fragment_bases,
                    max_exact: fast.max_exact,
                    covered_bins: fast.covered_bins,
                    terminal_mask: fast.terminal_mask,
                    left_extension: fast.left_extension,
                    right_extension: fast.right_extension,
                    aligned_mates: fast.aligned_mates,
                    locus_count,
                });
            }
        }
        ordinal += 1;
        if ordinal.is_multiple_of(1_048_576) {
            eprintln!("UCEFilter handled {} Mi fragments", ordinal / 1_048_576);
        }
    }
    let mut loci_written = 0_usize;
    let mut assignments = 0_usize;
    let mut routes: Vec<Vec<LocusId>> = (0..bank.len()).map(|_| Vec::new()).collect();
    let mut shadow_routes = config.alignment_shadow.then(|| {
        (0..bank.len())
            .map(|_| Vec::new())
            .collect::<Vec<Vec<LocusId>>>()
    });
    let mut router = OutputRouter::new(index.loci.len());
    let mut summary = BufWriter::new(
        File::create(config.output.join("uce_filter_summary.tsv")).map_err(|e| e.to_string())?,
    );
    if config.selection_auto {
        writeln!(
            summary,
            "locus\tcoarse_reads\trun_fragments\tselected_fragments\tminimum_exact\tthinning_interval\tselection_mode\teligible_fragments\tcovered_bins_64\ttarget_core_fragments\tleft_terminal_candidates\tright_terminal_candidates\tselected_left_terminal\tselected_right_terminal"
        )
        .map_err(|e| e.to_string())?;
    } else {
        writeln!(
            summary,
            "locus\tcoarse_reads\trun_fragments\tselected_fragments\tminimum_exact\tthinning_interval"
        )
        .map_err(|e| e.to_string())?;
    }
    for (locus_id, candidates) in locus_candidates.iter().enumerate() {
        let locus = &index.loci[locus_id];
        let decision = choose_legacy(
            candidates,
            locus.effective_length,
            config.kmer_size,
            config.min_depth,
            config.max_depth,
            config.max_size_mb,
        );
        let (selected_ids, auto_details) = if config.selection_auto {
            let automatic = choose_auto(
                candidates,
                &decision,
                locus.effective_length,
                config.reference_is_contig,
            );
            let details = (
                automatic.mode,
                automatic.eligible_fragments,
                automatic.covered_bins,
                automatic.target_core_fragments,
                automatic.left_terminal_candidates,
                automatic.right_terminal_candidates,
                automatic.selected_left_terminal,
                automatic.selected_right_terminal,
            );
            (automatic.selected_ids, Some(details))
        } else {
            let mut passing = 0_usize;
            let ids = candidates
                .iter()
                .filter(|candidate| selected(candidate, &decision, &mut passing))
                .map(|candidate| candidate.fragment_id)
                .collect();
            (ids, None)
        };
        let minimum_exact = decision
            .minimum_exact
            .map_or_else(|| "all".to_string(), |v| v.to_string());
        let thinning_interval = decision
            .thinning_interval
            .map_or_else(|| "none".to_string(), |v| v.to_string());
        if let Some((
            mode,
            eligible,
            covered_bins,
            target_core,
            left_terminal_candidates,
            right_terminal_candidates,
            selected_left_terminal,
            selected_right_terminal,
        )) = auto_details
        {
            writeln!(
                summary,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                locus.name,
                coarse_counts[locus_id],
                candidates.len(),
                selected_ids.len(),
                minimum_exact,
                thinning_interval,
                mode.as_str(),
                eligible,
                covered_bins,
                target_core,
                left_terminal_candidates,
                right_terminal_candidates,
                selected_left_terminal,
                selected_right_terminal,
            )
            .map_err(|e| e.to_string())?;
        } else {
            writeln!(
                summary,
                "{}\t{}\t{}\t{}\t{}\t{}",
                locus.name,
                coarse_counts[locus_id],
                candidates.len(),
                selected_ids.len(),
                minimum_exact,
                thinning_interval,
            )
            .map_err(|e| e.to_string())?;
        }
        if selected_ids.is_empty() {
            continue;
        }
        loci_written += 1;
        assignments += selected_ids.len();
        router.register(locus_id, filtered.join(format!("{}.fq", locus.name)))?;
        if let Some(shadow_routes) = shadow_routes.as_mut() {
            for id in evenly_sample(&selected_ids, config.shadow_per_locus) {
                shadow_routes[id as usize].push(locus_id as LocusId);
            }
        }
        for id in selected_ids {
            routes[id as usize].push(locus_id as LocusId);
        }
    }
    summary.flush().map_err(|e| e.to_string())?;
    let fragment_memory_bytes = bank.memory_bytes();
    let fragment_spill_bytes = bank.spill_bytes();
    let mut shadow_writer = if config.alignment_shadow {
        Some(ShadowWriter::new(
            &config.output.join("alignment_shadow.tsv"),
            config.shadow_band,
            config.terminal_window,
            index.loci.len(),
        )?)
    } else {
        None
    };
    bank.stream_in_order(|id, fragment| {
        for &locus in &routes[id as usize] {
            router.write_fragment(locus as usize, fragment)?;
        }
        if let (Some(shadow_routes), Some(shadow_writer)) =
            (shadow_routes.as_ref(), shadow_writer.as_mut())
        {
            for &locus in &shadow_routes[id as usize] {
                shadow_writer.write_fragment(
                    &index,
                    locus,
                    &index.loci[locus as usize].name,
                    fragment,
                    routes[id as usize].len(),
                )?;
            }
        }
        Ok(())
    })?;
    router.flush()?;
    if let Some(shadow_writer) = shadow_writer.as_mut() {
        shadow_writer.flush()?;
        shadow_writer.write_summary(&config.output.join("alignment_shadow_summary.tsv"), &index)?;
    }
    let mut counts = BufWriter::new(
        File::create(config.output.join("ref_reads_count_dict.txt")).map_err(|e| e.to_string())?,
    );
    for (locus, count) in index.loci.iter().zip(coarse_counts) {
        if count > 0 {
            writeln!(counts, "{},{}", locus.name, count).map_err(|e| e.to_string())?;
        }
    }
    counts.flush().map_err(|e| e.to_string())?;
    let (shadow_sampled_assignments, shadow_aligned_mates, shadow_seconds) =
        shadow_writer.map_or((0, 0, 0.0), |writer| {
            (
                writer.sampled_assignments,
                writer.aligned_mates,
                writer.elapsed_seconds,
            )
        });
    Ok(RunSummary {
        fragments_read: ordinal,
        fragments_retained_once: bank.len(),
        fragment_bases_retained_once: bank.bases(),
        assignments,
        loci_written,
        fragment_memory_bytes,
        fragment_spill_bytes,
        shadow_sampled_assignments,
        shadow_aligned_mates,
        shadow_seconds,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

pub fn output_exists(path: &Path) -> bool {
    path.join("filtered").is_dir() && path.join("ref_reads_count_dict.txt").is_file()
}

#[cfg(test)]
mod tests {
    use super::{add_exact_evidence, evenly_sample, keep_linked_pair, FastEvidence};
    use crate::index::UceIndex;
    use std::fs;
    use std::io::Write;

    #[test]
    fn paired_end_decision_preserves_legacy_linked_mate_semantics() {
        assert!(!keep_linked_pair(0, 0));
        assert!(keep_linked_pair(1, 0));
        assert!(keep_linked_pair(0, 2));
        assert!(!keep_linked_pair(1, 1));
        assert!(!keep_linked_pair(2, 2));
        assert!(keep_linked_pair(1, 2));
        assert!(keep_linked_pair(3, 0));
        assert!(keep_linked_pair(3, 3));
    }

    #[test]
    fn bounded_shadow_sampling_is_deterministic_and_spans_input_order() {
        let ids: Vec<u32> = (0..100).collect();
        assert_eq!(evenly_sample(&ids, 4), vec![0, 33, 66, 99]);
        assert_eq!(evenly_sample(&ids[..3], 4), vec![0, 1, 2]);
        assert_eq!(evenly_sample(&ids, 1), vec![50]);
    }

    #[test]
    fn exact_evidence_normalizes_reverse_coordinates_and_terminal_sides() {
        let root = std::env::temp_dir().join(format!("uce-filter-pipeline-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let sequence = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let mut out = fs::File::create(root.join("locus.fa")).unwrap();
        writeln!(out, ">ref\n{}", String::from_utf8_lossy(sequence)).unwrap();
        drop(out);
        let index = UceIndex::build(&root, 16).unwrap();
        let read = crate::index::reverse_complement(&sequence[..32]);
        let seed = index.best_exact(&read, 0).unwrap();
        let mut evidence = FastEvidence::default();
        add_exact_evidence(&index, seed, read.len() + 5, 150, &mut evidence);
        assert_ne!(evidence.covered_bins, 0);
        assert_eq!(evidence.aligned_mates, 1);
        assert_ne!(evidence.terminal_mask, 0);
        fs::remove_dir_all(root).unwrap();
    }
}
