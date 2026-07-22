//! Native Rust command dispatcher.
//!
//! The public CLI is implemented in Rust and does not require a Python runtime.

use sha2::{Digest, Sha256};
use std::env;
use std::fs;
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::sync::{Arc, Condvar, Mutex};
use std::thread;
use std::time::Instant;

const COMMANDS: &[&str] = &[
    "filter",
    "refilter",
    "assemble",
    "gene",
    "stats",
    "te",
    "population",
    "consensus",
    "trim",
    "combine",
    "tree",
    "gene-annotate",
    "gene-resolve",
    "gene-tree",
    "profiling",
    "mito",
    "rad",
    "rad-probe",
    "rad-validate",
];
const FLAG_OPTIONS: &[&str] = &[
    "--uce-alignment-shadow",
    "--uce-rescue-reads",
    "--stats-count-input-reads",
    "--stats-no-heatmap",
    "--population-panrefv2-include-low-confidence",
    "--population-skip-mark-duplicates",
    "--population-skip-plink",
    "--population-skip-admixture",
    "--strict-combine-errors",
    "--no-alignment",
    "--no-trimal",
    "--profile-force-rebuild",
    "--cleanup-intermediates",
    "--no-mito-adaptive-stop",
    "--reuse-reference-cache",
    "--legacy-uce-filter",
    "--workflow-profile",
    "--rad-denovo",
];
const VALUE_OPTIONS: &[&str] = &[
    "-f",
    "-r",
    "-o",
    "-p",
    "--assembly-mode",
    "-kf",
    "-s",
    "--step-size",
    "-ka",
    "--min-ka",
    "--max-ka",
    "-e",
    "--error-threshold",
    "-i",
    "--search-depth",
    "-sb",
    "--soft-boundary",
    "--min-coverage",
    "--depth-low-water-mark",
    "--depth-limit",
    "--file-size-limit",
    "--max-reads",
    "--assembler-graph-format",
    "--uce-side-candidates",
    "--uce-max-contig-length",
    "--uce-min-read-density",
    "--uce-density-check-min-length",
    "--uce-max-depth-cv",
    "--uce-max-depth-ratio",
    "--uce-shadow-per-locus",
    "--uce-shadow-band",
    "--uce-shadow-terminal-window",
    "--te-stage",
    "--te-kmer",
    "--te-min-kmer-count",
    "--te-catalog-pairs",
    "--te-read-ledger",
    "--te-library",
    "--te-annotate-min-fragment",
    "--te-annotate-max-fragment",
    "--te-annotate-min-support",
    "--te-annotate-min-identity",
    "--te-annotate-min-coverage",
    "--te-annotate-min-delta",
    "--te-assemble-min-kmer-count",
    "--te-assemble-branch-ratio",
    "--te-assemble-max-fragments",
    "--engine",
    "--population-reference-strategy",
    "--population-reference-fasta",
    "--population-min-mapq",
    "--population-min-baseq",
    "--population-min-dp",
    "--population-min-gq",
    "--population-min-qual",
    "--population-min-call-rate",
    "--population-min-mac",
    "--population-ld-window",
    "--population-ld-step",
    "--population-ld-r2",
    "--population-admixture-k-min",
    "--population-admixture-k-max",
    "--population-admixture-cv",
    "--population-start-at",
    "--population-stop-after",
    "--population-minibwa",
    "--population-samtools",
    "--population-bcftools",
    "--population-plink",
    "--population-admixture",
    "-c",
    "--consensus-threshold",
    "-ts",
    "--trim-source",
    "-tm",
    "--trim-mode",
    "-tr",
    "--trim-retention",
    "-cs",
    "--combine-source",
    "-cd",
    "--clean-difference",
    "-cn",
    "--clean-sequences",
    "--msa-program",
    "--msa-threads",
    "--alignment-filter",
    "--filter-processes",
    "--alifilter-model",
    "-m",
    "--tree-method",
    "-b",
    "--bootstrap",
    "--phylo-program",
    "--gene-protein-reference",
    "--gene-miniprot",
    "--gene-input",
    "--gene-mafft",
    "--gene-iqtree",
    "--gene-min-taxa",
    "--gene-min-aa-length",
    "--gene-min-effective-codon-sites",
    "--gene-outgroup",
    "--gene-ufboot",
    "--gene-taper",
    "--gene-julia",
    "--gene-species-mode",
    "--gene-aster",
    "--profile-kmer-size",
    "--profile-pseudoalign-threshold",
    "--profile-relevant-kmer-fraction",
    "--profile-group-map",
    "--profile-decoy",
    "--profile-index-dir",
    "--profile-index-memory-gb",
    "--profile-themisto",
    "--reference-cache-dir",
    "--mito-genbank",
    "--mito-flank",
    "--mito-tile-length",
    "--mito-tile-step",
    "--mito-min-overlap",
    "--mito-min-overlap-identity",
    "--mito-min-junction-support",
    "--mito-terminal-window",
    "--mito-link-kmer",
    "--mito-min-link-hits",
    "--mito-min-pair-support",
    "--mito-bridge-kmer",
    "--mito-bridge-min-depth",
    "--mito-max-bridge",
    "--mito-initial-reads",
    "--mito-max-reads",
    "--uce-rescue-rounds",
    "--uce-rescue-min-contig-length",
    "--uce-rescue-terminal-window",
    "--uce-rescue-min-density-ratio",
    "--assembler-implementation",
    "--assembler-read-chunk-size",
    "--uce-path-strategy",
    "--uce-backbone-lookahead",
    "--min-depth",
    "--max-depth",
    "--ipyrad-loci",
    "--rad-min-arm-breadth",
    "--rad-probe",
    "--ipyrad-params",
    "--ipyrad-executable",
    "--ipyrad-steps",
    "--rad-overhang",
    "--rad-overhang-r2",
    "--rad-kmer",
    "--rad-min-count",
    "--rad-min-samples",
    "--rad-min-length",
    "--rad-recovery",
    "--rad-validate-min-identity",
    "--rad-validate-min-breadth",
    "--rad-validate-min-delta",
];

#[derive(Clone, Debug)]
struct Sample {
    name: String,
    read1: String,
    read2: Option<String>,
}

#[derive(Clone, Debug)]
struct Options {
    raw: Vec<String>,
    commands: Vec<String>,
    reference: String,
    samples: String,
    output: String,
    assembly_mode: String,
    workers: usize,
    kf: String,
    step: String,
    ka: String,
    min_ka: String,
    max_ka: String,
    error_threshold: String,
    search_depth: String,
    soft_boundary: String,
    min_coverage: String,
    low_depth: String,
    depth_limit: String,
    size_limit: String,
    max_reads: String,
    graph_format: String,
    side_candidates: String,
    max_contig_length: String,
    min_density: String,
    density_min_length: String,
    max_depth_cv: String,
    max_depth_ratio: String,
    alignment_shadow: bool,
    shadow_per_locus: String,
    shadow_band: String,
    shadow_terminal_window: String,
    rescue: bool,
    stats_count_input_reads: bool,
    stats_no_heatmap: bool,
    cleanup_intermediates: bool,
    reuse_reference_cache: bool,
    legacy_uce_filter: bool,
    workflow_profile: bool,
}

fn value(args: &[String], names: &[&str], default: &str) -> Result<String, String> {
    for (index, arg) in args.iter().enumerate() {
        if names.contains(&arg.as_str()) {
            return args
                .get(index + 1)
                .cloned()
                .ok_or_else(|| format!("{arg} requires a value"));
        }
        for name in names {
            if let Some(value) = arg.strip_prefix(&format!("{name}=")) {
                return Ok(value.to_owned());
            }
        }
    }
    Ok(default.to_owned())
}

fn flag(args: &[String], name: &str) -> Result<bool, String> {
    for arg in args {
        if arg == name {
            return Ok(true);
        }
        if arg
            .strip_prefix(name)
            .is_some_and(|suffix| suffix.starts_with('='))
        {
            return Err(format!("{name} does not take a value"));
        }
    }
    Ok(false)
}

fn commands(args: &[String]) -> Result<Vec<String>, String> {
    let mut commands = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let arg = &args[index];
        if COMMANDS.contains(&arg.as_str()) {
            commands.push(arg.clone());
            index += 1;
        } else if arg.starts_with('-') {
            let option = arg.split_once('=').map_or(arg.as_str(), |(name, _)| name);
            if FLAG_OPTIONS.contains(&option) {
                if arg.contains('=') {
                    return Err(format!("{option} does not take a value"));
                }
                index += 1;
            } else if VALUE_OPTIONS.contains(&option) {
                if arg.contains('=') {
                    index += 1;
                } else {
                    if index + 1 == args.len() {
                        return Err(format!("{arg} requires a value"));
                    }
                    index += 2;
                }
            } else {
                return Err(format!("Rust CLI does not support option '{option}'"));
            }
        } else {
            return Err(format!("Rust CLI does not support command '{arg}'"));
        }
    }
    Ok(commands)
}

fn parse(args: &[String]) -> Result<Options, String> {
    let mut commands = commands(args)?;
    let assembly_mode = value(args, &["--assembly-mode"], "original")?;
    if !matches!(assembly_mode.as_str(), "original" | "uce") {
        return Err("--assembly-mode must be original or uce".into());
    }
    if commands == ["gene"] {
        commands = vec![
            "filter".into(),
            "refilter".into(),
            "assemble".into(),
            "gene".into(),
        ];
    } else if commands.is_empty() {
        commands = if assembly_mode == "uce" {
            vec![
                "filter".into(),
                "refilter".into(),
                "assemble".into(),
                "combine".into(),
                "tree".into(),
            ]
        } else {
            vec![
                "filter".into(),
                "refilter".into(),
                "assemble".into(),
                "trim".into(),
                "combine".into(),
                "tree".into(),
            ]
        };
    }
    if commands.iter().any(|command| command == "gene") {
        let missing = ["filter", "refilter", "assemble"]
            .into_iter()
            .filter(|stage| !commands.iter().any(|command| command == stage))
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            return Err(format!(
                "gene requires filter, refilter, and assemble; missing {}",
                missing.join(", ")
            ));
        }
    }
    Ok(Options {
        raw: args.to_vec(),
        commands,
        reference: value(args, &["-r"], "")?,
        samples: value(args, &["-f"], "")?,
        output: value(args, &["-o"], "")?,
        assembly_mode,
        workers: value(args, &["-p"], "1")?
            .parse()
            .map_err(|_| "-p must be a positive integer")?,
        kf: value(args, &["-kf"], "31")?,
        step: value(args, &["-s", "--step-size"], "4")?,
        ka: value(args, &["-ka"], "0")?,
        min_ka: value(args, &["--min-ka"], "21")?,
        max_ka: value(args, &["--max-ka"], "51")?,
        error_threshold: value(args, &["-e", "--error-threshold"], "2")?,
        search_depth: value(args, &["-i", "--search-depth"], "4096")?,
        soft_boundary: value(args, &["-sb", "--soft-boundary"], "auto")?,
        min_coverage: value(args, &["--min-coverage"], "0")?,
        low_depth: value(args, &["--depth-low-water-mark"], "50")?,
        depth_limit: value(args, &["--depth-limit"], "768")?,
        size_limit: value(args, &["--file-size-limit"], "6")?,
        max_reads: value(args, &["--max-reads"], "0")?,
        graph_format: value(args, &["--assembler-graph-format"], "none")?,
        side_candidates: value(args, &["--uce-side-candidates"], "8")?,
        max_contig_length: value(args, &["--uce-max-contig-length"], "0")?,
        min_density: value(args, &["--uce-min-read-density"], "0.003")?,
        density_min_length: value(args, &["--uce-density-check-min-length"], "1000")?,
        max_depth_cv: value(args, &["--uce-max-depth-cv"], "0")?,
        max_depth_ratio: value(args, &["--uce-max-depth-ratio"], "0")?,
        alignment_shadow: flag(args, "--uce-alignment-shadow")?,
        shadow_per_locus: value(args, &["--uce-shadow-per-locus"], "64")?,
        shadow_band: value(args, &["--uce-shadow-band"], "32")?,
        shadow_terminal_window: value(args, &["--uce-shadow-terminal-window"], "150")?,
        rescue: flag(args, "--uce-rescue-reads")?,
        stats_count_input_reads: flag(args, "--stats-count-input-reads")?,
        stats_no_heatmap: flag(args, "--stats-no-heatmap")?,
        cleanup_intermediates: flag(args, "--cleanup-intermediates")?,
        reuse_reference_cache: flag(args, "--reuse-reference-cache")?,
        legacy_uce_filter: flag(args, "--legacy-uce-filter")?,
        workflow_profile: flag(args, "--workflow-profile")?,
    })
}

fn sample_name(raw: &str) -> String {
    let value: String = raw
        .trim()
        .chars()
        .filter_map(|c| match c {
            ' ' | '-' => Some('_'),
            c if c.is_alphanumeric() || c == '_' || c == '.' => Some(c),
            _ => None,
        })
        .collect();
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return value;
    };
    let mut normalized = first.to_ascii_uppercase().to_string();
    normalized.extend(chars.map(|c| c.to_ascii_lowercase()));
    normalized
}

fn read_samples(path: &str, output: &Path) -> Result<Vec<Sample>, String> {
    let file =
        fs::File::open(path).map_err(|e| format!("Unable to read sample list '{path}': {e}"))?;
    let mut samples = Vec::new();
    for (index, line) in io::BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 2 {
            return Err(format!("Sample '{}' has no data files", fields[0]));
        }
        let name = sample_name(fields[0]);
        if name.is_empty() {
            return Err(format!("Invalid sample name '{}'", fields[0]));
        }
        let numbered = format!("{}_{}", index + 1, name);
        fs::create_dir_all(output.join(&numbered)).map_err(|e| e.to_string())?;
        // Keep the legacy two-column convention: a single supplied FASTX path is
        // deliberately used for both mates, rather than silently changing mode.
        let read2 = fields
            .get(2)
            .filter(|value| !value.is_empty())
            .map(|value| (*value).to_string())
            .or_else(|| Some(fields[1].to_string()));
        samples.push(Sample {
            name: numbered,
            read1: fields[1].to_string(),
            read2,
        });
    }
    if samples.is_empty() {
        return Err("Sample list is empty or invalid".into());
    }
    Ok(samples)
}

fn read_rad_samples(path: &str) -> Result<Vec<Sample>, String> {
    let file = fs::File::open(path)
        .map_err(|e| format!("Unable to read RAD sample list '{path}': {e}"))?;
    let mut samples = Vec::new();
    let mut names = std::collections::BTreeSet::new();
    for (index, line) in io::BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 || fields[1].is_empty() || fields[2].is_empty() {
            return Err(format!(
                "RAD sample row {} must be: sample<TAB>R1.fastq<TAB>R2.fastq",
                index + 1
            ));
        }
        let name = sample_name(fields[0]);
        if name.is_empty() {
            return Err(format!("Invalid RAD sample name '{}'", fields[0]));
        }
        if !names.insert(name.clone()) {
            return Err(format!(
                "Duplicate RAD sample name after normalization: {name}"
            ));
        }
        for read in [fields[1], fields[2]] {
            if !Path::new(read).is_file() {
                return Err(format!("RAD read file does not exist: {read}"));
            }
        }
        samples.push(Sample {
            name,
            read1: fields[1].to_owned(),
            read2: Some(fields[2].to_owned()),
        });
    }
    if samples.is_empty() {
        return Err("RAD sample list has no paired reads".into());
    }
    Ok(samples)
}

fn components() -> Result<PathBuf, String> {
    if let Some(path) = env::var_os("GM2_COMPONENT_DIR") {
        return Ok(PathBuf::from(path));
    }
    let exe = env::current_exe().map_err(|e| e.to_string())?;
    Ok(exe
        .parent()
        .ok_or("Cannot locate GeneMiner components")?
        .to_path_buf())
}

fn run(binary_dir: &Path, name: &str, args: &[String]) -> Result<(), String> {
    let program = binary_dir.join(name);
    let status = Command::new(&program)
        .args(args)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn soft_boundary(value: &str) -> Result<String, String> {
    match value {
        "auto" => Ok("-1".into()),
        "unlimited" => Ok("10000".into()),
        other => other
            .parse::<i32>()
            .map_err(|_| "invalid --soft-boundary".into())
            .map(|v| v.to_string()),
    }
}

fn uce_filter_args_for_recruit(
    opt: &Options,
    sample: &Sample,
    sample_dir: &Path,
    verify_reference: &Path,
    recruit_reference: &Path,
    role: &str,
) -> Vec<String> {
    let mut args = vec![
        "-r".into(),
        verify_reference.display().to_string(),
        "--recruit-references".into(),
        recruit_reference.display().to_string(),
        "-q1".into(),
        sample.read1.clone(),
    ];
    if let Some(read2) = &sample.read2 {
        args.extend(["-q2".into(), read2.clone()]);
    }
    args.extend([
        "-o".into(),
        sample_dir.display().to_string(),
        "-kf".into(),
        opt.kf.clone(),
        "-s".into(),
        opt.step.clone(),
        "--selection".into(),
        "auto".into(),
        "--reference-role".into(),
        role.into(),
        "--threads".into(),
        "1".into(),
        "--memory-limit-mib".into(),
        "256".into(),
        "--min-depth".into(),
        opt.low_depth.clone(),
        "--max-depth".into(),
        opt.depth_limit.clone(),
        "--max-size".into(),
        opt.size_limit.clone(),
    ]);
    if opt.max_reads != "0" {
        args.extend(["--max-fragments".into(), opt.max_reads.clone()]);
    }
    if opt.alignment_shadow {
        args.extend([
            "--alignment-shadow".into(),
            "--shadow-per-locus".into(),
            opt.shadow_per_locus.clone(),
            "--shadow-band".into(),
            opt.shadow_band.clone(),
            "--terminal-window".into(),
            opt.shadow_terminal_window.clone(),
        ]);
    }
    args
}

fn uce_filter_args_for(
    opt: &Options,
    sample: &Sample,
    sample_dir: &Path,
    reference: &Path,
    role: &str,
) -> Vec<String> {
    uce_filter_args_for_recruit(opt, sample, sample_dir, reference, reference, role)
}

fn uce_filter_args(opt: &Options, sample: &Sample, sample_dir: &Path) -> Vec<String> {
    uce_filter_args_for(opt, sample, sample_dir, Path::new(&opt.reference), "bait")
}

fn uce_assembler_args(opt: &Options, sample_dir: &Path) -> Result<Vec<String>, String> {
    Ok(vec![
        "-r".into(),
        opt.reference.clone(),
        "-o".into(),
        sample_dir.display().to_string(),
        "-ka".into(),
        opt.ka.clone(),
        "-k_min".into(),
        opt.min_ka.clone(),
        "-k_max".into(),
        opt.max_ka.clone(),
        "-limit_count".into(),
        opt.error_threshold.clone(),
        "-iteration".into(),
        opt.search_depth.clone(),
        "-sb".into(),
        soft_boundary(&opt.soft_boundary)?,
        "-cov_min".into(),
        opt.min_coverage.clone(),
        "-p".into(),
        "1".into(),
        "--assembly-mode".into(),
        "uce".into(),
        "--uce-side-candidates".into(),
        opt.side_candidates.clone(),
        "--uce-max-contig-length".into(),
        opt.max_contig_length.clone(),
        "--uce-min-read-density".into(),
        opt.min_density.clone(),
        "--uce-density-check-min-length".into(),
        opt.density_min_length.clone(),
        "--uce-max-depth-cv".into(),
        opt.max_depth_cv.clone(),
        "--uce-max-depth-ratio".into(),
        opt.max_depth_ratio.clone(),
        "--uce-path-strategy".into(),
        value(&opt.raw, &["--uce-path-strategy"], "backbone")?,
        "--uce-backbone-lookahead".into(),
        value(&opt.raw, &["--uce-backbone-lookahead"], "24")?,
        "--assembler-read-chunk-size".into(),
        value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
        "--assembler-kmer-count-threads".into(),
        "1".into(),
        "--assembler-graph-format".into(),
        opt.graph_format.clone(),
    ])
}

fn build_uce_rescue_reference(
    reference: &Path,
    sample: &Path,
    rescue: &Path,
    minimum: usize,
    active: Option<&std::collections::BTreeSet<String>>,
) -> Result<usize, String> {
    if rescue.exists() {
        fs::remove_dir_all(rescue).map_err(|e| e.to_string())?;
    }
    copy_tree(reference, rescue)?;
    let summary = read_uce_summary(&sample.join("uce_assembly_summary.csv"))?;
    let mut added = 0;
    for (locus, source) in reference_loci(reference)? {
        if active.is_some_and(|loci| !loci.contains(&locus)) {
            continue;
        }
        if !uce_row_accepted(summary.rows.get(&locus)) {
            continue;
        }
        let contig = sample.join("results").join(format!("{locus}.fasta"));
        if !contig.is_file() {
            continue;
        }
        let mut target = fs::OpenOptions::new()
            .append(true)
            .open(rescue.join(source.file_name().ok_or("invalid UCE reference filename")?))
            .map_err(|e| e.to_string())?;
        use std::io::Write;
        for (index, (_, sequence)) in fasta_records(&contig)?.into_iter().enumerate() {
            if sequence.len() >= minimum {
                writeln!(
                    target,
                    ">{locus}_gm2_rescue_contig_{}\n{sequence}",
                    index + 1
                )
                .map_err(|e| e.to_string())?;
                added += 1;
            }
        }
    }
    Ok(added)
}

fn build_uce_terminal_baits(
    sample: &Path,
    baits: &Path,
    active: &std::collections::BTreeSet<String>,
    window: usize,
    minimum: usize,
) -> Result<std::collections::BTreeSet<String>, String> {
    if baits.exists() {
        fs::remove_dir_all(baits).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(baits).map_err(|e| e.to_string())?;
    let summary = read_uce_summary(&sample.join("uce_assembly_summary.csv"))?;
    let mut written = std::collections::BTreeSet::new();
    for locus in active {
        if !uce_row_accepted(summary.rows.get(locus)) {
            continue;
        }
        let Some(sequence) =
            first_fasta_sequence(&sample.join("results").join(format!("{locus}.fasta")))?
        else {
            continue;
        };
        if sequence.len() < minimum {
            continue;
        }
        let flank = window.max(minimum).min(sequence.len());
        let left = &sequence[..flank];
        let right = &sequence[sequence.len() - flank..];
        let mut text = format!(">{locus}_gm2_left_terminal\n{left}\n");
        if left != right {
            text.push_str(&format!(">{locus}_gm2_right_terminal\n{right}\n"));
        }
        fs::write(baits.join(format!("{locus}.fasta")), text).map_err(|e| e.to_string())?;
        written.insert(locus.clone());
    }
    Ok(written)
}

fn restore_rescue_locus(sample: &Path, backup: &Path, locus: &str) -> Result<(), String> {
    for directory in ["results", "contigs_all", "contigs_all_low"] {
        let original = backup.join(directory).join(format!("{locus}.fasta"));
        let current = sample.join(directory).join(format!("{locus}.fasta"));
        if current.exists() {
            fs::remove_file(&current).map_err(|e| e.to_string())?;
        }
        if original.is_file() {
            if let Some(parent) = current.parent() {
                fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            fs::copy(original, current).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn row_density(row: Option<&std::collections::BTreeMap<String, String>>) -> Option<f64> {
    let length = uce_number(row, "selected_contig_length")? as f64;
    let reads =
        uce_number(row, "unique_read_count").or_else(|| uce_number(row, "read_count"))? as f64;
    (length > 0.0 && reads >= 0.0).then_some(reads / length)
}

fn write_rescue_reports(
    sample: &Sample,
    directory: &Path,
    initial: &UceSummary,
    final_rows: &UceSummary,
    rounds: &[(usize, String, UceSummary, UceSummary)],
) -> Result<(), String> {
    let mut round_csv = String::from("sample,round,locus,round_status,before_length,after_length,length_delta,before_unique_reads,after_unique_reads,unique_read_delta\n");
    for (round, status, before, after) in rounds {
        for locus in before
            .rows
            .keys()
            .chain(after.rows.keys())
            .collect::<std::collections::BTreeSet<_>>()
        {
            let left = before.rows.get(locus);
            let right = after.rows.get(locus);
            let length_delta = uce_number(right, "selected_contig_length")
                .zip(uce_number(left, "selected_contig_length"))
                .map(|(a, b)| (a - b).to_string())
                .unwrap_or_default();
            let read_delta = uce_number(right, "unique_read_count")
                .zip(uce_number(left, "unique_read_count"))
                .map(|(a, b)| (a - b).to_string())
                .unwrap_or_default();
            round_csv.push_str(&format!(
                "{},{},{},{},{},{},{},{},{},{}\n",
                sample.name,
                round,
                locus,
                status,
                left.and_then(|row| row.get("selected_contig_length"))
                    .cloned()
                    .unwrap_or_default(),
                right
                    .and_then(|row| row.get("selected_contig_length"))
                    .cloned()
                    .unwrap_or_default(),
                length_delta,
                left.and_then(|row| row.get("unique_read_count"))
                    .cloned()
                    .unwrap_or_default(),
                right
                    .and_then(|row| row.get("unique_read_count"))
                    .cloned()
                    .unwrap_or_default(),
                read_delta
            ));
        }
    }
    fs::write(directory.join("uce_rescue_rounds.csv"), round_csv).map_err(|e| e.to_string())?;
    let mut summary_csv = String::from("sample,locus,rescue_status,before_status,after_status,before_length,after_length,length_delta,before_read_density,after_read_density\n");
    for locus in initial
        .rows
        .keys()
        .chain(final_rows.rows.keys())
        .collect::<std::collections::BTreeSet<_>>()
    {
        let before = initial.rows.get(locus);
        let after = final_rows.rows.get(locus);
        let delta = uce_number(after, "selected_contig_length")
            .zip(uce_number(before, "selected_contig_length"))
            .map(|(a, b)| (a - b).to_string())
            .unwrap_or_default();
        summary_csv.push_str(&format!(
            "{},{},success,{},{},{},{},{},{:.6},{:.6}\n",
            sample.name,
            locus,
            before
                .and_then(|r| r.get("status"))
                .cloned()
                .unwrap_or_default(),
            after
                .and_then(|r| r.get("status"))
                .cloned()
                .unwrap_or_default(),
            before
                .and_then(|r| r.get("selected_contig_length"))
                .cloned()
                .unwrap_or_default(),
            after
                .and_then(|r| r.get("selected_contig_length"))
                .cloned()
                .unwrap_or_default(),
            delta,
            row_density(before).unwrap_or(0.0),
            row_density(after).unwrap_or(0.0)
        ));
    }
    fs::write(directory.join("uce_rescue_summary.csv"), summary_csv).map_err(|e| e.to_string())
}

fn execute_uce_rescue(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    sample_dir: &Path,
) -> Result<(), String> {
    let minimum = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-min-contig-length"],
        "60",
        "--uce-rescue-min-contig-length",
    )?
    .max(opt.kf.parse().unwrap_or(1));
    let maximum_rounds = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-rounds"],
        "2",
        "--uce-rescue-rounds",
    )?
    .clamp(1, 2);
    let terminal_window = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-terminal-window"],
        "350",
        "--uce-rescue-terminal-window",
    )?
    .max(minimum);
    let density_ratio = value(&opt.raw, &["--uce-rescue-min-density-ratio"], "0.5")?
        .parse::<f64>()
        .map_err(|_| "--uce-rescue-min-density-ratio must be a number")?;
    if !(0.0..=1.0).contains(&density_ratio) {
        return Err("--uce-rescue-min-density-ratio must be in [0, 1]".into());
    }
    let initial = read_uce_summary(&sample_dir.join("uce_assembly_summary.csv"))?;
    let mut current = initial.clone();
    let mut previous = initial.clone();
    let mut records: Vec<(usize, String, UceSummary)> = Vec::new();
    for round in 1..=maximum_rounds {
        let candidate = if round == 1 {
            None
        } else {
            Some(terminal_rescue_loci(&previous, &current))
        };
        if candidate.as_ref().is_some_and(|loci| loci.is_empty()) {
            break;
        }
        let root = sample_dir.join(format!("uce_rescue_round_{round}"));
        let reference = root.join("assembly_refs");
        let added = build_uce_rescue_reference(
            Path::new(&opt.reference),
            sample_dir,
            &reference,
            minimum,
            candidate.as_ref(),
        )?;
        if added == 0 {
            break;
        }
        let recruit = if let Some(active) = candidate.as_ref() {
            let terminal = root.join("terminal_baits");
            let baits =
                build_uce_terminal_baits(sample_dir, &terminal, active, terminal_window, minimum)?;
            if baits.is_empty() {
                break;
            }
            terminal
        } else {
            reference.clone()
        };
        let backup = Path::new(&opt.output)
            .join(".uce_rescue_backups")
            .join(format!("{}_r{round}", sample.name));
        if backup.exists() {
            fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
        }
        copy_tree(sample_dir, &backup)?;
        let before = current.clone();
        let result: Result<UceSummary, String> = (|| {
            let filtered = sample_dir.join("filtered");
            if filtered.exists() {
                fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
            }
            run(
                bins,
                "uce_filter",
                &uce_filter_args_for_recruit(
                    opt, sample, sample_dir, &reference, &recruit, "contig",
                ),
            )?;
            let mut rescue_opt = opt.clone();
            rescue_opt.reference = reference.display().to_string();
            run(
                bins,
                "main_assembler-rust",
                &uce_assembler_args(&rescue_opt, sample_dir)?,
            )?;
            let mut after = read_uce_summary(&sample_dir.join("uce_assembly_summary.csv"))?;
            for (locus, before_row) in &before.rows {
                let inactive = candidate
                    .as_ref()
                    .is_some_and(|active| !active.contains(locus));
                let failed =
                    uce_row_accepted(Some(before_row)) && !uce_row_accepted(after.rows.get(locus));
                let density_drop = row_density(Some(before_row))
                    .zip(row_density(after.rows.get(locus)))
                    .is_some_and(|(old, new)| old > 0.0 && new / old < density_ratio);
                if inactive || failed || density_drop {
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                }
            }
            write_uce_summary(&sample_dir.join("uce_assembly_summary.csv"), &after)?;
            Ok(after)
        })();
        match result {
            Ok(after) => {
                records.push((
                    round,
                    if round == 1 {
                        "whole-contig"
                    } else {
                        "terminal-only"
                    }
                    .into(),
                    before,
                ));
                records.push((round, "accepted".into(), after.clone()));
                previous = current;
                current = after;
            }
            Err(error) => {
                if sample_dir.exists() {
                    fs::remove_dir_all(sample_dir).map_err(|e| e.to_string())?;
                }
                copy_tree(&backup, sample_dir)?;
                eprintln!(
                    "Warning: UCE rescue round {round} rolled back for {}: {error}",
                    sample.name
                );
                break;
            }
        }
        if backup.exists() {
            fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
        }
    }
    // Pair each before/after state for compact reports.
    let pairs = records
        .chunks(2)
        .filter(|pair| pair.len() == 2)
        .map(|pair| {
            (
                pair[0].0,
                pair[0].1.clone(),
                pair[0].2.clone(),
                pair[1].2.clone(),
            )
        })
        .collect::<Vec<_>>();
    write_rescue_reports(sample, sample_dir, &initial, &current, &pairs)
}

#[derive(Clone, Default)]
struct WorkflowProfile {
    rows: Arc<Mutex<Vec<WorkflowProfileRow>>>,
}

#[derive(Clone)]
struct WorkflowProfileRow {
    sample: String,
    round: u32,
    stage: String,
    wall_ms: u128,
    input_bytes: u64,
    output_bytes: u64,
    status: &'static str,
}

fn profile_path_size(path: &Path) -> u64 {
    if path.is_file() {
        return fs::metadata(path)
            .map(|metadata| metadata.len())
            .unwrap_or(0);
    }
    if path.is_dir() {
        return directory_size(path).unwrap_or(0);
    }
    0
}

fn record_profile_named(
    profile: Option<&WorkflowProfile>,
    sample: &str,
    stage: &str,
    started: Instant,
    input_bytes: u64,
    output: &Path,
    result: &Result<(), String>,
) {
    let Some(profile) = profile else { return };
    profile
        .rows
        .lock()
        .expect("workflow profile poisoned")
        .push(WorkflowProfileRow {
            sample: sample.into(),
            round: 0,
            stage: stage.into(),
            wall_ms: started.elapsed().as_millis(),
            input_bytes,
            output_bytes: profile_path_size(output),
            status: if result.is_ok() { "ok" } else { "failed" },
        });
}

fn record_profile_stage(
    profile: Option<&WorkflowProfile>,
    sample: &Sample,
    stage: &str,
    started: Instant,
    input_bytes: u64,
    output: &Path,
    result: &Result<(), String>,
) {
    record_profile_named(
        profile,
        &sample.name,
        stage,
        started,
        input_bytes,
        output,
        result,
    );
}

fn run_profiled_action<F>(
    profile: Option<&WorkflowProfile>,
    sample: &str,
    stage: &str,
    input: &Path,
    output: &Path,
    action: F,
) -> Result<(), String>
where
    F: FnOnce() -> Result<(), String>,
{
    let input_bytes = profile_path_size(input);
    let started = Instant::now();
    let result = action();
    record_profile_named(
        profile,
        sample,
        stage,
        started,
        input_bytes,
        output,
        &result,
    );
    result
}

#[allow(clippy::too_many_arguments)]
fn run_profiled(
    profile: Option<&WorkflowProfile>,
    sample: &Sample,
    stage: &str,
    input: &Path,
    output: &Path,
    bins: &Path,
    binary: &str,
    args: &[String],
) -> Result<(), String> {
    run_profiled_action(profile, &sample.name, stage, input, output, || {
        run(bins, binary, args)
    })
}

fn execute_uce(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    profile: Option<&WorkflowProfile>,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    if opt.commands.iter().any(|c| c == "filter") {
        let args = uce_filter_args(opt, sample, &sample_dir);
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "uce_filter",
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        let args = uce_assembler_args(opt, &sample_dir)?;
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            "main_assembler-rust",
            &args,
        )?;
        if opt.rescue {
            let rescue_input_bytes = profile_path_size(&sample_dir);
            let started = Instant::now();
            let result = execute_uce_rescue(opt, bins, sample, &sample_dir);
            record_profile_stage(
                profile,
                sample,
                "rescue",
                started,
                rescue_input_bytes,
                &sample_dir,
                &result,
            );
            result?;
        }
    }
    Ok(())
}

fn execute_uce_legacy(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    dictionary: &Path,
    profile: Option<&WorkflowProfile>,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    let candidates = sample_dir.join("filtered_pe");
    if opt.commands.iter().any(|command| command == "filter") {
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-q1".into(),
            sample.read1.clone(),
        ];
        if let Some(read2) = &sample.read2 {
            args.extend(["-q2".into(), read2.clone()]);
        }
        args.extend([
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            "5".into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
        ]);
        if opt.max_reads != "0" {
            args.extend(["-m_reads".into(), opt.max_reads.clone()]);
        }
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "MainFilterNew",
            &args,
        )?;
    }
    if opt.commands.iter().any(|command| command == "refilter") {
        if !candidates.is_dir() {
            return Err("No successful filter run, cannot re-filter".into());
        }
        let args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-qd".into(),
            candidates.display().to_string(),
            "-o".into(),
            sample_dir.join("filtered").display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-p".into(),
            "1".into(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            opt.low_depth.clone(),
            "--max-depth".into(),
            opt.depth_limit.clone(),
            "--max-size".into(),
            opt.size_limit.clone(),
            "--use-gm2-format".into(),
            "--keep-linked-mates".into(),
        ];
        run_profiled(
            profile,
            sample,
            "refilter",
            &candidates,
            &sample_dir.join("filtered"),
            bins,
            "main_refilter_new",
            &args,
        )?;
    }
    if opt.commands.iter().any(|command| command == "assemble") {
        let args = uce_assembler_args(opt, &sample_dir)?;
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            "main_assembler-rust",
            &args,
        )?;
        if opt.rescue {
            let rescue_input_bytes = profile_path_size(&sample_dir);
            let started = Instant::now();
            let result = execute_uce_rescue(opt, bins, sample, &sample_dir);
            record_profile_stage(
                profile,
                sample,
                "rescue",
                started,
                rescue_input_bytes,
                &sample_dir,
                &result,
            );
            result?;
        }
    }
    Ok(())
}

fn execute_gene(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    dictionary: &Path,
    profile: Option<&WorkflowProfile>,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    let filtered_reads = sample_dir.join("filtered_pe");
    if opt.commands.iter().any(|c| c == "filter") {
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-q1".into(),
            sample.read1.clone(),
        ];
        if let Some(read2) = &sample.read2 {
            args.extend(["-q2".into(), read2.clone()]);
        }
        args.extend([
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            if sample.read2.is_some() { "5" } else { "0" }.into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
        ]);
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "MainFilterNew",
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "refilter") {
        let input_flag = if sample.read2.is_some() { "-qd" } else { "-qs" };
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            input_flag.into(),
            filtered_reads.display().to_string(),
            "-o".into(),
            sample_dir.join("filtered").display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-p".into(),
            "1".into(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            opt.low_depth.clone(),
            "--max-depth".into(),
            opt.depth_limit.clone(),
            "--max-size".into(),
            opt.size_limit.clone(),
        ];
        if sample.read2.is_some() {
            args.push("--use-gm2-format".into());
        }
        run_profiled(
            profile,
            sample,
            "refilter",
            &filtered_reads,
            &sample_dir.join("filtered"),
            bins,
            "main_refilter_new",
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
        let binary =
            match implementation.as_str() {
                "auto" | "original-rust" => "main_assembler-original-rust",
                "original" => "main_assembler-original",
                "uce-rust" => "main_assembler-rust",
                _ => return Err(
                    "--assembler-implementation must be auto, uce-rust, original, or original-rust"
                        .into(),
                ),
            };
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-o".into(),
            sample_dir.display().to_string(),
            "-ka".into(),
            opt.ka.clone(),
            "-k_min".into(),
            opt.min_ka.clone(),
            "-k_max".into(),
            opt.max_ka.clone(),
            "-limit_count".into(),
            opt.error_threshold.clone(),
            "-iteration".into(),
            opt.search_depth.clone(),
            "-sb".into(),
            soft_boundary(&opt.soft_boundary)?,
            "-cov_min".into(),
            opt.min_coverage.clone(),
            "-p".into(),
            "1".into(),
        ];
        if implementation == "uce-rust" {
            args.extend([
                "--assembly-mode".into(),
                "original".into(),
                "--assembler-read-chunk-size".into(),
                value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
            ]);
        }
        if implementation != "original" {
            if let Some(cache) = assembler_cache_directory(opt)? {
                fs::create_dir_all(&cache).map_err(|e| e.to_string())?;
                args.extend([
                    "--assembler-reference-cache-dir".into(),
                    cache.display().to_string(),
                ]);
            }
        }
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            binary,
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "gene") {
        let classify_input_bytes = profile_path_size(&sample_dir.join("contigs_all"));
        let started = Instant::now();
        let result = run(
            bins,
            "gene_workflow",
            &[
                "classify".into(),
                "--reference".into(),
                opt.reference.clone(),
                "--contigs".into(),
                sample_dir.join("contigs_all").display().to_string(),
                "--sample".into(),
                sample.name.clone(),
                "--out".into(),
                Path::new(&opt.output).join("gene").display().to_string(),
            ],
        );
        record_profile_stage(
            profile,
            sample,
            "gene-classify",
            started,
            classify_input_bytes,
            &Path::new(&opt.output).join("gene"),
            &result,
        );
        result?;
    }
    Ok(())
}

fn optional_value(args: &[String], names: &[&str]) -> Result<Option<String>, String> {
    let value = value(args, names, "")?;
    Ok((!value.is_empty()).then_some(value))
}

fn reference_loci(reference: &Path) -> Result<Vec<(String, PathBuf)>, String> {
    let entries = fs::read_dir(reference).map_err(|e| {
        format!(
            "Unable to read reference directory '{}': {e}",
            reference.display()
        )
    })?;
    let mut loci = entries
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.is_file())
        .filter_map(|path| {
            let extension = path.extension()?.to_str()?.to_ascii_lowercase();
            if matches!(extension.as_str(), "fa" | "fas" | "fasta") {
                let name = path.file_stem()?.to_str()?.to_owned();
                Some((name, path))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    loci.sort_by(|a, b| a.1.cmp(&b.1));
    Ok(loci)
}

fn reference_cache_directory(opt: &Options) -> Result<Option<PathBuf>, String> {
    let configured = optional_value(&opt.raw, &["--reference-cache-dir"])?;
    if configured.is_some() && !opt.reuse_reference_cache {
        return Err("--reference-cache-dir requires --reuse-reference-cache".into());
    }
    if !opt.reuse_reference_cache {
        return Ok(None);
    }
    Ok(Some(configured.map(PathBuf::from).unwrap_or_else(|| {
        Path::new(&opt.output).join(".gm2_reference_cache")
    })))
}

fn reference_dictionary_path(opt: &Options) -> Result<PathBuf, String> {
    let Some(cache) = reference_cache_directory(opt)? else {
        return Ok(Path::new(&opt.output).join(format!("kmer_dict_k{}.dict", opt.kf)));
    };
    let mut digest = Sha256::new();
    digest.update(opt.reference.as_bytes());
    digest.update([0]);
    digest.update(opt.kf.as_bytes());
    digest.update([0]);
    digest.update(opt.step.as_bytes());
    for (_, path) in reference_loci(Path::new(&opt.reference))? {
        digest.update([0]);
        digest.update(
            path.file_name()
                .and_then(|name| name.to_str())
                .unwrap_or_default()
                .as_bytes(),
        );
        digest.update([0]);
        let metadata = fs::metadata(&path).map_err(|e| e.to_string())?;
        digest.update(metadata.len().to_le_bytes());
        digest.update(fs::read(&path).map_err(|e| e.to_string())?);
    }
    let hex = format!("{:x}", digest.finalize());
    Ok(cache.join(format!(
        "reference_kmer_k{}_s{}_{}.dict",
        opt.kf,
        opt.step,
        &hex[..16]
    )))
}

fn assembler_cache_directory(opt: &Options) -> Result<Option<PathBuf>, String> {
    Ok(reference_cache_directory(opt)?.map(|root| root.join("assembler")))
}

fn fastx_output_extension(read: &str) -> &'static str {
    let path = Path::new(read);
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default();
    let extension = if extension.eq_ignore_ascii_case("gz") {
        path.file_stem()
            .and_then(|value| Path::new(value).extension())
            .and_then(|value| value.to_str())
            .unwrap_or_default()
    } else {
        extension
    };
    if extension.eq_ignore_ascii_case("fq") || extension.eq_ignore_ascii_case("fastq") {
        ".fq"
    } else {
        ".fasta"
    }
}

fn run_program(program: &Path, args: &[String]) -> Result<(), String> {
    let status = Command::new(program)
        .args(args)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn run_program_in(program: &Path, args: &[String], directory: &Path) -> Result<(), String> {
    let status = Command::new(program)
        .args(args)
        .current_dir(directory)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn first_fasta_sequence(path: &Path) -> Result<Option<String>, String> {
    let content = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut sequence = String::new();
    let mut active = false;
    for line in content.lines() {
        if line.starts_with('>') {
            if active {
                break;
            }
            active = true;
        } else if active {
            sequence.push_str(line.trim());
        }
    }
    Ok((!sequence.is_empty()).then_some(sequence))
}

#[derive(Clone, Debug, Default)]
struct UceSummary {
    headers: Vec<String>,
    rows: std::collections::BTreeMap<String, std::collections::BTreeMap<String, String>>,
}

fn read_uce_summary(path: &Path) -> Result<UceSummary, String> {
    if !path.is_file() {
        return Ok(UceSummary::default());
    }
    let content = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = content.lines();
    let headers = lines
        .next()
        .map(|line| line.split(',').map(str::to_owned).collect::<Vec<_>>())
        .unwrap_or_default();
    let Some(locus_index) = headers.iter().position(|field| field == "locus") else {
        return Ok(UceSummary {
            headers,
            rows: Default::default(),
        });
    };
    let mut rows = std::collections::BTreeMap::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let values = line.split(',').collect::<Vec<_>>();
        let Some(locus) = values.get(locus_index).filter(|value| !value.is_empty()) else {
            continue;
        };
        let row = headers
            .iter()
            .enumerate()
            .map(|(index, header)| {
                (
                    header.clone(),
                    values.get(index).copied().unwrap_or_default().to_owned(),
                )
            })
            .collect();
        rows.insert((*locus).to_owned(), row);
    }
    Ok(UceSummary { headers, rows })
}

fn uce_row_accepted(row: Option<&std::collections::BTreeMap<String, String>>) -> bool {
    let Some(row) = row else {
        return false;
    };
    match row
        .get("accepted")
        .map(|value| value.trim().to_ascii_lowercase())
    {
        Some(value) if !value.is_empty() => matches!(value.as_str(), "1" | "true" | "yes"),
        _ => {
            row.get("status").is_some_and(|status| status == "success")
                && !row.get("low_quality").is_some_and(|value| {
                    matches!(
                        value.trim().to_ascii_lowercase().as_str(),
                        "1" | "true" | "yes"
                    )
                })
        }
    }
}

fn uce_accepted_loci(path: &Path) -> Result<Option<std::collections::HashSet<String>>, String> {
    let summary = read_uce_summary(&path.join("uce_assembly_summary.csv"))?;
    Ok(Some(
        summary
            .rows
            .iter()
            .filter(|(_, row)| uce_row_accepted(Some(row)))
            .map(|(locus, _)| locus.clone())
            .collect(),
    ))
}

fn uce_number(row: Option<&std::collections::BTreeMap<String, String>>, key: &str) -> Option<i64> {
    row.and_then(|row| row.get(key))
        .and_then(|value| value.parse().ok())
}

fn terminal_rescue_loci(
    before: &UceSummary,
    after: &UceSummary,
) -> std::collections::BTreeSet<String> {
    after
        .rows
        .iter()
        .filter_map(|(locus, row)| {
            if !uce_row_accepted(Some(row)) {
                return None;
            }
            let previous = before.rows.get(locus);
            let length_gain = uce_number(Some(row), "selected_contig_length")
                .zip(uce_number(previous, "selected_contig_length"))
                .is_some_and(|(next, prior)| next - prior >= 30);
            let read_gain = uce_number(Some(row), "unique_read_count")
                .zip(uce_number(previous, "unique_read_count"))
                .is_some_and(|(next, prior)| next - prior >= 2);
            (previous.is_none() || !uce_row_accepted(previous) || length_gain || read_gain)
                .then(|| locus.clone())
        })
        .collect()
}

fn write_uce_summary(path: &Path, summary: &UceSummary) -> Result<(), String> {
    if summary.headers.is_empty() {
        return Ok(());
    }
    let mut text = String::new();
    text.push_str(&summary.headers.join(","));
    text.push('\n');
    for row in summary.rows.values() {
        text.push_str(
            &summary
                .headers
                .iter()
                .map(|field| row.get(field).cloned().unwrap_or_default())
                .collect::<Vec<_>>()
                .join(","),
        );
        text.push('\n');
    }
    fs::write(path, text).map_err(|e| e.to_string())
}

fn write_combined_locus(
    locus: &str,
    input_dir: &str,
    output: &Path,
    samples: &[Sample],
    uce: bool,
) -> Result<bool, String> {
    let mut records = String::new();
    for sample in samples {
        let sample_dir = output.join(&sample.name);
        if uce && uce_accepted_loci(&sample_dir)?.is_some_and(|accepted| !accepted.contains(locus))
        {
            continue;
        }
        let source = sample_dir.join(input_dir).join(format!("{locus}.fasta"));
        if !source.is_file() {
            continue;
        }
        if let Some(sequence) = first_fasta_sequence(&source)? {
            records.push('>');
            records.push_str(&sample.name);
            records.push('\n');
            records.push_str(&sequence);
            records.push('\n');
        }
    }
    if records.is_empty() {
        return Ok(false);
    }
    fs::write(
        output
            .join("combined_results")
            .join(format!("{locus}.fasta")),
        records,
    )
    .map_err(|e| e.to_string())?;
    Ok(true)
}

fn alignment_filter(raw: &[String]) -> Result<String, String> {
    if flag(raw, "--no-trimal")? {
        Ok("none".into())
    } else {
        value(raw, &["--alignment-filter"], "trimal")
    }
}

fn phylogeny_binary(program: &str) -> Result<PathBuf, String> {
    let (env_name, default) = match program {
        "raxmlng" => ("GM2_RAXMLNG", "raxml-ng"),
        "iqtree" => ("GM2_IQTREE", "iqtree"),
        "veryfasttree" => ("GM2_VERYFASTTREE", "VeryFastTree"),
        "fasttree" => ("GM2_FASTTREE", "FastTree"),
        _ => {
            return Err("--phylo-program must be raxmlng, iqtree, fasttree, or veryfasttree".into())
        }
    };
    Ok(PathBuf::from(
        env::var(env_name).unwrap_or_else(|_| default.into()),
    ))
}

fn build_tree(
    program: &str,
    binary: &Path,
    input: &Path,
    bootstrap: usize,
    threads: usize,
    quiet: bool,
) -> Result<PathBuf, String> {
    let input_text = input.display().to_string();
    let output = match program {
        "raxmlng" => format!("{input_text}.raxml.bestTree"),
        "iqtree" => format!("{input_text}.treefile"),
        "veryfasttree" => format!("{input_text}.veryfasttree.tre"),
        _ => format!("{input_text}.fasttree.tre"),
    };
    let output = PathBuf::from(output);
    if output.exists() {
        fs::remove_file(&output).map_err(|e| e.to_string())?;
    }
    let mut command = Command::new(binary);
    match program {
        "raxmlng" => {
            command.args([
                "--msa",
                &input_text,
                "--msa-format",
                "FASTA",
                "--model",
                "GTR+G",
                "--redo",
            ]);
            if bootstrap > 0 {
                command.args(["--all", "--bs-trees", &bootstrap.to_string()]);
            } else {
                command.arg("--search");
            }
            if threads > 1 {
                command.args([
                    "--threads",
                    &format!("auto{{{threads}}}"),
                    "--workers",
                    "auto",
                ]);
            } else {
                command.args(["--threads", "1"]);
            }
        }
        "iqtree" => {
            command.args(["-s", &input_text, "-redo"]);
            if bootstrap > 0 {
                command.args(["-B", &bootstrap.to_string()]);
            }
            if threads > 1 {
                command.args(["-T", "AUTO", "-ntmax", &threads.to_string()]);
            } else {
                command.args(["-T", "1"]);
            }
        }
        "veryfasttree" => {
            command.args(["-out", &output.display().to_string(), "-gtr"]);
            if bootstrap > 0 {
                command.args(["-boot", &bootstrap.to_string()]);
            } else {
                command.arg("-nosupport");
            }
            if threads > 1 {
                command.args(["-threads", &threads.to_string()]);
            }
            command.args(["-nt", &input_text]);
        }
        _ => {
            command.args(["-out", &output.display().to_string(), "-gtr"]);
            if bootstrap > 0 {
                command.args(["-boot", &bootstrap.to_string()]);
            } else {
                command.arg("-nosupport");
            }
            command.args(["-nt", &input_text]);
        }
    }
    if quiet {
        command
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null());
    }
    let status = command
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", binary.display()))?;
    if !status.success() {
        return Err(format!("{} exited with {status}", binary.display()));
    }
    if output.is_file() {
        Ok(output)
    } else {
        Err(format!(
            "{} did not create {}",
            binary.display(),
            output.display()
        ))
    }
}

fn file_sha256(path: &Path) -> Result<String, String> {
    let mut file = fs::File::open(path).map_err(|e| e.to_string())?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 65536];
    loop {
        let count = file.read(&mut buffer).map_err(|e| e.to_string())?;
        if count == 0 {
            break;
        }
        digest.update(&buffer[..count]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn raw_number<T: std::str::FromStr>(
    raw: &[String],
    names: &[&str],
    default: &str,
    label: &str,
) -> Result<T, String> {
    value(raw, names, default)?
        .parse()
        .map_err(|_| format!("{label} must be numeric"))
}

fn mito_reference(opt: &Options, bins: &Path) -> Result<PathBuf, String> {
    let raw = &opt.raw;
    let input = PathBuf::from(value(raw, &["--mito-genbank"], "")?);
    if !input.is_file() {
        return Err("mito requires a readable --mito-genbank".into());
    }
    let flank: usize = raw_number(raw, &["--mito-flank"], "150", "--mito-flank")?;
    let length: usize = raw_number(raw, &["--mito-tile-length"], "1200", "--mito-tile-length")?;
    let step: usize = raw_number(raw, &["--mito-tile-step"], "600", "--mito-tile-step")?;
    if step == 0 || length == 0 || step > length {
        return Err("mito requires 0 < --mito-tile-step <= --mito-tile-length".into());
    }
    let reference = Path::new(&opt.output).join(".gm2_mito_reference");
    if reference.exists() {
        fs::remove_dir_all(&reference).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "mito_workflow",
        &[
            "prepare-reference".into(),
            "--input".into(),
            input.display().to_string(),
            "--out-dir".into(),
            reference.display().to_string(),
            "--flank".into(),
            flank.to_string(),
            "--tile-length".into(),
            length.to_string(),
            "--tile-step".into(),
            step.to_string(),
        ],
    )?;
    Ok(reference)
}

fn mito_assembler_args(
    opt: &Options,
    reference: &Path,
    sample: &Path,
) -> Result<Vec<String>, String> {
    let raw = &opt.raw;
    let ka = value(raw, &["-ka"], "0")?;
    let ka = if ka == "0" { "31".into() } else { ka };
    let args = vec![
        "-r".into(),
        reference.display().to_string(),
        "-o".into(),
        sample.display().to_string(),
        "-ka".into(),
        ka,
        "-k_min".into(),
        value(raw, &["--min-ka"], "21")?,
        "-k_max".into(),
        value(raw, &["--max-ka"], "51")?,
        "-limit_count".into(),
        value(raw, &["-e", "--error-threshold"], "2")?,
        "-iteration".into(),
        raw_number::<usize>(raw, &["-i", "--search-depth"], "4096", "--search-depth")?
            .max(30000)
            .to_string(),
        "-sb".into(),
        "10000".into(),
        "-cov_min".into(),
        value(raw, &["--min-coverage"], "0")?,
        "-p".into(),
        "1".into(),
        "--assembly-mode".into(),
        "uce".into(),
        "--uce-side-candidates".into(),
        value(raw, &["--uce-side-candidates"], "8")?,
        "--uce-max-contig-length".into(),
        value(raw, &["--uce-max-contig-length"], "0")?,
        "--uce-min-read-density".into(),
        "0".into(),
        "--uce-density-check-min-length".into(),
        value(raw, &["--uce-density-check-min-length"], "1000")?,
        "--uce-max-depth-cv".into(),
        value(raw, &["--uce-max-depth-cv"], "0")?,
        "--uce-max-depth-ratio".into(),
        value(raw, &["--uce-max-depth-ratio"], "0")?,
        "--uce-path-strategy".into(),
        value(&opt.raw, &["--uce-path-strategy"], "backbone")?,
        "--uce-backbone-lookahead".into(),
        value(&opt.raw, &["--uce-backbone-lookahead"], "24")?,
        "--assembler-read-chunk-size".into(),
        value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
        "--assembler-kmer-count-threads".into(),
        "1".into(),
        "--assembler-graph-format".into(),
        "gfa".into(),
    ];
    Ok(args)
}

fn mito_recruit_refilter_assemble(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    sample: &Sample,
    sample_dir: &Path,
    dictionary: &Path,
    max_reads: usize,
) -> Result<(), String> {
    let raw = &opt.raw;
    let paired = sample.read2.as_ref().unwrap_or(&sample.read1);
    let candidates = sample_dir.join("filtered_pe");
    if candidates.exists() {
        fs::remove_dir_all(&candidates).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-q1".into(),
            sample.read1.clone(),
            "-q2".into(),
            paired.clone(),
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            value(raw, &["-kf"], "31")?,
            "-s".into(),
            value(raw, &["-s", "--step-size"], "4")?,
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            "4".into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m_reads".into(),
            max_reads.to_string(),
        ],
    )?;
    let collapsed = sample_dir.join("filtered_pe_collapsed");
    if collapsed.exists() {
        fs::remove_dir_all(&collapsed).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "mito_workflow",
        &[
            "collapse-baits".into(),
            "--input-dir".into(),
            candidates.display().to_string(),
            "--out-dir".into(),
            collapsed.display().to_string(),
            "--output-name".into(),
            "mitochondrion".into(),
        ],
    )?;
    fs::remove_dir_all(&candidates).map_err(|e| e.to_string())?;
    fs::rename(&collapsed, &candidates).map_err(|e| e.to_string())?;
    let filtered = sample_dir.join("filtered");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "main_refilter_new",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-qd".into(),
            candidates.display().to_string(),
            "-o".into(),
            filtered.display().to_string(),
            "-kf".into(),
            value(raw, &["-kf"], "31")?,
            "-p".into(),
            "1".into(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            value(raw, &["--depth-low-water-mark"], "50")?,
            "--max-depth".into(),
            value(raw, &["--depth-limit"], "768")?,
            "--max-size".into(),
            value(raw, &["--file-size-limit"], "6")?,
        ],
    )?;
    run(
        bins,
        "main_assembler-rust",
        &mito_assembler_args(opt, reference, sample_dir)?,
    )
}

fn copy_tree(source: &Path, destination: &Path) -> Result<(), String> {
    fs::create_dir_all(destination).map_err(|e| e.to_string())?;
    for entry in fs::read_dir(source).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let from = entry.path();
        let to = destination.join(entry.file_name());
        if from.is_dir() {
            copy_tree(&from, &to)?;
        } else {
            fs::copy(&from, &to).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn fasta_records(path: &Path) -> Result<Vec<(String, String)>, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut records = Vec::new();
    let mut id = String::new();
    let mut sequence = String::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !id.is_empty() && !sequence.is_empty() {
                records.push((std::mem::take(&mut id), std::mem::take(&mut sequence)));
            }
            id = header
                .split_whitespace()
                .next()
                .unwrap_or("sequence")
                .to_owned();
        } else if !id.is_empty() {
            sequence.push_str(line.trim());
        }
    }
    if !id.is_empty() && !sequence.is_empty() {
        records.push((id, sequence));
    }
    Ok(records)
}

fn build_mito_rescue_reference(reference: &Path, sample: &Path) -> Result<Option<PathBuf>, String> {
    let contigs = sample.join("contigs_all/mitochondrion.fasta");
    if !contigs.is_file() {
        return Ok(None);
    }
    let seeds = fasta_records(&contigs)?;
    if seeds.is_empty() {
        return Ok(None);
    }
    let rescue = sample.join("mito_rescue_round_1/assembly_refs");
    if rescue.exists() {
        fs::remove_dir_all(&rescue).map_err(|e| e.to_string())?;
    }
    copy_tree(reference, &rescue)?;
    let mut bait = fs::OpenOptions::new()
        .append(true)
        .open(rescue.join("mitochondrion.fasta"))
        .map_err(|e| e.to_string())?;
    use std::io::Write;
    for (index, (_, sequence)) in seeds.iter().enumerate() {
        if sequence.len() >= 31 {
            writeln!(bait, ">sample_seed_{index}\n{sequence}").map_err(|e| e.to_string())?;
        }
    }
    Ok(Some(rescue))
}

fn build_mito_dictionary(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    dictionary: &Path,
    output: &Path,
) -> Result<(), String> {
    if let Some(parent) = dictionary.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-o".into(),
            output.display().to_string(),
            "-kf".into(),
            value(&opt.raw, &["-kf"], "31")?,
            "-s".into(),
            value(&opt.raw, &["-s", "--step-size"], "4")?,
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ],
    )
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            _ => 'N',
        })
        .collect()
}
fn minimal_rotation(sequence: &str) -> String {
    let bytes = sequence.as_bytes();
    let n = bytes.len();
    if n == 0 {
        return String::new();
    }
    let mut doubled = Vec::with_capacity(n * 2);
    doubled.extend_from_slice(bytes);
    doubled.extend_from_slice(bytes);
    let (mut left, mut right, mut offset) = (0usize, 1usize, 0usize);
    while left < n && right < n && offset < n {
        let (a, b) = (doubled[left + offset], doubled[right + offset]);
        if a == b {
            offset += 1;
            continue;
        }
        if a > b {
            left += offset + 1;
            if left == right {
                left += 1;
            }
        } else {
            right += offset + 1;
            if left == right {
                right += 1;
            }
        }
        offset = 0;
    }
    String::from_utf8(doubled[left.min(right)..left.min(right) + n].to_vec()).expect("DNA UTF-8")
}
fn canonical_circular(sequence: &str) -> String {
    let sequence = sequence.to_ascii_uppercase();
    minimal_rotation(&sequence).min(minimal_rotation(&reverse_complement(&sequence)))
}
fn mito_observation(sample: &Path) -> Result<(String, String), String> {
    let summary = fs::read_to_string(sample.join("mito/mitochondrial_assembly_summary.tsv"))
        .unwrap_or_default();
    let status = summary
        .lines()
        .find_map(|line| line.strip_prefix("status\t"))
        .unwrap_or("missing")
        .to_owned();
    let fasta = sample.join("mito/mitochondrial_assembly.fasta");
    let evidence = if status == "circular" {
        fasta_records(&fasta)?
            .first()
            .map(|(_, s)| canonical_circular(s))
            .unwrap_or_default()
    } else if fasta.is_file() {
        file_sha256(&fasta)?
    } else {
        String::new()
    };
    Ok((status, evidence))
}

fn finalize_mito_sample(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    sample_dir: &Path,
    require_circular: bool,
) -> Result<(), String> {
    let raw = &opt.raw;
    let mut args = vec![
        "finalize".into(),
        "--reference-genome".into(),
        reference
            .join("metadata/mitochondrial_reference.fasta")
            .display()
            .to_string(),
        "--gene-metadata".into(),
        reference
            .join("metadata/mitochondrial_genes.tsv")
            .display()
            .to_string(),
        "--contigs".into(),
        sample_dir
            .join("contigs_all/mitochondrion.fasta")
            .display()
            .to_string(),
        "--paired-reads".into(),
        sample_dir
            .join("filtered/mitochondrion.fq")
            .display()
            .to_string(),
        "--out-dir".into(),
        sample_dir.join("mito").display().to_string(),
        "--minimum-overlap".into(),
        value(raw, &["--mito-min-overlap"], "41")?,
        "--minimum-identity".into(),
        value(raw, &["--mito-min-overlap-identity"], "0.98")?,
        "--terminal-window".into(),
        value(raw, &["--mito-terminal-window"], "500")?,
        "--link-kmer".into(),
        value(raw, &["--mito-link-kmer"], "31")?,
        "--minimum-link-hits".into(),
        value(raw, &["--mito-min-link-hits"], "2")?,
        "--minimum-pair-support".into(),
        value(raw, &["--mito-min-pair-support"], "3")?,
        "--bridge-kmer".into(),
        value(raw, &["--mito-bridge-kmer"], "31")?,
        "--bridge-minimum-depth".into(),
        value(raw, &["--mito-bridge-min-depth"], "2")?,
        "--maximum-bridge".into(),
        value(raw, &["--mito-max-bridge"], "1000")?,
        "--minimum-junction-support".into(),
        value(raw, &["--mito-min-junction-support"], "3")?,
        "--require-circular".into(),
        require_circular.to_string(),
    ];
    let graph = sample_dir.join("assembly_graphs/mitochondrion.gfa");
    if graph.is_file() {
        args.extend(["--graph".into(), graph.display().to_string()]);
    }
    run(bins, "mito_workflow", &args)
}

fn execute_mito_single_stage(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    reference: &Path,
    output: &Path,
    max_reads: usize,
    require_circular: bool,
) -> Result<(), String> {
    let dictionary = output.join(format!(
        "mito_kmer_dict_k{}.dict",
        value(&opt.raw, &["-kf"], "31")?
    ));
    fs::create_dir_all(output).map_err(|e| e.to_string())?;
    build_mito_dictionary(opt, bins, reference, &dictionary, output)?;
    let failures = Arc::new(Mutex::new(Vec::new()));
    let queued_samples = samples.to_vec();
    let next = Arc::new(Mutex::new(queued_samples.into_iter()));
    let stage_options = opt.clone();
    let mut handles = Vec::new();
    for _ in 0..opt.workers.min(samples.len()).max(1) {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let bins = bins.to_path_buf();
        let reference = reference.to_path_buf();
        let dictionary = dictionary.clone();
        let staged = output.to_path_buf();
        let mut stage_opt = stage_options.clone();
        stage_opt.output = staged.display().to_string();
        stage_opt.workers = 1;
        handles.push(thread::spawn(move || {
            let Some(sample) = next.lock().expect("mito queue poisoned").next() else {
                return;
            };
            let sample_dir = staged.join(&sample.name);
            let result = mito_recruit_refilter_assemble(
                &stage_opt,
                &bins,
                &reference,
                &sample,
                &sample_dir,
                &dictionary,
                max_reads,
            )
            .and_then(|_| {
                let backup = staged.join(".mito_seed_backups").join(&sample.name);
                if backup.exists() {
                    fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
                }
                copy_tree(&sample_dir, &backup)?;
                let rescue_result = build_mito_rescue_reference(&reference, &sample_dir)?.map_or(
                    Ok(()),
                    |rescue| {
                        let rescue_root =
                            rescue.parent().ok_or("rescue reference has no parent")?;
                        let rescue_dict = rescue_root.join("filter.dict");
                        build_mito_dictionary(
                            &stage_opt,
                            &bins,
                            &rescue,
                            &rescue_dict,
                            rescue_root,
                        )?;
                        mito_recruit_refilter_assemble(
                            &stage_opt,
                            &bins,
                            &rescue,
                            &sample,
                            &sample_dir,
                            &rescue_dict,
                            max_reads,
                        )
                    },
                );
                if rescue_result.is_err() {
                    if sample_dir.exists() {
                        fs::remove_dir_all(&sample_dir).map_err(|e| e.to_string())?;
                    }
                    copy_tree(&backup, &sample_dir)?;
                }
                if backup.exists() {
                    fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
                }
                finalize_mito_sample(&stage_opt, &bins, &reference, &sample_dir, require_circular)
            });
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("mito failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust mito worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "mito failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} mitochondrial sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn inferred_ipyrad_loci(params: &Path) -> Result<PathBuf, String> {
    let text = fs::read_to_string(params)
        .map_err(|e| format!("Unable to read ipyrad params '{}': {e}", params.display()))?;
    let values = text
        .lines()
        .filter_map(|line| line.split("##").next())
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if values.len() < 2 {
        return Err("ipyrad params must contain assembly_name [0] and project_dir [1]".into());
    }
    let project = PathBuf::from(&values[1]);
    let project = if project.is_absolute() {
        project
    } else {
        params
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(project)
    };
    Ok(project
        .join(format!("{}_outfiles", values[0]))
        .join(format!("{}.loci", values[0])))
}

fn rad_loci_input(opt: &Options) -> Result<(PathBuf, String), String> {
    let supplied = value(&opt.raw, &["--ipyrad-loci"], "")?;
    let params = value(&opt.raw, &["--ipyrad-params"], "")?;
    if !params.is_empty() {
        let params_path = PathBuf::from(&params);
        if !params_path.is_file() {
            return Err("--ipyrad-params must name a readable params file".into());
        }
        let executable = value(&opt.raw, &["--ipyrad-executable"], "ipyrad")?;
        let steps = value(&opt.raw, &["--ipyrad-steps"], "1234567")?;
        if steps.is_empty()
            || !steps
                .bytes()
                .all(|byte| byte.is_ascii_digit() && byte != b'0')
        {
            return Err("--ipyrad-steps must be a non-empty sequence of steps 1-7".into());
        }
        if !steps.bytes().all(|byte| matches!(byte, b'1'..=b'7')) {
            return Err("--ipyrad-steps may contain only digits 1 through 7".into());
        }
        let status = Command::new(&executable)
            .arg("-p")
            .arg(&params)
            .arg("-s")
            .arg(&steps)
            .status()
            .map_err(|e| format!("Unable to start ipyrad executable '{executable}': {e}"))?;
        if !status.success() {
            return Err(format!("ipyrad assembly failed with status {status}"));
        }
        let loci = if supplied.is_empty() {
            inferred_ipyrad_loci(&params_path)?
        } else {
            PathBuf::from(supplied)
        };
        if !loci.is_file() {
            return Err(format!(
                "ipyrad completed but no .loci file was found at '{}'; pass --ipyrad-loci explicitly if its output was relocated",
                loci.display()
            ));
        }
        return Ok((loci, format!("ipyrad params={params} steps={steps}")));
    }
    let loci = PathBuf::from(supplied);
    if !loci.is_file() {
        return Err("provide a readable --ipyrad-loci FILE, or --ipyrad-params FILE to assemble raw RAD reads with ipyrad".into());
    }
    Ok((loci, "existing ipyrad .loci".into()))
}

fn build_rad_reference(opt: &Options, bins: &Path, reference: &Path) -> Result<PathBuf, String> {
    let (loci, source) = rad_loci_input(opt)?;
    run(
        bins,
        "rad_workflow",
        &[
            "reference".into(),
            "--loci".into(),
            loci.display().to_string(),
            "--out".into(),
            reference.display().to_string(),
        ],
    )?;
    fs::write(
        reference.join("PROVENANCE.txt"),
        format!("source\t{source}\nloci\t{}\n", loci.display()),
    )
    .map_err(|e| e.to_string())?;
    Ok(reference.to_path_buf())
}

fn execute_rad_probe(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad-probe"] {
        return Err("rad-probe cannot be combined with other subcommands".into());
    }
    let root = Path::new(&opt.output);
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    let reference = root.join("rad_reference");
    if flag(&opt.raw, "--rad-denovo")? {
        if !value(&opt.raw, &["--ipyrad-loci", "--ipyrad-params"], "")?.is_empty() {
            return Err(
                "--rad-denovo cannot be combined with --ipyrad-loci or --ipyrad-params".into(),
            );
        }
        if opt.samples.is_empty() {
            return Err("rad-probe --rad-denovo requires -f paired_rad_samples.tsv".into());
        }
        let samples = read_rad_samples(&opt.samples)?;
        let mut args = vec![
            "denovo".into(),
            "--out".into(),
            reference.display().to_string(),
        ];
        for sample in &samples {
            args.extend([
                "--sample".into(),
                sample.name.clone(),
                "--read1".into(),
                sample.read1.clone(),
                "--read2".into(),
                sample.read2.clone().expect("paired RAD samples validated"),
            ]);
        }
        let options = [
            ("--rad-overhang", "--overhang"),
            ("--rad-overhang-r2", "--overhang-r2"),
            ("--rad-kmer", "--kmer"),
            ("--rad-min-count", "--min-count"),
            ("--rad-min-samples", "--min-samples"),
            ("--rad-min-length", "--min-length"),
        ];
        for (source, target) in options {
            if let Some(value) = optional_value(&opt.raw, &[source])? {
                args.extend([target.into(), value]);
            }
        }
        run(bins, "rad_workflow", &args)?;
        fs::write(
            reference.join("PROVENANCE.txt"),
            "source\tdenovo_candidate_probe\nmode\tcanonical_solid_kmer_paired_arms\n",
        )
        .map_err(|e| e.to_string())?;
    } else {
        if !opt.samples.is_empty() {
            return Err("rad-probe uses no -f unless --rad-denovo is selected".into());
        }
        build_rad_reference(opt, bins, &reference)?;
    }
    Ok(())
}

fn execute_rad_validate(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad-validate"] {
        return Err("rad-validate cannot be combined with other subcommands".into());
    }
    if !opt.samples.is_empty() {
        return Err("rad-validate discovers samples from --rad-recovery; do not pass -f".into());
    }
    let reference = value(&opt.raw, &["--rad-probe"], "")?;
    if reference.is_empty() || !Path::new(&reference).join("arms").is_dir() {
        return Err("rad-validate requires --rad-probe DIR containing arms/".into());
    }
    let recovery = value(&opt.raw, &["--rad-recovery"], "")?;
    if recovery.is_empty() || !Path::new(&recovery).is_dir() {
        return Err("rad-validate requires --rad-recovery DIR from a completed rad run".into());
    }
    let mut args = vec![
        "validate".into(),
        "--reference".into(),
        reference,
        "--recovery".into(),
        recovery,
        "--out".into(),
        Path::new(&opt.output)
            .join("rad_validated")
            .display()
            .to_string(),
    ];
    for (source, target) in [
        ("--rad-validate-min-identity", "--min-identity"),
        ("--rad-validate-min-breadth", "--min-breadth"),
        ("--rad-validate-min-delta", "--min-delta"),
    ] {
        if let Some(value) = optional_value(&opt.raw, &[source])? {
            args.extend([target.into(), value]);
        }
    }
    run(bins, "rad_workflow", &args)
}

fn execute_rad(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad"] {
        return Err(
            "rad is a complete workflow and cannot be combined with other subcommands".into(),
        );
    }
    let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
    if !matches!(
        implementation.as_str(),
        "auto" | "original" | "original-rust"
    ) {
        return Err(
            "rad requires --assembler-implementation auto, original, or original-rust".into(),
        );
    }
    let min_arm_breadth = value(&opt.raw, &["--rad-min-arm-breadth"], "0.80")?;
    let breadth = min_arm_breadth
        .parse::<f64>()
        .map_err(|_| "--rad-min-arm-breadth must be a number")?;
    if !(0.0..=1.0).contains(&breadth) {
        return Err("--rad-min-arm-breadth must be in [0, 1]".into());
    }
    let root = Path::new(&opt.output);
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    let provided_reference = value(&opt.raw, &["--rad-probe"], "")?;
    if !provided_reference.is_empty()
        && (!value(&opt.raw, &["--ipyrad-loci"], "")?.is_empty()
            || !value(&opt.raw, &["--ipyrad-params"], "")?.is_empty())
    {
        return Err("rad accepts either --rad-probe or an ipyrad input, not both".into());
    }
    let reference = if provided_reference.is_empty() {
        build_rad_reference(opt, bins, &root.join("rad_reference"))?
    } else {
        let path = PathBuf::from(provided_reference);
        if !path.join("arms").is_dir() {
            return Err("--rad-probe must name a rad_reference directory containing arms/".into());
        }
        path
    };
    let recovery = root.join("rad_recovery");
    fs::create_dir_all(&recovery).map_err(|e| e.to_string())?;
    let samples = read_samples(&opt.samples, &recovery)?;
    let mut stage = opt.clone();
    stage.reference = reference.join("arms").display().to_string();
    stage.output = recovery.display().to_string();
    stage.assembly_mode = "original".into();
    stage.commands = vec!["filter".into(), "refilter".into(), "assemble".into()];
    let dictionary = recovery.join(format!("rad_kmer_dict_k{}.dict", stage.kf));
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            stage.reference.clone(),
            "-o".into(),
            stage.output.clone(),
            "-kf".into(),
            stage.kf.clone(),
            "-s".into(),
            stage.step.clone(),
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ],
    )?;
    let failures = Arc::new(Mutex::new(Vec::new()));
    let pending = Arc::new(Mutex::new(samples.clone().into_iter()));
    let mut handles = Vec::new();
    for _ in 0..stage.workers.min(samples.len()).max(1) {
        let pending = Arc::clone(&pending);
        let failures = Arc::clone(&failures);
        let stage = stage.clone();
        let bins = bins.to_path_buf();
        let dictionary = dictionary.clone();
        handles.push(thread::spawn(move || loop {
            let Some(sample) = pending.lock().expect("rad sample queue poisoned").next() else {
                break;
            };
            if let Err(error) = execute_gene(&stage, &bins, &sample, &dictionary, None) {
                failures
                    .lock()
                    .expect("rad failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust rad worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "rad failures poisoned")?;
    if !failures.is_empty() {
        return Err(format!(
            "{} RAD sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ));
    }
    let mut finalize = vec![
        "finalize".into(),
        "--reference".into(),
        reference.display().to_string(),
        "--recovery".into(),
        recovery.display().to_string(),
        "--out".into(),
        root.join("rad_matrix").display().to_string(),
        "--min-arm-breadth".into(),
        min_arm_breadth,
    ];
    for sample in &samples {
        finalize.extend(["--sample".into(), sample.name.clone()]);
    }
    run(bins, "rad_workflow", &finalize)?;
    if stage.cleanup_intermediates {
        cleanup_native_intermediates(&stage, &samples)?;
    }
    Ok(())
}

fn execute_mito(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    if opt.commands != ["mito"] {
        return Err(
            "mito is a complete workflow and cannot be combined with other subcommands".into(),
        );
    }
    let reference = mito_reference(opt, bins)?;
    let initial: usize = raw_number(
        &opt.raw,
        &["--mito-initial-reads"],
        "10",
        "--mito-initial-reads",
    )?;
    let maximum: usize = raw_number(&opt.raw, &["--mito-max-reads"], "320", "--mito-max-reads")?;
    if initial == 0 || maximum < initial {
        return Err("--mito-max-reads must be at least --mito-initial-reads".into());
    }
    if flag(&opt.raw, "--no-mito-adaptive-stop")? {
        return execute_mito_single_stage(
            opt,
            bins,
            samples,
            &reference,
            Path::new(&opt.output),
            initial,
            true,
        );
    }
    let root = Path::new(&opt.output);
    let stages = root.join(".mito_adaptive");
    let mut previous: Option<std::collections::BTreeMap<String, (String, String)>> = None;
    let mut limit = initial;
    loop {
        let stage = stages.join(format!("{limit}m"));
        if stage.exists() {
            fs::remove_dir_all(&stage).map_err(|e| e.to_string())?;
        }
        execute_mito_single_stage(opt, bins, samples, &reference, &stage, limit, false)?;
        let mut current = std::collections::BTreeMap::new();
        for sample in samples {
            current.insert(
                sample.name.clone(),
                mito_observation(&stage.join(&sample.name))?,
            );
        }
        let stable = previous
            .as_ref()
            .is_some_and(|previous| previous == &current);
        if stable || limit >= maximum {
            for sample in samples {
                let destination = root.join(&sample.name);
                if destination.exists() {
                    fs::remove_dir_all(&destination).map_err(|e| e.to_string())?;
                }
                copy_tree(&stage.join(&sample.name), &destination)?;
            }
            if stable && current.values().all(|(status, _)| status == "circular") {
                return Ok(());
            }
            let statuses = current
                .iter()
                .map(|(sample, (status, _))| format!("{sample}={status}"))
                .collect::<Vec<_>>()
                .join(", ");
            return Err(if stable {
                format!("mito adaptive stop reached a stable non-circular assembly; preserved the final partial result ({statuses})")
            } else {
                format!("mito adaptive stop did not confirm a stable circular assembly by {limit}M reads; {statuses}")
            });
        }
        previous = Some(current);
        limit = (limit.saturating_mul(2)).min(maximum);
    }
}

fn profile_cache_key(paths: &[&str], kmer: &str) -> Result<String, String> {
    let mut digest = Sha256::new();
    digest.update(kmer.as_bytes());
    for path in paths.iter().filter(|path| !path.is_empty()) {
        let resolved = fs::canonicalize(path).unwrap_or_else(|_| PathBuf::from(path));
        digest.update(b"\0");
        digest.update(resolved.as_os_str().as_encoded_bytes());
        if resolved.is_file() {
            let mut file = fs::File::open(&resolved).map_err(|e| e.to_string())?;
            let mut buffer = [0u8; 65536];
            loop {
                let n = file.read(&mut buffer).map_err(|e| e.to_string())?;
                if n == 0 {
                    break;
                }
                digest.update(&buffer[..n]);
            }
        }
    }
    Ok(format!("{:x}", digest.finalize())[..16].to_owned())
}

fn materialize_profile_reference(opt: &Options) -> Result<(PathBuf, PathBuf), String> {
    let input = PathBuf::from(&opt.reference);
    if !input.is_file() {
        return Err("profiling requires -r to be exactly one marker .fa/.fasta file".into());
    }
    let extension = input
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if !matches!(extension.as_str(), "fa" | "fasta") {
        return Err("profiling reference must use the .fa or .fasta extension".into());
    }
    let directory = Path::new(&opt.output).join(".marker_profile_reference");
    fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let target = directory.join(
        input
            .file_name()
            .ok_or("invalid profiling reference filename")?,
    );
    if target.exists() {
        fs::remove_file(&target).map_err(|e| e.to_string())?;
    }
    if fs::hard_link(&input, &target).is_err() {
        fs::copy(&input, &target).map_err(|e| e.to_string())?;
    }
    Ok((target, directory))
}

fn execute_profiling(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    if opt.commands != ["profiling"] {
        return Err(
            "profiling is a complete marker workflow and cannot be combined with other subcommands"
                .into(),
        );
    }
    let raw = &opt.raw;
    let kmer = value(raw, &["--profile-kmer-size"], "21")?;
    let kmer_number = kmer
        .parse::<usize>()
        .map_err(|_| "--profile-kmer-size must be an odd integer from 15 to 31")?;
    let threshold = value(raw, &["--profile-pseudoalign-threshold"], "0.8")?
        .parse::<f64>()
        .map_err(|_| "--profile-pseudoalign-threshold must be a number")?;
    let relevant = value(raw, &["--profile-relevant-kmer-fraction"], "0.5")?
        .parse::<f64>()
        .map_err(|_| "--profile-relevant-kmer-fraction must be a number")?;
    let memory = value(raw, &["--profile-index-memory-gb"], "2")?
        .parse::<usize>()
        .map_err(|_| "--profile-index-memory-gb must be an integer")?;
    let valid_parameters = (15..=31).contains(&kmer_number)
        && !kmer_number.is_multiple_of(2)
        && 0.0 < threshold
        && threshold <= 1.0
        && (0.0..=1.0).contains(&relevant)
        && memory > 0;
    if !valid_parameters {
        return Err("invalid profiling parameters".into());
    }
    let group = optional_value(raw, &["--profile-group-map"])?;
    if group
        .as_ref()
        .is_some_and(|path| !Path::new(path).is_file())
    {
        return Err("--profile-group-map must be a readable TSV file".into());
    }
    let decoy = optional_value(raw, &["--profile-decoy"])?;
    let (reference, reference_dir) = materialize_profile_reference(opt)?;
    let themisto = optional_value(raw, &["--profile-themisto"])?
        .or_else(|| env::var("GM2_THEMISTO").ok())
        .unwrap_or_else(|| "themisto".into());
    let key = profile_cache_key(
        &[
            &reference.display().to_string(),
            group.as_deref().unwrap_or(""),
            decoy.as_deref().unwrap_or(""),
            &themisto,
        ],
        &kmer,
    )?;
    let cache_root = optional_value(raw, &["--profile-index-dir"])?
        .or_else(|| {
            optional_value(raw, &["--reference-cache-dir"])
                .ok()
                .flatten()
        })
        .unwrap_or_else(|| {
            Path::new(&opt.output)
                .join(".gm2_reference_cache")
                .display()
                .to_string()
        });
    let cache = Path::new(&cache_root).join(format!("profile_themisto_k{kmer}_{key}"));
    let failures = Arc::new(Mutex::new(Vec::new()));
    let queued_samples = samples.to_vec();
    let next = Arc::new(Mutex::new(queued_samples.into_iter()));
    let mut handles = Vec::new();
    for _ in 0..opt.workers.min(samples.len()).max(1) {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let reference_dir = reference_dir.clone();
        let reference = reference.clone();
        let cache = cache.clone();
        let themisto = themisto.clone();
        let group = group.clone();
        let decoy = decoy.clone();
        let kmer = kmer.clone();
        let output_root = opt.output.clone();
        let step = opt.step.clone();
        let low_depth = opt.low_depth.clone();
        let depth_limit = opt.depth_limit.clone();
        let size_limit = opt.size_limit.clone();
        let max_reads = opt.max_reads.clone();
        let bins = bins.to_path_buf();
        let force = flag(raw, "--profile-force-rebuild")?;
        handles.push(thread::spawn(move || {
            let Some(sample) = next.lock().expect("profiling queue poisoned").next() else {
                return;
            };
            let sample_dir = Path::new(&output_root).join(&sample.name);
            let mut filter = vec![
                "-r".into(),
                reference_dir.display().to_string(),
                "--recruit-references".into(),
                reference_dir.display().to_string(),
                "-q1".into(),
                sample.read1.clone(),
            ];
            if let Some(read2) = &sample.read2 {
                filter.extend(["-q2".into(), read2.clone()]);
            }
            filter.extend([
                "-o".into(),
                sample_dir.display().to_string(),
                "-kf".into(),
                kmer.clone(),
                "-s".into(),
                step,
                "--selection".into(),
                "auto".into(),
                "--reference-role".into(),
                "bait".into(),
                "--threads".into(),
                "1".into(),
                "--memory-limit-mib".into(),
                "256".into(),
                "--min-depth".into(),
                low_depth,
                "--max-depth".into(),
                depth_limit,
                "--max-size".into(),
                size_limit,
            ]);
            if max_reads != "0" {
                filter.extend(["--max-fragments".into(), max_reads]);
            }
            let result = run(&bins, "uce_filter", &filter).and_then(|_| {
                let filtered = sample_dir.join("filtered");
                let reads = fs::read_dir(&filtered)
                    .map_err(|e| e.to_string())?
                    .filter_map(Result::ok)
                    .map(|entry| entry.path())
                    .filter(|path| path.is_file())
                    .filter(|path| {
                        matches!(
                            path.extension()
                                .and_then(|x| x.to_str())
                                .unwrap_or_default()
                                .to_ascii_lowercase()
                                .as_str(),
                            "fq" | "fastq" | "fasta" | "fa"
                        )
                    })
                    .collect::<Vec<_>>();
                if reads.len() != 1 {
                    return Err("profiling requires exactly one merged recruited-read file".into());
                }
                let profile = sample_dir.join("marker_profile");
                if profile.exists() {
                    fs::remove_dir_all(&profile).map_err(|e| e.to_string())?;
                }
                let mut args = vec![
                    "--reference".into(),
                    reference.display().to_string(),
                    "--reads".into(),
                    reads[0].display().to_string(),
                    "--output".into(),
                    profile.display().to_string(),
                    "--cache".into(),
                    cache.display().to_string(),
                    "--themisto".into(),
                    themisto,
                    "--threads".into(),
                    "1".into(),
                    "--kmer-size".into(),
                    kmer,
                    "--threshold".into(),
                    threshold.to_string(),
                    "--relevant-kmer-fraction".into(),
                    relevant.to_string(),
                    "--index-memory-gb".into(),
                    memory.to_string(),
                ];
                if let Some(group) = group {
                    args.extend(["--groups".into(), group]);
                }
                if let Some(decoy) = decoy {
                    args.extend(["--decoy".into(), decoy]);
                }
                if force {
                    args.push("--force-rebuild".into());
                }
                run(&bins, "marker_profile", &args)?;
                if profile.join("marker_reference_support.tsv").is_file() {
                    Ok(())
                } else {
                    Err("profiling failed to produce marker_reference_support.tsv".into())
                }
            });
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("profiling failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle
            .join()
            .map_err(|_| "Rust profiling worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "profiling failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_gene_annotate(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let input = value(raw, &["--gene-input"], "")?;
    let proteins = value(raw, &["--gene-protein-reference"], "")?;
    if !Path::new(&input).is_dir() {
        return Err("--gene-input must be a gene output directory".into());
    }
    if !Path::new(&proteins).is_dir() {
        return Err("gene-annotate requires --gene-protein-reference".into());
    }
    run(
        bins,
        "gene_workflow",
        &[
            "annotate".into(),
            "--input".into(),
            input,
            "--protein-reference".into(),
            proteins,
            "--out".into(),
            opt.output.clone(),
            "--miniprot".into(),
            value(raw, &["--gene-miniprot"], "miniprot")?,
            "--threads".into(),
            opt.workers.to_string(),
        ],
    )
}

fn execute_gene_resolve(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let input = value(raw, &["--gene-input"], "")?;
    if !Path::new(&input).is_dir() {
        return Err("--gene-input must be an annotation directory".into());
    }
    let mut args = vec![
        "resolve".into(),
        "--input".into(),
        input,
        "--out".into(),
        opt.output.clone(),
        "--mafft".into(),
        value(raw, &["--gene-mafft"], "mafft")?,
        "--iqtree".into(),
        value(raw, &["--gene-iqtree"], "iqtree")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--min-taxa".into(),
        value(raw, &["--gene-min-taxa"], "4")?,
        "--min-aa-length".into(),
        value(raw, &["--gene-min-aa-length"], "30")?,
        "--min-effective-codon-sites".into(),
        value(raw, &["--gene-min-effective-codon-sites"], "30")?,
    ];
    if let Some(path) = optional_value(raw, &["--gene-outgroup"])? {
        if !Path::new(&path).is_file() {
            return Err("--gene-outgroup must be a readable file".into());
        }
        args.extend(["--outgroup".into(), path]);
    }
    let ufboot = value(raw, &["--gene-ufboot"], "0")?;
    if ufboot != "0" {
        args.extend(["--ufboot".into(), ufboot]);
    }
    if let Some(path) = optional_value(raw, &["--gene-taper"])? {
        if !Path::new(&path).is_file() {
            return Err("--gene-taper must be a readable correction_multi.jl script".into());
        }
        args.extend([
            "--taper-script".into(),
            path,
            "--julia".into(),
            value(raw, &["--gene-julia"], "julia")?,
        ]);
    }
    run(bins, "gene_workflow", &args)
}

fn execute_gene_tree(opt: &Options) -> Result<(), String> {
    let raw = &opt.raw;
    let input = PathBuf::from(value(raw, &["--gene-input"], "")?);
    if !input.is_dir() {
        return Err("--gene-input must be a gene-resolve output directory".into());
    }
    let mode = value(raw, &["--gene-species-mode"], "strict")?;
    if !matches!(mode.as_str(), "strict" | "multicopy") {
        return Err("--gene-species-mode must be strict or multicopy".into());
    }
    let (trees, mapping, output_name) = if mode == "strict" {
        (
            input.join("astral_input/resolved_1to1.trees"),
            None,
            "gene_strict_aster.tree",
        )
    } else {
        (
            input.join("astralpro_input/multicopy.trees"),
            Some(input.join("astralpro_input/leaf_to_species.tsv")),
            "gene_multicopy_aster.tree",
        )
    };
    if !trees.is_file() || fs::metadata(&trees).map_err(|e| e.to_string())?.len() == 0 {
        return Err(format!(
            "No usable {mode} gene trees found: {}",
            trees.display()
        ));
    }
    if let Some(mapping) = &mapping {
        if !mapping.is_file() {
            return Err(format!("Missing multicopy leaf map: {}", mapping.display()));
        }
    }
    fs::create_dir_all(&opt.output).map_err(|e| e.to_string())?;
    let aster = PathBuf::from(value(raw, &["--gene-aster"], "astral")?);
    let output = Path::new(&opt.output).join(output_name);
    let log = Path::new(&opt.output).join(format!("{output_name}.log"));
    if output.exists() {
        fs::remove_file(&output).map_err(|e| e.to_string())?;
    }
    let mut args = vec![
        "-i".into(),
        trees.display().to_string(),
        "-o".into(),
        output.display().to_string(),
        "-t".into(),
        opt.workers.to_string(),
    ];
    if let Some(mapping) = &mapping {
        args.extend(["-a".into(), mapping.display().to_string()]);
    }
    let file = fs::File::create(&log).map_err(|e| e.to_string())?;
    let status = Command::new(&aster)
        .args(&args)
        .stdout(file.try_clone().map_err(|e| e.to_string())?)
        .stderr(file)
        .status()
        .map_err(|e| format!("Cannot find ASTER2 executable: {}: {e}", aster.display()))?;
    if !status.success() {
        return Err(format!(
            "ASTER2 exited with {status}; inspect {}",
            log.display()
        ));
    }
    let tree = fs::read_to_string(&output)
        .map_err(|_| {
            format!(
                "ASTER2 completed without a species tree; inspect {}",
                log.display()
            )
        })?
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .unwrap_or("")
        .to_owned();
    if !tree.starts_with('(') || !tree.ends_with(';') {
        return Err(format!(
            "ASTER2 output is not a Newick tree; inspect {}",
            log.display()
        ));
    }
    let mut provenance = format!("field\tvalue\nmode\t{mode}\naster_executable\t{}\ncommand\t{}\ngene_trees\t{}\ngene_trees_sha256\t{}\nspecies_tree\t{}\nspecies_tree_sha256\t{}\n", aster.display(), args.join(" "), trees.display(), file_sha256(&trees)?, output.display(), file_sha256(&output)?);
    if let Some(mapping) = mapping {
        provenance.push_str(&format!(
            "leaf_to_species\t{}\nleaf_to_species_sha256\t{}\n",
            mapping.display(),
            file_sha256(&mapping)?
        ));
    }
    fs::write(
        Path::new(&opt.output).join("gene_tree_provenance.tsv"),
        provenance,
    )
    .map_err(|e| e.to_string())
}

fn execute_tree(opt: &Options) -> Result<(), String> {
    let raw = &opt.raw;
    let method = value(raw, &["-m", "--tree-method"], "coalescent")?;
    if !matches!(method.as_str(), "coalescent" | "concatenation") {
        return Err("--tree-method must be coalescent or concatenation".into());
    }
    let program = value(raw, &["--phylo-program"], "fasttree")?;
    let binary = phylogeny_binary(&program)?;
    let filter = alignment_filter(raw)?;
    if !matches!(filter.as_str(), "trimal" | "alifilter" | "none") {
        return Err("--alignment-filter must be trimal, alifilter, or none".into());
    }
    let output = Path::new(&opt.output);
    if method == "concatenation" {
        let alignment = output.join(if filter == "none" {
            "combined_results.fasta"
        } else {
            "combined_trimed.fasta"
        });
        if !alignment.is_file() {
            return Err(format!(
                "Unable to find the concatenated alignment at '{}'",
                alignment.display()
            ));
        }
        let bootstrap = value(raw, &["-b", "--bootstrap"], "1000")?
            .parse::<usize>()
            .map_err(|_| "--bootstrap must be an integer")?;
        let tree = build_tree(&program, &binary, &alignment, bootstrap, opt.workers, false)?;
        fs::copy(tree, output.join("Concatenation.tree")).map_err(|e| e.to_string())?;
        return Ok(());
    }
    let alignment_dir = output.join(if filter == "none" {
        "combined_results/aligned"
    } else {
        "combined_trimed"
    });
    let loci = reference_loci(Path::new(&opt.reference))?;
    let mut trees = Vec::new();
    let mut failures = Vec::new();
    for (locus, _) in loci {
        let alignment = alignment_dir.join(format!("{locus}.fasta"));
        if !alignment.is_file() {
            continue;
        }
        match build_tree(&program, &binary, &alignment, 0, 1, true) {
            Ok(path) => trees.push(path),
            Err(error) => failures.push((locus, alignment, error)),
        }
    }
    let failure_path = output.join("failed_gene_trees.tsv");
    if failures.is_empty() {
        if failure_path.exists() {
            fs::remove_file(&failure_path).map_err(|e| e.to_string())?;
        }
    } else {
        let mut text = "locus\talignment\terror\n".to_owned();
        for (locus, alignment, error) in failures {
            text.push_str(&format!(
                "{locus}\t{}\t{}\n",
                alignment.display(),
                error.replace(['\t', '\n'], " ")
            ));
        }
        fs::write(failure_path, text).map_err(|e| e.to_string())?;
    }
    trees.sort();
    let mut all = String::new();
    for tree in trees {
        let content = fs::read_to_string(&tree).map_err(|e| e.to_string())?;
        if let Some(line) = content.lines().map(str::trim).find(|line| !line.is_empty()) {
            all.push_str(line);
            all.push('\n');
        }
    }
    if all.is_empty() {
        return Err(
            "Unable to reconstruct coalescent trees because no gene tree is available".into(),
        );
    }
    let trees_path = output.join("combined_genes.trees");
    fs::write(&trees_path, all).map_err(|e| e.to_string())?;
    let coalescent = output.join("Coalescent.tree");
    if coalescent.exists() {
        fs::remove_file(&coalescent).map_err(|e| e.to_string())?;
    }
    let astral = PathBuf::from(env::var("GM2_ASTRAL").unwrap_or_else(|_| "astral".into()));
    run_program(
        &astral,
        &[
            "-i".into(),
            trees_path.display().to_string(),
            "-o".into(),
            coalescent.display().to_string(),
            "-t".into(),
            opt.workers.to_string(),
        ],
    )
}

#[derive(Clone)]
struct PermitPool {
    state: Arc<(Mutex<usize>, Condvar)>,
}

struct Permit {
    state: Arc<(Mutex<usize>, Condvar)>,
}

impl PermitPool {
    fn new(limit: usize) -> Self {
        Self {
            state: Arc::new((Mutex::new(limit), Condvar::new())),
        }
    }

    fn acquire(&self) -> Permit {
        let (available, changed) = &*self.state;
        let mut remaining = available.lock().expect("permit pool poisoned");
        while *remaining == 0 {
            remaining = changed.wait(remaining).expect("permit pool poisoned");
        }
        *remaining -= 1;
        Permit {
            state: Arc::clone(&self.state),
        }
    }
}

impl Drop for Permit {
    fn drop(&mut self) {
        let (available, changed) = &*self.state;
        *available.lock().expect("permit pool poisoned") += 1;
        changed.notify_one();
    }
}

fn execute_combine(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    default_source: &str,
) -> Result<(), String> {
    let raw = &opt.raw;
    let source = value(raw, &["-cs", "--combine-source"], default_source)?;
    let input_dir = match source.as_str() {
        "assembly" => "results",
        "consensus" => "consensus",
        "trimmed" => "blast",
        _ => return Err("--combine-source must be assembly, consensus, or trimmed".into()),
    };
    let no_alignment = flag(raw, "--no-alignment")?;
    let filter = if flag(raw, "--no-trimal")? {
        "none".into()
    } else {
        value(raw, &["--alignment-filter"], "trimal")?
    };
    if !matches!(filter.as_str(), "trimal" | "alifilter" | "none") {
        return Err("--alignment-filter must be trimal, alifilter, or none".into());
    }
    let strict = flag(raw, "--strict-combine-errors")?;
    let clean_difference = value(raw, &["-cd", "--clean-difference"], "1")?
        .parse::<f64>()
        .map_err(|_| "--clean-difference must be a number")?;
    let clean_sequences = value(raw, &["-cn", "--clean-sequences"], "0")?
        .parse::<usize>()
        .map_err(|_| "--clean-sequences must be an integer")?;
    if !(0.0..=1.0).contains(&clean_difference) || clean_sequences > samples.len() {
        return Err("invalid combine cleanup thresholds".into());
    }
    let combined = Path::new(&opt.output).join("combined_results");
    if combined.exists() {
        fs::remove_dir_all(&combined).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&combined).map_err(|e| e.to_string())?;
    let loci = reference_loci(Path::new(&opt.reference))?;
    let names = loci
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let uce = opt.assembly_mode == "uce";
    for locus in &names {
        write_combined_locus(locus, input_dir, Path::new(&opt.output), samples, uce)?;
    }
    if no_alignment {
        return Ok(());
    }
    let msa_program = value(raw, &["--msa-program"], "mafft")?;
    if !matches!(msa_program.as_str(), "mafft" | "clustalo") {
        return Err("--msa-program must be mafft or clustalo".into());
    }
    let msa_threads = value(raw, &["--msa-threads"], "1")?
        .parse::<usize>()
        .map_err(|_| "--msa-threads must be an integer")?;
    if msa_threads == 0 {
        return Err("--msa-threads must be at least 1".into());
    }
    if msa_threads > opt.workers {
        return Err("--msa-threads cannot be greater than -p".into());
    }
    let aligned = combined.join("aligned");
    fs::create_dir_all(&aligned).map_err(|e| e.to_string())?;
    let filtered = Path::new(&opt.output).join("combined_trimed");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    if filter != "none" {
        fs::create_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    let mafft = PathBuf::from(env::var("GM2_MAFFT").unwrap_or_else(|_| "mafft".into()));
    let clustalo = PathBuf::from(env::var("GM2_CLUSTALO").unwrap_or_else(|_| "clustalo".into()));
    let trimal = PathBuf::from(env::var("GM2_TRIMAL").unwrap_or_else(|_| "trimal".into()));
    let alifilter = PathBuf::from(env::var("GM2_ALIFILTER").unwrap_or_else(|_| "AliFilter".into()));
    let model = optional_value(raw, &["--alifilter-model"])?;
    if model.is_some() && filter != "alifilter" {
        return Err("--alifilter-model requires --alignment-filter alifilter".into());
    }
    let filter_processes = value(raw, &["--filter-processes"], &opt.workers.to_string())?
        .parse::<usize>()
        .map_err(|_| "--filter-processes must be an integer")?;
    if filter_processes == 0 {
        return Err("--filter-processes must be at least 1".into());
    }
    // Preserve the original scheduler: up to -p loci can make progress, while
    // MSA and column-filter subprocesses have independent resource caps.
    let msa_pool = PermitPool::new((opt.workers / msa_threads).max(1));
    let filter_pool = (filter != "none").then(|| PermitPool::new(filter_processes));
    let pending = Arc::new(Mutex::new(names.into_iter()));
    let failures = Arc::new(Mutex::new(Vec::<(String, String)>::new()));
    thread::scope(|scope| {
        for _ in 0..opt.workers {
            let pending = Arc::clone(&pending);
            let failures = Arc::clone(&failures);
            let combined = combined.clone();
            let aligned = aligned.clone();
            let filtered = filtered.clone();
            let msa_program = msa_program.clone();
            let mafft = mafft.clone();
            let clustalo = clustalo.clone();
            let trimal = trimal.clone();
            let alifilter = alifilter.clone();
            let filter = filter.clone();
            let model = model.clone();
            let msa_pool = msa_pool.clone();
            let filter_pool = filter_pool.clone();
            scope.spawn(move || loop {
                let Some(locus) = pending.lock().expect("combine queue poisoned").next() else {
                    break;
                };
                let input = combined.join(format!("{locus}.fasta"));
                if !input.is_file() {
                    continue;
                }
                let output = aligned.join(format!("{locus}.fasta"));
                let msa_permit = msa_pool.acquire();
                let result = if msa_program == "mafft" {
                    let file = match fs::File::create(&output) {
                        Ok(file) => file,
                        Err(error) => {
                            failures
                                .lock()
                                .expect("combine failures poisoned")
                                .push((locus, error.to_string()));
                            continue;
                        }
                    };
                    let status = Command::new(&mafft)
                        .args([
                            "--auto",
                            "--quiet",
                            "--nuc",
                            "--thread",
                            &msa_threads.to_string(),
                            &input.display().to_string(),
                        ])
                        .stdout(file)
                        .status()
                        .map_err(|error| error.to_string());
                    match status {
                        Ok(status) if status.success() => Ok(()),
                        Ok(status) => Err(format!("mafft exited with {status}")),
                        Err(error) => Err(error),
                    }
                } else {
                    run_program(
                        &clustalo,
                        &[
                            "-i".into(),
                            input.display().to_string(),
                            "-o".into(),
                            output.display().to_string(),
                            "--auto".into(),
                            "--force".into(),
                            "--seqtype=DNA".into(),
                            format!("--threads={msa_threads}"),
                        ],
                    )
                };
                drop(msa_permit);
                let result = result.and_then(|_| {
                    run_program(
                        &bins.join("fix_alignment"),
                        &[
                            "-f".into(),
                            output.display().to_string(),
                            "-n".into(),
                            clean_sequences.to_string(),
                            "-p".into(),
                            clean_difference.to_string(),
                        ],
                    )
                });
                let result = result.and_then(|_| {
                    let filter_permit = filter_pool.as_ref().map(PermitPool::acquire);
                    let filtered_result = if filter == "trimal" {
                        run_program(
                            &trimal,
                            &[
                                "-in".into(),
                                output.display().to_string(),
                                "-out".into(),
                                filtered
                                    .join(format!("{locus}.fasta"))
                                    .display()
                                    .to_string(),
                                "-automated1".into(),
                            ],
                        )
                    } else if filter == "alifilter" {
                        let mut args = vec![
                            "-i".into(),
                            output.display().to_string(),
                            "-o".into(),
                            filtered
                                .join(format!("{locus}.fasta"))
                                .display()
                                .to_string(),
                        ];
                        if let Some(model) = &model {
                            args.extend(["-m".into(), model.clone()]);
                        }
                        run_program(&alifilter, &args)
                    } else {
                        Ok(())
                    };
                    drop(filter_permit);
                    filtered_result
                });
                if let Err(error) = result {
                    let _ = fs::remove_file(&output);
                    failures
                        .lock()
                        .expect("combine failures poisoned")
                        .push((locus, error));
                }
            });
        }
    });
    let mut failures = failures
        .lock()
        .map_err(|_| "combine failures poisoned")?
        .clone();
    failures.sort_by(|left, right| left.0.cmp(&right.0));
    if strict && !failures.is_empty() {
        let (locus, error) = &failures[0];
        return Err(format!("combine failed on {locus}: {error}"));
    }
    for (locus, error) in failures {
        eprintln!("Warning: combine failed on {locus}: {error}");
    }
    run_program(
        &bins.join("merge_seq"),
        &[
            "-input".into(),
            aligned.display().to_string(),
            "-exts".into(),
            ".fasta".into(),
            "-missing".into(),
            "-".into(),
            "-output".into(),
            Path::new(&opt.output)
                .join("combined_results.fasta")
                .display()
                .to_string(),
        ],
    )?;
    if filter != "none" {
        run_program(
            &bins.join("merge_seq"),
            &[
                "-input".into(),
                filtered.display().to_string(),
                "-exts".into(),
                ".fasta".into(),
                "-missing".into(),
                "-".into(),
                "-output".into(),
                Path::new(&opt.output)
                    .join("combined_trimed.fasta")
                    .display()
                    .to_string(),
            ],
        )?;
    }
    Ok(())
}

fn execute_trim(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    default_source: &str,
) -> Result<(), String> {
    let raw = &opt.raw;
    let source = value(raw, &["-ts", "--trim-source"], default_source)?;
    if !matches!(source.as_str(), "assembly" | "consensus") {
        return Err("--trim-source must be assembly or consensus".into());
    }
    let mode_name = value(raw, &["-tm", "--trim-mode"], "terminal")?;
    let mode = match mode_name.as_str() {
        "all" => "0",
        "longest" => "1",
        "terminal" => "2",
        "isoform" => "3",
        _ => return Err("--trim-mode must be all, longest, terminal, or isoform".into()),
    };
    let retention = value(raw, &["-tr", "--trim-retention"], "0")?
        .parse::<f64>()
        .map_err(|_| "--trim-retention must be a number")?;
    if !(0.0..=1.0).contains(&retention) {
        return Err("--trim-retention must be in [0, 1]".into());
    }
    let makeblastdb =
        PathBuf::from(env::var("GM2_MAKEBLASTDB").unwrap_or_else(|_| "makeblastdb".into()));
    let blast = PathBuf::from(if mode_name == "isoform" {
        env::var("GM2_MAGICBLAST").unwrap_or_else(|_| "magicblast".into())
    } else {
        env::var("GM2_BLASTN").unwrap_or_else(|_| "blastn".into())
    });
    let database_dir = Path::new(&opt.output).join("blast_db");
    fs::create_dir_all(&database_dir).map_err(|e| e.to_string())?;
    let loci = reference_loci(Path::new(&opt.reference))?;
    for (locus, reference) in &loci {
        run_program_in(
            &makeblastdb,
            &[
                "-in".into(),
                fs::canonicalize(reference)
                    .map_err(|e| e.to_string())?
                    .display()
                    .to_string(),
                "-dbtype".into(),
                "nucl".into(),
                "-out".into(),
                locus.clone(),
            ],
            &database_dir,
        )?;
    }
    let mut tasks = Vec::new();
    for sample in samples {
        let sample_dir = Path::new(&opt.output).join(&sample.name);
        let input = sample_dir.join(if source == "consensus" {
            "consensus"
        } else {
            "results"
        });
        if !input.is_dir() {
            continue;
        }
        let output = sample_dir.join("blast");
        if output.exists() {
            fs::remove_dir_all(&output).map_err(|e| e.to_string())?;
        }
        fs::create_dir_all(&output).map_err(|e| e.to_string())?;
        for (locus, reference) in &loci {
            let query = input.join(format!("{locus}.fasta"));
            if query.is_file() {
                tasks.push((
                    locus.clone(),
                    query,
                    reference.clone(),
                    output.join(format!("{locus}.fasta")),
                ));
            }
        }
    }
    let trim = bins.join("build_trimed");
    let next = Arc::new(Mutex::new(tasks.into_iter()));
    let failures = Arc::new(Mutex::new(Vec::new()));
    let mut handles = Vec::new();
    for _ in 0..opt.workers {
        let next = Arc::clone(&next);
        let failures = Arc::clone(&failures);
        let trim = trim.clone();
        let blast = blast.clone();
        let database_dir = database_dir.clone();
        handles.push(thread::spawn(move || loop {
            let Some((locus, query, reference, output)) =
                next.lock().expect("trim queue poisoned").next()
            else {
                break;
            };
            let result = run_program(
                &trim,
                &[
                    "-i".into(),
                    query.display().to_string(),
                    "-r".into(),
                    reference.display().to_string(),
                    "-o".into(),
                    output.display().to_string(),
                    "-b".into(),
                    database_dir.join(&locus).display().to_string(),
                    "-m".into(),
                    mode.into(),
                    "-p".into(),
                    (retention * 100.0).to_string(),
                    "--executable".into(),
                    blast.display().to_string(),
                ],
            );
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("trim failures poisoned")
                    .push(format!("{locus}: {error}"));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust trim worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "trim failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} trim task(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_consensus(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    let threshold = value(&opt.raw, &["-c", "--consensus-threshold"], "0.75")?
        .parse::<f64>()
        .map_err(|_| "--consensus-threshold must be a number")?;
    if !(0.0 < threshold && threshold <= 1.0) {
        return Err("--consensus-threshold must be in (0, 1]".into());
    }
    let loci = reference_loci(Path::new(&opt.reference))?;
    let minimap2 = PathBuf::from(env::var("GM2_MINIMAP2").unwrap_or_else(|_| "minimap2".into()));
    let consensus = bins.join("build_consensus");
    let mut tasks = Vec::new();
    for sample in samples {
        let sample_dir = Path::new(&opt.output).join(&sample.name);
        let results = sample_dir.join("results");
        if !results.is_dir() {
            continue;
        }
        let out = sample_dir.join("consensus");
        if out.exists() {
            fs::remove_dir_all(&out).map_err(|e| e.to_string())?;
        }
        fs::create_dir_all(&out).map_err(|e| e.to_string())?;
        let read_extension = fastx_output_extension(&sample.read1);
        for (locus, assembly) in &loci {
            let assembled = results.join(assembly.file_name().ok_or("invalid reference filename")?);
            let reads = sample_dir
                .join("filtered")
                .join(format!("{locus}{read_extension}"));
            if assembled.is_file() && reads.is_file() {
                tasks.push((assembled, reads, out.join(format!("{locus}.sam"))));
            }
        }
    }
    let failures = Arc::new(Mutex::new(Vec::new()));
    let next = Arc::new(Mutex::new(tasks.into_iter()));
    let workers = opt.workers.min(1.max(loci.len() * samples.len()));
    let mut handles = Vec::new();
    for _ in 0..workers {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let minimap2 = minimap2.clone();
        let consensus = consensus.clone();
        handles.push(thread::spawn(move || loop {
            let Some((assembly, reads, sam)) =
                next.lock().expect("consensus queue poisoned").next()
            else {
                break;
            };
            let mapped = run_program(
                &minimap2,
                &[
                    "-ax".into(),
                    "sr".into(),
                    "-t".into(),
                    "1".into(),
                    "--sam-hit-only".into(),
                    "--secondary=no".into(),
                    "-o".into(),
                    sam.display().to_string(),
                    assembly.display().to_string(),
                    reads.display().to_string(),
                ],
            );
            let result = mapped
                .and_then(|_| {
                    run_program(
                        &consensus,
                        &[
                            "-i".into(),
                            sam.display().to_string(),
                            "-c".into(),
                            threshold.to_string(),
                            "-o".into(),
                            sam.parent()
                                .ok_or("SAM has no parent")?
                                .display()
                                .to_string(),
                            "-s".into(),
                            "0".into(),
                        ],
                    )
                })
                .and_then(|_| fs::remove_file(&sam).map_err(|e| e.to_string()));
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("consensus failures poisoned")
                    .push(format!("{}: {error}", sam.display()));
            }
        }));
    }
    for handle in handles {
        handle
            .join()
            .map_err(|_| "Rust consensus worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "consensus failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} consensus task(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_te(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let mut args = vec![
        "--samples".into(),
        opt.samples.clone(),
        "--output".into(),
        opt.output.clone(),
        "--stage".into(),
        value(raw, &["--te-stage"], "all")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--kmer".into(),
        value(raw, &["--te-kmer"], "25")?,
        "--min-kmer-count".into(),
        value(raw, &["--te-min-kmer-count"], "8")?,
        "--catalog-pairs".into(),
        value(raw, &["--te-catalog-pairs"], "10000")?,
        "--mainfilter".into(),
        bins.join("MainFilterNew").display().to_string(),
        "--annotation-min-fragment".into(),
        value(raw, &["--te-annotate-min-fragment"], "80")?,
        "--annotation-max-fragment".into(),
        value(raw, &["--te-annotate-max-fragment"], "800")?,
        "--annotation-min-support".into(),
        value(raw, &["--te-annotate-min-support"], "5")?,
        "--annotation-min-identity".into(),
        value(raw, &["--te-annotate-min-identity"], "0.8")?,
        "--annotation-min-coverage".into(),
        value(raw, &["--te-annotate-min-coverage"], "0.6")?,
        "--annotation-min-delta".into(),
        value(raw, &["--te-annotate-min-delta"], "0.1")?,
        "--assemble-min-kmer-count".into(),
        value(raw, &["--te-assemble-min-kmer-count"], "3")?,
        "--assemble-branch-ratio".into(),
        value(raw, &["--te-assemble-branch-ratio"], "1.5")?,
        "--assemble-max-fragments".into(),
        value(raw, &["--te-assemble-max-fragments"], "3")?,
    ];
    if let Some(path) = optional_value(raw, &["--te-read-ledger"])? {
        args.extend(["--read-ledger".into(), path]);
    }
    if let Some(path) = optional_value(raw, &["--te-library"])? {
        args.extend(["--te-library".into(), path]);
    }
    run(bins, "main_repeat", &args)
}

fn execute_population(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let engine = value(raw, &["--engine"], "pseudoref")?;
    if !matches!(engine.as_str(), "pseudoref" | "panref" | "panrefv2") {
        return Err("--engine must be pseudoref, panref, or panrefv2".into());
    }
    if matches!(engine.as_str(), "panref" | "panrefv2") && opt.reference.is_empty() {
        return Err("-r is required with --engine panref or panrefv2".into());
    }
    let mut args = vec![
        "--output".into(),
        opt.output.clone(),
        "--samples".into(),
        opt.samples.clone(),
        "--engine".into(),
        engine.clone(),
        "--reference-strategy".into(),
        value(raw, &["--population-reference-strategy"], "sqcl-longest")?,
        "--start-at".into(),
        value(raw, &["--population-start-at"], "reference")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--min-mapq".into(),
        value(raw, &["--population-min-mapq"], "20")?,
        "--min-baseq".into(),
        value(raw, &["--population-min-baseq"], "20")?,
        "--min-dp".into(),
        value(raw, &["--population-min-dp"], "5")?,
        "--min-gq".into(),
        value(raw, &["--population-min-gq"], "20")?,
        "--min-qual".into(),
        value(raw, &["--population-min-qual"], "20")?,
        "--min-call-rate".into(),
        value(raw, &["--population-min-call-rate"], "0.8")?,
        "--min-mac".into(),
        value(raw, &["--population-min-mac"], "2")?,
        "--ld-window".into(),
        value(raw, &["--population-ld-window"], "50")?,
        "--ld-step".into(),
        value(raw, &["--population-ld-step"], "5")?,
        "--ld-r2".into(),
        value(raw, &["--population-ld-r2"], "0.2")?,
        "--admixture-k-min".into(),
        value(raw, &["--population-admixture-k-min"], "2")?,
        "--admixture-k-max".into(),
        value(raw, &["--population-admixture-k-max"], "6")?,
        "--admixture-cv".into(),
        value(raw, &["--population-admixture-cv"], "10")?,
        "--stop-after".into(),
        value(raw, &["--population-stop-after"], "selection")?,
        "--minibwa".into(),
        value(raw, &["--population-minibwa"], "minibwa")?,
        "--samtools".into(),
        value(raw, &["--population-samtools"], "samtools")?,
        "--bcftools".into(),
        value(raw, &["--population-bcftools"], "bcftools")?,
        "--plink".into(),
        value(raw, &["--population-plink"], "plink")?,
        "--admixture".into(),
        value(raw, &["--population-admixture"], "admixture")?,
    ];
    if matches!(engine.as_str(), "panref" | "panrefv2") {
        args.extend(["--panref-baits".into(), opt.reference.clone()]);
    }
    if flag(raw, "--population-panrefv2-include-low-confidence")? {
        args.push("--panrefv2-include-low-confidence".into());
    }
    if let Some(path) = optional_value(raw, &["--population-reference-fasta"])? {
        args.extend(["--reference-fasta".into(), path]);
    }
    if flag(raw, "--population-skip-mark-duplicates")? {
        args.push("--skip-mark-duplicates".into());
    }
    if flag(raw, "--population-skip-plink")? {
        args.push("--skip-plink".into());
    }
    if flag(raw, "--population-skip-admixture")? {
        args.push("--skip-admixture".into());
    }
    run(bins, "main_population", &args)
}

fn execute_stats(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    let mut args = vec![
        "--output".into(),
        opt.output.clone(),
        "--reference".into(),
        opt.reference.clone(),
    ];
    for sample in samples {
        args.extend([
            "--sample".into(),
            sample.name.clone(),
            sample.read1.clone(),
            sample.read2.clone().unwrap_or_default(),
        ]);
    }
    if opt.stats_count_input_reads {
        args.push("--count-input-reads".into());
    }
    if opt.stats_no_heatmap {
        args.push("--no-heatmap".into());
    }
    run(bins, "gm2_stats", &args)
}

fn cleanup_native_intermediates(opt: &Options, samples: &[Sample]) -> Result<(), String> {
    if !opt.cleanup_intermediates {
        return Ok(());
    }
    if !opt.commands.iter().any(|command| command == "filter")
        || !opt.commands.iter().any(|command| command == "assemble")
    {
        return Err(
            "--cleanup-intermediates requires filter and assemble in the same invocation".into(),
        );
    }
    let root = fs::canonicalize(&opt.output).map_err(|e| e.to_string())?;
    let mut rows = String::from("path\tbytes\treason\n");
    for sample in samples {
        let sample_dir = root.join(&sample.name);
        for (name, reason) in [
            ("filtered", "reproducible filtered reads"),
            ("filtered_pe", "reproducible filter candidates"),
        ] {
            let path = sample_dir.join(name);
            if path.is_dir() && !path.is_symlink() {
                let bytes = directory_size(&path)?;
                fs::remove_dir_all(&path).map_err(|e| e.to_string())?;
                rows.push_str(&format!("{}\t{bytes}\t{reason}\n", path.display()));
            }
        }
    }
    fs::write(root.join("cleanup_manifest.tsv"), rows).map_err(|e| e.to_string())
}

fn directory_size(path: &Path) -> Result<u64, String> {
    let mut bytes = 0;
    for entry in fs::read_dir(path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let child = entry.path();
        if child.is_symlink() {
            continue;
        }
        if child.is_dir() {
            bytes += directory_size(&child)?;
        } else {
            bytes += fs::metadata(&child).map_err(|e| e.to_string())?.len();
        }
    }
    Ok(bytes)
}

fn write_native_workflow_profile(
    output: &Path,
    profile: &WorkflowProfile,
    elapsed_ms: u128,
) -> Result<(), String> {
    let mut rows = profile
        .rows
        .lock()
        .map_err(|_| "workflow profile poisoned")?
        .clone();
    rows.sort_by(|left, right| {
        (&left.sample, left.round, &left.stage).cmp(&(&right.sample, right.round, &right.stage))
    });
    let mut text =
        String::from("sample\tround\tstage\twall_ms\tinput_bytes\toutput_bytes\tstatus\n");
    for row in rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.sample,
            row.round,
            row.stage,
            row.wall_ms,
            row.input_bytes,
            row.output_bytes,
            row.status
        ));
    }
    text.push_str(&format!(
        "__workflow__\t0\tnative_dispatch\t{elapsed_ms}\t0\t0\tok\n"
    ));
    let path = output.join("workflow_profile.tsv");
    let temporary = output.join("workflow_profile.tsv.tmp");
    fs::write(&temporary, text).map_err(|e| e.to_string())?;
    fs::rename(&temporary, path).map_err(|e| e.to_string())
}

fn validate_parallelism(opt: &Options) -> Result<(), String> {
    let msa_threads = value(&opt.raw, &["--msa-threads"], "1")?
        .parse::<usize>()
        .map_err(|_| "--msa-threads must be an integer")?;
    if msa_threads == 0 {
        return Err("--msa-threads must be at least 1".into());
    }
    if msa_threads > opt.workers {
        return Err("--msa-threads cannot be greater than -p".into());
    }
    if let Some(filter_processes) = optional_value(&opt.raw, &["--filter-processes"])? {
        let filter_processes = filter_processes
            .parse::<usize>()
            .map_err(|_| "--filter-processes must be an integer")?;
        if filter_processes == 0 {
            return Err("--filter-processes must be at least 1".into());
        }
    }
    Ok(())
}

fn execute_native(opt: Options) -> Result<(), String> {
    let workflow_started = Instant::now();
    if opt.output.is_empty() {
        return Err("-o is required".into());
    }
    if opt.workers == 0 {
        return Err("-p must be at least 1".into());
    }
    validate_parallelism(&opt)?;
    let bins = components()?;
    let standalone = [
        "te",
        "gene-annotate",
        "gene-resolve",
        "gene-tree",
        "profiling",
        "mito",
        "rad",
        "rad-probe",
        "rad-validate",
    ];
    if opt.commands.len() > 1
        && opt
            .commands
            .iter()
            .any(|command| standalone.contains(&command.as_str()))
    {
        return Err("this Rust migration route currently requires the selected post-processing command to run alone".into());
    }
    if opt.commands == ["gene-annotate"] {
        return execute_gene_annotate(&opt, &bins);
    }
    if opt.commands == ["gene-resolve"] {
        return execute_gene_resolve(&opt, &bins);
    }
    if opt.commands == ["gene-tree"] {
        return execute_gene_tree(&opt);
    }
    if opt.commands == ["mito"] {
        if opt.samples.is_empty() {
            return Err("-f is required for mito".into());
        }
        let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
        return execute_mito(&opt, &bins, &samples);
    }
    if opt.commands == ["rad-probe"] {
        return execute_rad_probe(&opt, &bins);
    }
    if opt.commands == ["rad-validate"] {
        return execute_rad_validate(&opt, &bins);
    }
    if opt.commands == ["rad"] {
        if opt.samples.is_empty() {
            return Err("-f is required for rad".into());
        }
        return execute_rad(&opt, &bins);
    }
    if opt.commands == ["profiling"] {
        if opt.samples.is_empty() {
            return Err("-f is required for profiling".into());
        }
        let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
        return execute_profiling(&opt, &bins, &samples);
    }
    if opt.samples.is_empty() {
        return Err("-f is required for this command".into());
    }
    if opt.commands == ["te"] {
        return execute_te(&opt, &bins);
    }
    if opt.commands == ["population"] {
        return execute_population(&opt, &bins);
    }
    if opt.reference.is_empty() {
        return Err("-r is required for this command".into());
    }
    fs::create_dir_all(&opt.output).map_err(|e| e.to_string())?;
    let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
    if opt.commands == ["stats"] {
        return execute_stats(&opt, &bins, &samples);
    }
    if opt.commands == ["consensus"] {
        return execute_consensus(&opt, &bins, &samples);
    }
    if opt.commands == ["trim"] {
        return execute_trim(&opt, &bins, &samples, "assembly");
    }
    if opt.commands == ["combine"] {
        return execute_combine(&opt, &bins, &samples, "assembly");
    }
    if opt.commands == ["tree"] {
        return execute_tree(&opt);
    }
    let cohort_samples = samples
        .iter()
        .map(|sample| sample.name.clone())
        .collect::<Vec<_>>();
    let is_uce = opt.assembly_mode == "uce";
    if is_uce {
        let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
        if !matches!(implementation.as_str(), "auto" | "uce-rust") {
            return Err("UCE assembly requires --assembler-implementation auto or uce-rust".into());
        }
    }
    if is_uce && samples.iter().any(|sample| sample.read2.is_none()) {
        return Err(
            "UCE workflow requires paired input; two-column sample lists retain the legacy duplicated-mate convention"
                .into(),
        );
    }
    if is_uce
        && !opt.legacy_uce_filter
        && opt.commands.iter().any(|command| command == "refilter")
        && !opt.commands.iter().any(|command| command == "filter")
    {
        return Err("UCE refilter is fused into the filter stage; run filter first".into());
    }
    let profiler = opt.workflow_profile.then(WorkflowProfile::default);
    let dictionary = reference_dictionary_path(&opt)?;
    if let Some(parent) = dictionary.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    if (!is_uce || opt.legacy_uce_filter) && opt.commands.iter().any(|c| c == "filter") {
        let index_args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-o".into(),
            opt.output.clone(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ];
        let result = run_profiled_action(
            profiler.as_ref(),
            "__reference__",
            "mainfilter_index",
            Path::new(&opt.reference),
            &dictionary,
            || run(&bins, "MainFilterNew", &index_args),
        );
        if let Err(error) = result {
            if let Some(profile) = profiler.as_ref() {
                write_native_workflow_profile(
                    Path::new(&opt.output),
                    profile,
                    workflow_started.elapsed().as_millis(),
                )?;
            }
            return Err(error);
        }
    }
    let failures = Arc::new(Mutex::new(Vec::new()));
    let shared = Arc::new(opt);
    macro_rules! profile_try {
        ($result:expr) => {
            if let Err(error) = $result {
                if let Some(profile) = profiler.as_ref() {
                    write_native_workflow_profile(
                        Path::new(&shared.output),
                        profile,
                        workflow_started.elapsed().as_millis(),
                    )?;
                }
                return Err(error);
            }
        };
    }
    let next = Arc::new(Mutex::new(samples.clone().into_iter()));
    let mut handles = Vec::new();
    for _ in 0..shared.workers {
        let bins = bins.clone();
        let opt = Arc::clone(&shared);
        let next = Arc::clone(&next);
        let failures = Arc::clone(&failures);
        let dictionary = dictionary.clone();
        let profiler = profiler.clone();
        handles.push(thread::spawn(move || loop {
            let Some(sample) = next.lock().expect("sample queue poisoned").next() else {
                break;
            };
            let result = if opt.assembly_mode == "uce" {
                if opt.legacy_uce_filter {
                    execute_uce_legacy(&opt, &bins, &sample, &dictionary, profiler.as_ref())
                } else {
                    execute_uce(&opt, &bins, &sample, profiler.as_ref())
                }
            } else {
                execute_gene(&opt, &bins, &sample, &dictionary, profiler.as_ref())
            };
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("failure list poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust workflow worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "failure list poisoned")?;
    if !failures.is_empty() {
        if let Some(profile) = profiler.as_ref() {
            write_native_workflow_profile(
                Path::new(&shared.output),
                profile,
                workflow_started.elapsed().as_millis(),
            )?;
        }
        return Err(format!(
            "{} sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ));
    }
    if !is_uce && shared.commands.iter().any(|c| c == "gene") {
        let mut cohort = vec![
            "cohort".into(),
            "--reference".into(),
            shared.reference.clone(),
            "--out".into(),
            Path::new(&shared.output).join("gene").display().to_string(),
        ];
        for name in cohort_samples {
            cohort.extend(["--sample".into(), name]);
        }
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__cohort__",
            "gene-cohort",
            Path::new(&shared.output),
            &Path::new(&shared.output).join("gene"),
            || run(&bins, "gene_workflow", &cohort),
        ));
    }
    if shared.commands.iter().any(|command| command == "consensus") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "consensus",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_consensus(&shared, &bins, &samples),
        ));
    }
    if shared.commands.iter().any(|command| command == "trim") {
        let source = if shared.commands.iter().any(|command| command == "consensus") {
            "consensus"
        } else {
            "assembly"
        };
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "trim",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_trim(&shared, &bins, &samples, source),
        ));
    }
    if shared.commands.iter().any(|command| command == "combine") {
        let source = if shared.commands.iter().any(|command| command == "trim") {
            "trimmed"
        } else if shared.commands.iter().any(|command| command == "consensus") {
            "consensus"
        } else {
            "assembly"
        };
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "combine",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_combine(&shared, &bins, &samples, source),
        ));
    }
    if shared.commands.iter().any(|command| command == "tree") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "tree",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_tree(&shared),
        ));
    }
    if shared.commands.iter().any(|command| command == "stats") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "stats",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_stats(&shared, &bins, &samples),
        ));
    }
    if shared
        .commands
        .iter()
        .any(|command| command == "population")
    {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "population",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_population(&shared, &bins),
        ));
    }
    profile_try!(run_profiled_action(
        profiler.as_ref(),
        "__workflow__",
        "cleanup",
        Path::new(&shared.output),
        Path::new(&shared.output),
        || cleanup_native_intermediates(&shared, &samples),
    ));
    if let Some(profile) = profiler.as_ref() {
        write_native_workflow_profile(
            Path::new(&shared.output),
            profile,
            workflow_started.elapsed().as_millis(),
        )?;
    }
    Ok(())
}

fn print_help() {
    println!(
        "GeneMiner2 Rust CLI\n\nNative Rust command dispatcher; no Python runtime is required."
    );
}

fn main() -> ExitCode {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        return ExitCode::SUCCESS;
    }
    match parse(&args).and_then(execute_native) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("Error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn consensus_uses_legacy_filtered_fastx_extension() {
        assert_eq!(fastx_output_extension("reads.fastq.gz"), ".fq");
        assert_eq!(fastx_output_extension("reads.fa"), ".fasta");
        let consensus = parse(&[
            "consensus".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
            "-c".into(),
            "0.9".into(),
        ])
        .unwrap();
        assert_eq!(consensus.commands, ["consensus"]);
    }

    #[test]
    fn te_and_population_are_native_standalone_commands() {
        let te = parse(&[
            "te".into(),
            "-f".into(),
            "taxa.tsv".into(),
            "-o".into(),
            "out".into(),
            "--te-stage".into(),
            "discover".into(),
        ])
        .unwrap();
        assert_eq!(te.commands, ["te"]);
        let population = parse(&[
            "population".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-o".into(),
            "out".into(),
            "--engine".into(),
            "pseudoref".into(),
        ])
        .unwrap();
        assert_eq!(population.commands, ["population"]);
    }

    #[test]
    fn stats_is_a_native_standalone_command() {
        let parsed = parse(&[
            "stats".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
            "--stats-count-input-reads".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["stats"]);
        assert!(parsed.stats_count_input_reads);
    }

    #[test]
    fn gene_expands_to_recovery_stages() {
        let parsed = parse(&[
            "gene".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble", "gene"]);
    }
    #[test]
    fn uce_default_stages_are_complete() {
        let parsed = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap();
        assert_eq!(
            parsed.commands,
            ["filter", "refilter", "assemble", "combine", "tree"]
        );
    }
    #[test]
    fn sample_names_match_legacy_rule() {
        assert_eq!(sample_name("foo bar-1"), "Foo_bar_1");
    }

    #[test]
    fn commands_can_follow_options() {
        let parsed = parse(&[
            "--assembly-mode".into(),
            "original".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
            "gene".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble", "gene"]);
    }

    #[test]
    fn unicode_sample_names_do_not_panic() {
        assert_eq!(sample_name("样本-A"), "样本_a");
    }

    #[test]
    fn equals_style_options_are_parsed() {
        let parsed = parse(&[
            "--assembly-mode=uce".into(),
            "--max-reads=7".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert_eq!(
            parsed.commands,
            ["filter", "refilter", "assemble", "combine", "tree"]
        );
        assert_eq!(parsed.max_reads, "7");
    }

    #[test]
    fn invalid_assembly_mode_is_rejected() {
        let error = parse(&[
            "gene".into(),
            "--assembly-mode".into(),
            "typo".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("--assembly-mode must be original or uce"));
    }

    #[test]
    fn boolean_options_reject_explicit_values() {
        let error = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "--uce-rescue-reads=true".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("--uce-rescue-reads does not take a value"));
    }

    #[test]
    fn unknown_options_are_rejected() {
        let error = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "--max-read=7".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("does not support option '--max-read'"));
    }

    #[test]
    fn incomplete_gene_stage_set_is_rejected() {
        let error = parse(&[
            "gene".into(),
            "filter".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("gene requires filter, refilter, and assemble"));
    }

    #[test]
    fn python_compatibility_options_are_all_accepted() {
        let parsed = parse(&[
            "filter".into(),
            "refilter".into(),
            "assemble".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "--assembler-implementation".into(),
            "uce-rust".into(),
            "--assembler-read-chunk-size".into(),
            "4096".into(),
            "--uce-path-strategy".into(),
            "search".into(),
            "--uce-backbone-lookahead".into(),
            "12".into(),
            "--min-depth".into(),
            "1".into(),
            "--max-depth".into(),
            "2".into(),
            "--reuse-reference-cache".into(),
            "--legacy-uce-filter".into(),
            "--workflow-profile".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert!(parsed.reuse_reference_cache);
        assert!(parsed.legacy_uce_filter);
        assert!(parsed.workflow_profile);
    }

    #[test]
    fn unsupported_stage_is_not_silently_ignored() {
        let error = parse(&[
            "rescue".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap_err();
        assert!(error.contains("does not support command 'rescue'"));
    }
}
