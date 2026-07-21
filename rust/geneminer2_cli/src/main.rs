//! Compatibility-first Rust command dispatcher.
//!
//! Python remains the default command engine during migration.  Setting
//! `GENEMINER2_ENGINE=rust` runs the native UCE and `gene` recovery paths.

use std::env;
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::sync::{Arc, Mutex};
use std::thread;

const COMMANDS: &[&str] = &["filter", "refilter", "assemble", "gene"];
const FLAG_OPTIONS: &[&str] = &["--uce-alignment-shadow", "--uce-rescue-reads"];
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
];

#[derive(Clone, Debug)]
struct Sample {
    name: String,
    read1: String,
    read2: Option<String>,
}

#[derive(Debug)]
struct Options {
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
                return Err(format!(
                    "Rust migration engine does not support option '{option}'; use GENEMINER2_ENGINE=legacy"
                ));
            }
        } else {
            return Err(format!(
                "Rust migration engine does not support command '{arg}'; use GENEMINER2_ENGINE=legacy"
            ));
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
            vec!["filter".into(), "refilter".into(), "assemble".into()]
        } else {
            return Err(
                "the Rust migration engine needs an explicit filter/refilter/assemble command"
                    .into(),
            );
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
        let read2 = fields
            .get(2)
            .filter(|value| !value.is_empty())
            .map(|value| (*value).to_string());
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

fn uce_filter_args(opt: &Options, sample: &Sample, sample_dir: &Path) -> Vec<String> {
    let mut args = vec![
        "-r".into(),
        opt.reference.clone(),
        "--recruit-references".into(),
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
        "--selection".into(),
        "auto".into(),
        "--reference-role".into(),
        "bait".into(),
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
        "backbone".into(),
        "--uce-backbone-lookahead".into(),
        "24".into(),
        "--assembler-read-chunk-size".into(),
        "8192".into(),
        "--assembler-kmer-count-threads".into(),
        "1".into(),
        "--assembler-graph-format".into(),
        opt.graph_format.clone(),
    ])
}

fn execute_uce(opt: &Options, bins: &Path, sample: &Sample) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    if opt.commands.iter().any(|c| c == "filter") {
        run(
            bins,
            "uce_filter",
            &uce_filter_args(opt, sample, &sample_dir),
        )?;
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        run(
            bins,
            "main_assembler-rust",
            &uce_assembler_args(opt, &sample_dir)?,
        )?;
    }
    Ok(())
}

fn execute_gene(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    dictionary: &Path,
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
        run(bins, "MainFilterNew", &args)?;
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
        run(bins, "main_refilter_new", &args)?;
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        run(
            bins,
            "main_assembler-original-rust",
            &[
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
            ],
        )?;
    }
    if opt.commands.iter().any(|c| c == "gene") {
        run(
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
        )?;
    }
    Ok(())
}

fn execute_native(opt: Options) -> Result<(), String> {
    if opt.reference.is_empty() || opt.samples.is_empty() || opt.output.is_empty() {
        return Err("-f, -r and -o are required".into());
    }
    if opt.workers == 0 {
        return Err("-p must be at least 1".into());
    }
    if opt.rescue {
        return Err("Rust migration engine does not yet own --uce-rescue-reads; use GENEMINER2_ENGINE=legacy for that path".into());
    }
    let bins = components()?;
    fs::create_dir_all(&opt.output).map_err(|e| e.to_string())?;
    let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
    let cohort_samples = samples
        .iter()
        .map(|sample| sample.name.clone())
        .collect::<Vec<_>>();
    let is_uce = opt.assembly_mode == "uce";
    if is_uce && samples.iter().any(|sample| sample.read2.is_none()) {
        return Err(
            "Rust UCE workflow currently requires paired-end input; use GENEMINER2_ENGINE=legacy for single-end data"
                .into(),
        );
    }
    if is_uce
        && opt.commands.iter().any(|command| command == "refilter")
        && !opt.commands.iter().any(|command| command == "filter")
    {
        return Err(
            "UCE refilter is fused into the native filter stage; run filter or use GENEMINER2_ENGINE=legacy"
                .into(),
        );
    }
    if !is_uce && !opt.commands.iter().any(|c| c == "gene") {
        return Err(
            "Rust migration engine currently supports --assembly-mode uce and the gene command"
                .into(),
        );
    }
    let dictionary = Path::new(&opt.output).join(format!("kmer_dict_k{}.dict", opt.kf));
    if !is_uce && opt.commands.iter().any(|c| c == "filter") {
        run(
            &bins,
            "MainFilterNew",
            &[
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
            ],
        )?;
    }
    let failures = Arc::new(Mutex::new(Vec::new()));
    let shared = Arc::new(opt);
    let next = Arc::new(Mutex::new(samples.into_iter()));
    let mut handles = Vec::new();
    for _ in 0..shared.workers {
        let bins = bins.clone();
        let opt = Arc::clone(&shared);
        let next = Arc::clone(&next);
        let failures = Arc::clone(&failures);
        let dictionary = dictionary.clone();
        handles.push(thread::spawn(move || loop {
            let Some(sample) = next.lock().expect("sample queue poisoned").next() else {
                break;
            };
            let result = if opt.assembly_mode == "uce" {
                execute_uce(&opt, &bins, &sample)
            } else {
                execute_gene(&opt, &bins, &sample, &dictionary)
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
        run(&bins, "gene_workflow", &cohort)?;
    }
    Ok(())
}

fn legacy(args: &[String]) -> Result<ExitCode, String> {
    let python = env::var("GM2_LEGACY_PYTHON").unwrap_or_else(|_| "python3".into());
    let script = env::var("GM2_LEGACY_SCRIPT").unwrap_or_else(|_| "scripts/unix_command.py".into());
    let status = Command::new(python)
        .arg(script)
        .args(args)
        .status()
        .map_err(|e| e.to_string())?;
    Ok(if status.success() {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(status.code().unwrap_or(1) as u8)
    })
}
fn print_help() {
    println!("GeneMiner2 Rust migration launcher\n\nDefault: GENEMINER2_ENGINE=legacy (existing Python CLI).\nPreview: GENEMINER2_ENGINE=rust for UCE and gene recovery.\n\nAll existing commands and options remain available through the default legacy engine.");
}
fn main() -> ExitCode {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|arg| arg == "-h" || arg == "--help")
        && env::var("GENEMINER2_ENGINE").ok().as_deref() == Some("rust")
    {
        print_help();
        return ExitCode::SUCCESS;
    }
    match env::var("GENEMINER2_ENGINE")
        .unwrap_or_else(|_| "legacy".into())
        .as_str()
    {
        "legacy" => legacy(&args).unwrap_or_else(|error| {
            eprintln!("Error: {error}");
            ExitCode::from(1)
        }),
        "rust" => parse(&args)
            .and_then(execute_native)
            .map(|_| ExitCode::SUCCESS)
            .unwrap_or_else(|error| {
                eprintln!("Error: {error}");
                ExitCode::from(1)
            }),
        other => {
            eprintln!("Error: GENEMINER2_ENGINE must be legacy or rust, got '{other}'");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble"]);
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
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble"]);
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
    fn unsupported_stage_is_not_silently_ignored() {
        let error = parse(&[
            "stats".into(),
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
        assert!(error.contains("does not support command 'stats'"));
    }
}
