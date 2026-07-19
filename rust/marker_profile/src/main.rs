use std::collections::{BTreeMap, BTreeSet};

use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;

type CountMap = BTreeMap<String, usize>;

#[derive(Debug)]
struct SupportCounts {
    evidence: CountMap,
    exclusive: CountMap,
    queries: usize,
    positive: usize,
    target_decoy_shared: usize,
}

#[derive(Debug)]
struct Args {
    reference: PathBuf,
    reads: PathBuf,
    output: PathBuf,
    cache: PathBuf,
    themisto: PathBuf,
    msweep: PathBuf,
    groups: PathBuf,
    decoy: Option<PathBuf>,
    threads: usize,
    kmer_size: usize,
    threshold: f64,
    relevant_kmer_fraction: f64,
    index_memory_gb: usize,
    min_evidence: usize,
    force_rebuild: bool,
}

struct CacheLock(PathBuf);

impl CacheLock {
    fn acquire(path: PathBuf) -> Result<Self, String> {
        for _ in 0..600 {
            match fs::create_dir(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {
                    thread::sleep(Duration::from_millis(100))
                }
                Err(error) => {
                    return Err(format!(
                        "cannot acquire marker reference-cache lock: {error}"
                    ))
                }
            }
        }
        Err("timed out waiting for marker reference-cache construction".to_string())
    }
}

impl Drop for CacheLock {
    fn drop(&mut self) {
        let _ = fs::remove_dir(&self.0);
    }
}

fn usage() -> &'static str {
    "Usage: marker_profile --reference REF.fa --reads recruited.fq --output DIR --cache DIR \\
  --themisto PATH --msweep PATH [options]\n\n\
Options:\n\
  --groups FILE                    Required TSV: reference ID, then reporting group\n\
  --decoy FILE                     Optional non-target decoy FASTA\n\
  --threads INT                    Worker threads (default: 1)\n\
  --kmer-size INT                  Themisto k-mer size (default: 21)\n\
  --threshold FLOAT                Themisto pseudoalignment threshold (default: 0.80)\n\
  --relevant-kmer-fraction FLOAT   Minimum fraction of query k-mers found in any target (default: 0.50)\n\
  --index-memory-gb INT            Themisto build memory limit (default: 2)\n\
  --min-evidence INT               Minimum exclusive group-supporting queries required for detection (default: 3)\n\
  --force-rebuild                  Rebuild the cached reference index\n"
}

fn next_value(args: &[String], pos: &mut usize, option: &str) -> Result<String, String> {
    *pos += 1;
    args.get(*pos)
        .cloned()
        .ok_or_else(|| format!("{option} requires a value"))
}

fn parse_args() -> Result<Args, String> {
    let raw: Vec<String> = env::args().skip(1).collect();
    if raw.iter().any(|arg| arg == "--help" || arg == "-h") {
        return Err(usage().to_string());
    }
    let mut reference = None;
    let mut reads = None;
    let mut output = None;
    let mut cache = None;
    let mut themisto = None;
    let mut msweep = None;
    let mut groups = None;
    let mut decoy = None;
    let mut threads = 1usize;
    let mut kmer_size = 21usize;
    let mut threshold = 0.80f64;
    let mut relevant_kmer_fraction = 0.50f64;
    let mut index_memory_gb = 2usize;
    let mut min_evidence = 3usize;
    let mut force_rebuild = false;
    let mut pos = 0;
    while pos < raw.len() {
        match raw[pos].as_str() {
            "--reference" => {
                reference = Some(PathBuf::from(next_value(&raw, &mut pos, "--reference")?))
            }
            "--reads" => reads = Some(PathBuf::from(next_value(&raw, &mut pos, "--reads")?)),
            "--output" => output = Some(PathBuf::from(next_value(&raw, &mut pos, "--output")?)),
            "--cache" => cache = Some(PathBuf::from(next_value(&raw, &mut pos, "--cache")?)),
            "--themisto" => {
                themisto = Some(PathBuf::from(next_value(&raw, &mut pos, "--themisto")?))
            }
            "--msweep" => msweep = Some(PathBuf::from(next_value(&raw, &mut pos, "--msweep")?)),
            "--groups" => groups = Some(PathBuf::from(next_value(&raw, &mut pos, "--groups")?)),
            "--decoy" => decoy = Some(PathBuf::from(next_value(&raw, &mut pos, "--decoy")?)),
            "--threads" => {
                threads = next_value(&raw, &mut pos, "--threads")?
                    .parse()
                    .map_err(|_| "invalid --threads".to_string())?
            }
            "--kmer-size" => {
                kmer_size = next_value(&raw, &mut pos, "--kmer-size")?
                    .parse()
                    .map_err(|_| "invalid --kmer-size".to_string())?
            }
            "--threshold" => {
                threshold = next_value(&raw, &mut pos, "--threshold")?
                    .parse()
                    .map_err(|_| "invalid --threshold".to_string())?
            }
            "--relevant-kmer-fraction" => {
                relevant_kmer_fraction = next_value(&raw, &mut pos, "--relevant-kmer-fraction")?
                    .parse()
                    .map_err(|_| "invalid --relevant-kmer-fraction".to_string())?
            }
            "--index-memory-gb" => {
                index_memory_gb = next_value(&raw, &mut pos, "--index-memory-gb")?
                    .parse()
                    .map_err(|_| "invalid --index-memory-gb".to_string())?
            }
            "--min-evidence" => {
                min_evidence = next_value(&raw, &mut pos, "--min-evidence")?
                    .parse()
                    .map_err(|_| "invalid --min-evidence".to_string())?
            }
            "--force-rebuild" => force_rebuild = true,
            other => return Err(format!("unknown option: {other}\n\n{}", usage())),
        }
        pos += 1;
    }
    let args = Args {
        reference: reference.ok_or_else(|| "--reference is required".to_string())?,
        reads: reads.ok_or_else(|| "--reads is required".to_string())?,
        output: output.ok_or_else(|| "--output is required".to_string())?,
        cache: cache.ok_or_else(|| "--cache is required".to_string())?,
        themisto: themisto.ok_or_else(|| "--themisto is required".to_string())?,
        msweep: msweep.ok_or_else(|| "--msweep is required".to_string())?,
        groups: groups.ok_or_else(|| "--groups is required".to_string())?,
        decoy,
        threads,
        kmer_size,
        threshold,
        relevant_kmer_fraction,
        index_memory_gb,
        min_evidence,
        force_rebuild,
    };
    if args.threads == 0
        || args.kmer_size < 15
        || args.kmer_size > 31
        || args.kmer_size.is_multiple_of(2)
        || args.min_evidence == 0
        || !(0.0..=1.0).contains(&args.threshold)
        || !(0.0..=1.0).contains(&args.relevant_kmer_fraction)
    {
        return Err("invalid marker quantification parameter".to_string());
    }
    Ok(args)
}

fn load_group_map(path: &Path) -> Result<BTreeMap<String, String>, String> {
    let input =
        File::open(path).map_err(|e| format!("cannot open group map {}: {e}", path.display()))?;
    let mut map = BTreeMap::new();
    for (line_no, line) in BufReader::new(input).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.is_empty() || line.starts_with("#") {
            continue;
        }
        let mut fields = line.split("\t");
        let id = fields.next().unwrap_or("").trim();
        let group = fields.next().unwrap_or("").trim();
        if id.is_empty() || group.is_empty() || fields.next().is_some() {
            return Err(format!(
                "invalid group-map line {}: expected reference_id<TAB>group",
                line_no + 1
            ));
        }
        if let Some(existing) = map.get(id) {
            if existing != group {
                return Err(format!("conflicting groups for reference ID: {id}"));
            }
        } else {
            map.insert(id.to_string(), group.to_string());
        }
    }
    if map.is_empty() {
        return Err("group map is empty".to_string());
    }
    Ok(map)
}

fn validate_group_map_coverage(
    reference: &Path,
    map: &BTreeMap<String, String>,
) -> Result<(), String> {
    let input =
        File::open(reference).map_err(|e| format!("cannot open {}: {e}", reference.display()))?;
    let mut reference_ids = BTreeSet::new();
    for line in BufReader::new(input).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = header.split_whitespace().next() {
                reference_ids.insert(id.to_string());
            }
        }
    }
    let mapped_ids: BTreeSet<String> = map.keys().cloned().collect();
    let missing: Vec<_> = reference_ids
        .difference(&mapped_ids)
        .take(10)
        .cloned()
        .collect();
    let unused: Vec<_> = mapped_ids
        .difference(&reference_ids)
        .take(10)
        .cloned()
        .collect();
    if !missing.is_empty() || !unused.is_empty() {
        return Err(format!(
            "group-map coverage mismatch; missing reference IDs: {:?}; unused mapped IDs: {:?}",
            missing, unused
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_record(
    dir: &Path,
    id: usize,
    header: &str,
    sequence: &str,
    group: &str,
    list: &mut BufWriter<File>,
    groups: &mut BufWriter<File>,
    metadata: &mut BufWriter<File>,
) -> io::Result<()> {
    let record_id = format!("marker_{id:06}");
    let record_path = dir.join(format!("{record_id}.fasta"));
    let mut record = BufWriter::new(File::create(&record_path)?);
    writeln!(record, ">{record_id}")?;
    writeln!(record, "{}", sequence)?;
    writeln!(list, "{}", record_path.display())?;
    writeln!(groups, "{group}")?;
    writeln!(
        metadata,
        "{id}\t{record_id}\t{group}\t{}",
        header.replace('\t', " ")
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_fasta(
    path: &Path,
    records_dir: &Path,
    next_id: &mut usize,
    list: &mut BufWriter<File>,
    groups: &mut BufWriter<File>,
    metadata: &mut BufWriter<File>,
    group_map: Option<&BTreeMap<String, String>>,
    decoy: bool,
) -> Result<(usize, usize), String> {
    let input = File::open(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let mut header = String::new();
    let mut sequence = String::new();
    let mut kept = 0usize;
    let skipped = 0usize;
    let mut flush = |header: &str, sequence: &str| -> Result<(), String> {
        if header.is_empty() || sequence.is_empty() {
            return Ok(());
        }
        let group = if decoy {
            "DECOY".to_string()
        } else {
            let id = header.split_whitespace().next().unwrap_or("");
            group_map
                .and_then(|map| map.get(id))
                .cloned()
                .ok_or_else(|| format!("reference ID {id} is missing from --groups"))?
        };
        write_record(
            records_dir,
            *next_id,
            header,
            sequence,
            &group,
            list,
            groups,
            metadata,
        )
        .map_err(|e| e.to_string())?;
        *next_id += 1;
        kept += 1;
        Ok(())
    };
    for line in BufReader::new(input).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(next_header) = line.strip_prefix('>') {
            flush(&header, &sequence)?;
            header = next_header.trim().to_string();
            sequence.clear();
        } else {
            sequence.push_str(line.trim());
        }
    }
    flush(&header, &sequence)?;
    Ok((kept, skipped))
}

fn run(command: &mut Command, label: &str) -> Result<(), String> {
    let status = command
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .map_err(|e| format!("cannot launch {label}: {e}"))?;
    status
        .success()
        .then_some(())
        .ok_or_else(|| format!("{label} failed with {status}"))
}

fn prepare_reference(args: &Args) -> Result<(PathBuf, PathBuf, PathBuf), String> {
    let _lock = CacheLock::acquire(args.cache.with_extension("lock"))?;
    let records = args.cache.join("records");
    let list_path = args.cache.join("reference_files.txt");
    let groups_path = args.cache.join("groups.txt");
    let metadata_path = args.cache.join("marker_reference_metadata.tsv");
    let index_prefix = args.cache.join("themisto_index");
    let index_ready = index_prefix.with_extension("tdbg").is_file()
        && index_prefix.with_extension("tcolors").is_file();
    if !args.force_rebuild
        && index_ready
        && list_path.is_file()
        && groups_path.is_file()
        && metadata_path.is_file()
    {
        return Ok((index_prefix, groups_path, metadata_path));
    }
    if args.cache.exists() {
        fs::remove_dir_all(&args.cache).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&records).map_err(|e| e.to_string())?;
    let mut list = BufWriter::new(File::create(&list_path).map_err(|e| e.to_string())?);
    let mut groups = BufWriter::new(File::create(&groups_path).map_err(|e| e.to_string())?);
    let mut metadata = BufWriter::new(File::create(&metadata_path).map_err(|e| e.to_string())?);
    writeln!(metadata, "color\treference_id\tgroup\toriginal_header").map_err(|e| e.to_string())?;
    let mut next_id = 0usize;
    let group_map = load_group_map(&args.groups)?;
    validate_group_map_coverage(&args.reference, &group_map)?;
    let (kept, _) = append_fasta(
        &args.reference,
        &records,
        &mut next_id,
        &mut list,
        &mut groups,
        &mut metadata,
        Some(&group_map),
        false,
    )?;
    if let Some(decoy) = &args.decoy {
        append_fasta(
            decoy,
            &records,
            &mut next_id,
            &mut list,
            &mut groups,
            &mut metadata,
            None,
            true,
        )?;
    }
    list.flush().map_err(|e| e.to_string())?;
    groups.flush().map_err(|e| e.to_string())?;
    metadata.flush().map_err(|e| e.to_string())?;
    if kept == 0 {
        return Err("no marker reference records were loaded".to_string());
    }
    let temp_dir = args.cache.join("temp");
    fs::create_dir_all(&temp_dir).map_err(|e| e.to_string())?;
    let mut command = Command::new(&args.themisto);
    command
        .args(["build", "-k", &args.kmer_size.to_string(), "-i"])
        .arg(&list_path)
        .args(["--index-prefix"])
        .arg(&index_prefix)
        .args(["--temp-dir"])
        .arg(&temp_dir)
        .args([
            "--mem-gigas",
            &args.index_memory_gb.to_string(),
            "--n-threads",
            &args.threads.to_string(),
            "--file-colors",
        ]);
    run(&mut command, "Themisto index construction")?;
    Ok((index_prefix, groups_path, metadata_path))
}

fn support_counts(pseudoalignments: &Path, groups_file: &Path) -> Result<SupportCounts, String> {
    let groups: Vec<String> = BufReader::new(File::open(groups_file).map_err(|e| e.to_string())?)
        .lines()
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
    let mut evidence = BTreeMap::<String, usize>::new();
    let mut exclusive = BTreeMap::<String, usize>::new();
    let mut queries = 0usize;
    let mut positive = 0usize;
    let mut target_decoy_shared = 0usize;
    for line in BufReader::new(File::open(pseudoalignments).map_err(|e| e.to_string())?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        queries += 1;
        let mut labels = Vec::<String>::new();
        for field in line.split_whitespace().skip(1) {
            let color: usize = field
                .parse()
                .map_err(|_| format!("invalid Themisto color: {field}"))?;
            let label = groups
                .get(color)
                .ok_or_else(|| format!("Themisto color {color} has no group label"))?
                .clone();
            if !labels.contains(&label) {
                labels.push(label);
            }
        }
        if !labels.is_empty() {
            positive += 1;
        }
        if labels.len() > 1 && labels.iter().any(|label| label == "DECOY") {
            target_decoy_shared += 1;
        }
        for label in &labels {
            *evidence.entry(label.clone()).or_insert(0) += 1;
        }
        if labels.len() == 1 {
            *exclusive.entry(labels[0].clone()).or_insert(0) += 1;
        }
    }
    Ok(SupportCounts {
        evidence,
        exclusive,
        queries,
        positive,
        target_decoy_shared,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_results(
    output: &Path,
    abundance_file: &Path,
    pseudoalignments: &Path,
    groups_file: &Path,
    metadata: &Path,
    min_evidence: usize,
    kmer_size: usize,
    threshold: f64,
    relevant_kmer_fraction: f64,
) -> Result<(), String> {
    let mut abundance = BTreeMap::<String, f64>::new();
    let mut msweep_reads = None;
    let mut msweep_aligned = None;
    for line in BufReader::new(File::open(abundance_file).map_err(|e| e.to_string())?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(value) = line.strip_prefix("#num_reads:\t") {
            msweep_reads = value.parse::<usize>().ok();
        }
        if let Some(value) = line.strip_prefix("#num_aligned:\t") {
            msweep_aligned = value.parse::<usize>().ok();
        }
        if line.starts_with('#') || line.starts_with("#c_id") || line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        if let (Some(group), Some(value)) = (fields.next(), fields.next()) {
            abundance.insert(
                group.to_string(),
                value
                    .parse()
                    .map_err(|_| format!("invalid mSWEEP abundance: {line}"))?,
            );
        }
    }
    let SupportCounts {
        evidence,
        exclusive,
        queries: pseudo_lines,
        positive: pseudo_positive,
        target_decoy_shared,
    } = support_counts(pseudoalignments, groups_file)?;
    let groups: BTreeSet<String> =
        BufReader::new(File::open(groups_file).map_err(|e| e.to_string())?)
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| e.to_string())?
            .into_iter()
            .filter(|group| group != "DECOY")
            .collect();
    let total: f64 = groups
        .iter()
        .map(|group| {
            if exclusive.get(group).copied().unwrap_or(0) >= min_evidence {
                abundance.get(group).copied().unwrap_or(0.0)
            } else {
                0.0
            }
        })
        .sum();
    let mut out = BufWriter::new(
        File::create(output.join("marker_group_abundance.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(out, "group\traw_abundance\tevidence_queries\texclusive_queries\tdetection_status\trelative_proportion\tcalibration_status").map_err(|e| e.to_string())?;
    for label in &groups {
        let value = abundance.get(label).copied().unwrap_or(0.0);
        let exclusive_count = exclusive.get(label).copied().unwrap_or(0);
        let detected = exclusive_count >= min_evidence;
        let proportion = if detected && total > 0.0 {
            value / total
        } else {
            0.0
        };
        writeln!(
            out,
            "{label}\t{value:.8}\t{}\t{}\t{}\t{proportion:.8}\tuncalibrated",
            evidence.get(label).copied().unwrap_or(0),
            exclusive_count,
            if detected { "detected" } else { "not_detected" }
        )
        .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())?;
    let mut qc =
        BufWriter::new(File::create(output.join("marker_qc.tsv")).map_err(|e| e.to_string())?);
    writeln!(qc, "metric\tvalue").map_err(|e| e.to_string())?;
    writeln!(qc, "pseudoaligned_queries\t{pseudo_lines}").map_err(|e| e.to_string())?;
    writeln!(qc, "queries_with_reference_hits\t{pseudo_positive}").map_err(|e| e.to_string())?;
    writeln!(qc, "target_decoy_shared_queries\t{target_decoy_shared}")
        .map_err(|e| e.to_string())?;
    writeln!(
        qc,
        "decoy_evidence_queries\t{}",
        evidence.get("DECOY").copied().unwrap_or(0)
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        qc,
        "decoy_exclusive_queries\t{}",
        exclusive.get("DECOY").copied().unwrap_or(0)
    )
    .map_err(|e| e.to_string())?;
    writeln!(qc, "kmer_size\t{kmer_size}").map_err(|e| e.to_string())?;
    writeln!(qc, "pseudoalign_threshold\t{threshold}").map_err(|e| e.to_string())?;
    writeln!(qc, "relevant_kmer_fraction\t{relevant_kmer_fraction}").map_err(|e| e.to_string())?;
    writeln!(qc, "min_exclusive_evidence\t{min_evidence}").map_err(|e| e.to_string())?;
    if let Some(value) = msweep_reads {
        writeln!(qc, "msweep_queries\t{value}").map_err(|e| e.to_string())?;
    }
    if let Some(value) = msweep_aligned {
        writeln!(qc, "msweep_aligned_queries\t{value}").map_err(|e| e.to_string())?;
    }
    if let Some(value) = abundance.get("DECOY") {
        writeln!(qc, "decoy_abundance\t{value:.8}").map_err(|e| e.to_string())?;
    }
    qc.flush().map_err(|e| e.to_string())?;
    fs::copy(metadata, output.join("marker_reference_metadata.tsv")).map_err(|e| e.to_string())?;
    Ok(())
}

fn main() -> Result<(), String> {
    if env::args()
        .skip(1)
        .any(|arg| arg == "--help" || arg == "-h")
    {
        print!("{}", usage());
        return Ok(());
    }
    let args = parse_args()?;
    if !args.themisto.is_file()
        || !args.msweep.is_file()
        || !args.reference.is_file()
        || !args.reads.is_file()
        || !args.groups.is_file()
    {
        return Err(
            "reference, group map, reads, Themisto and mSWEEP paths must be existing files"
                .to_string(),
        );
    }
    let (index_prefix, groups, metadata) = prepare_reference(&args)?;
    if args.output.exists()
        && fs::read_dir(&args.output)
            .map_err(|e| e.to_string())?
            .next()
            .is_some()
    {
        return Err(format!(
            "output directory is not empty: {}",
            args.output.display()
        ));
    }
    fs::create_dir_all(&args.output).map_err(|e| e.to_string())?;
    let pseudoalignments = args.output.join("themisto_pseudoalignments.txt");
    let temp_dir = args.output.join("themisto_temp");
    fs::create_dir_all(&temp_dir).map_err(|e| e.to_string())?;
    let mut pseudoalign = Command::new(&args.themisto);
    pseudoalign
        .args(["pseudoalign", "--query-file"])
        .arg(&args.reads)
        .args(["--index-prefix"])
        .arg(&index_prefix)
        .args(["--temp-dir"])
        .arg(&temp_dir)
        .args(["--out-file"])
        .arg(&pseudoalignments)
        .args([
            "--sort-output-lines",
            "--threshold",
            &args.threshold.to_string(),
            "--relevant-kmers-fraction",
            &args.relevant_kmer_fraction.to_string(),
            "--n-threads",
            &args.threads.to_string(),
        ]);
    run(&mut pseudoalign, "Themisto pseudoalignment")?;
    let prefix = args.output.join("msweep");
    let mut msweep = Command::new(&args.msweep);
    msweep
        .args(["--themisto"])
        .arg(&pseudoalignments)
        .args(["-i"])
        .arg(&groups)
        .args(["-o"])
        .arg(&prefix)
        .args(["-t", &args.threads.to_string()]);
    run(&mut msweep, "mSWEEP abundance estimation")?;
    write_results(
        &args.output,
        &args.output.join("msweep_abundances.txt"),
        &pseudoalignments,
        &groups,
        &metadata,
        args.min_evidence,
        args.kmer_size,
        args.threshold,
        args.relevant_kmer_fraction,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_file(name: &str, contents: &str) -> PathBuf {
        let path = env::temp_dir().join(format!(
            "marker_profile_{}_{}_{}",
            std::process::id(),
            name,
            contents.len()
        ));
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn group_map_accepts_consistent_duplicates() {
        let path = temp_file("duplicates.tsv", "r1\tA\nr1\tA\nr2\tB\n");
        let map = load_group_map(&path).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("r1").map(String::as_str), Some("A"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn group_map_rejects_conflicting_duplicates() {
        let path = temp_file("conflict.tsv", "r1\tA\nr1\tB\n");
        assert!(load_group_map(&path)
            .unwrap_err()
            .contains("conflicting groups"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn group_map_requires_exact_reference_coverage() {
        let reference = temp_file("reference.fa", ">r1 first\nACGT\n>r2\nACGT\n");
        let mut map = BTreeMap::new();
        map.insert("r1".to_string(), "A".to_string());
        map.insert("unused".to_string(), "B".to_string());
        let error = validate_group_map_coverage(&reference, &map).unwrap_err();
        assert!(error.contains("r2"));
        assert!(error.contains("unused"));
        let _ = fs::remove_file(reference);
    }

    #[test]
    fn results_use_dynamic_group_names() {
        let root = env::temp_dir().join(format!("marker_profile_results_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let abundance = root.join("abundance.tsv");
        let pseudo = root.join("pseudo.txt");
        let groups = root.join("groups.txt");
        let metadata = root.join("metadata.tsv");
        fs::write(&abundance, "Alpha\t0.6\nBeta\t0.4\n").unwrap();
        fs::write(&pseudo, "0 0\n1 1\n2 0 1\n").unwrap();
        fs::write(&groups, "Alpha\nBeta\n").unwrap();
        fs::write(&metadata, "color\treference_id\tgroup\toriginal_header\n").unwrap();
        write_results(
            &root, &abundance, &pseudo, &groups, &metadata, 1, 21, 0.8, 0.5,
        )
        .unwrap();
        let result = fs::read_to_string(root.join("marker_group_abundance.tsv")).unwrap();
        assert!(result.contains("Alpha\t"));
        assert!(result.contains("Beta\t"));
        let _ = fs::remove_dir_all(root);
    }
}
