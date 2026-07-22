//! RAD matrix augmentation helpers for GeneMiner2-UCE.
//!
//! The workflow deliberately treats the two sequenced RAD arms as independent
//! observations.  It never invents the unsequenced insert between them.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use gm2_tools::fastx::{FastxFormat, FastxReader};

type ArmReferences = (Vec<(String, String)>, Vec<(String, String)>);

#[derive(Clone)]
struct ArmRecord {
    locus: String,
    raw_locus: String,
    sample: String,
    r1: String,
    r2: String,
}

fn usage() -> ! {
    eprintln!("Usage:\n  rad_workflow reference --loci FILE --out DIR\n  rad_workflow denovo --out DIR --sample NAME --read1 FILE --read2 FILE [--sample NAME --read1 FILE --read2 FILE ...] [--overhang DNA] [--overhang-r2 DNA] [--kmer N] [--min-count N] [--min-samples N] [--min-length N]\n  rad_workflow validate --reference DIR --recovery DIR --out DIR [--sample NAME ...] [--min-identity FLOAT] [--min-breadth FLOAT] [--min-delta FLOAT]\n  rad_workflow finalize --reference DIR --recovery DIR --out DIR --sample NAME [--sample NAME ...] [--min-arm-breadth FLOAT]");
    std::process::exit(2);
}

fn options(args: &[String]) -> HashMap<String, Vec<String>> {
    let mut values = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let flag = &args[index];
        if !flag.starts_with("--") || index + 1 >= args.len() {
            usage();
        }
        values
            .entry(flag.clone())
            .or_insert_with(Vec::new)
            .push(args[index + 1].clone());
        index += 2;
    }
    values
}

fn required_path(values: &HashMap<String, Vec<String>>, flag: &str) -> PathBuf {
    values
        .get(flag)
        .and_then(|v| v.first())
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            eprintln!("Missing {flag}");
            usage();
        })
}

fn required_values(values: &HashMap<String, Vec<String>>, flag: &str) -> Vec<String> {
    values
        .get(flag)
        .cloned()
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| {
            eprintln!("Missing {flag}");
            usage();
        })
}

fn safe_name(raw: &str) -> String {
    let value = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>();
    if value.is_empty() {
        "unnamed".into()
    } else {
        value
    }
}

fn normalize(sequence: &str) -> String {
    sequence
        .bytes()
        .filter_map(|base| match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => Some(base.to_ascii_uppercase() as char),
            b'U' => Some('T'),
            _ => None,
        })
        .collect()
}

fn split_arms(sequence: &str) -> Option<(String, String)> {
    let raw = sequence.as_bytes();
    let mut index = 0;
    while index + 2 < raw.len() {
        if raw[index].eq_ignore_ascii_case(&b'N')
            && raw[index + 1].eq_ignore_ascii_case(&b'N')
            && raw[index + 2].eq_ignore_ascii_case(&b'N')
        {
            let start = index;
            index += 3;
            while index < raw.len()
                && (raw[index].eq_ignore_ascii_case(&b'N') || raw[index] == b'-')
            {
                index += 1;
            }
            let left = normalize(&sequence[..start]);
            let right = normalize(&sequence[index..]);
            return (!left.is_empty() && !right.is_empty()).then_some((left, right));
        }
        index += 1;
    }
    None
}

fn locus_id(separator: &str, ordinal: usize) -> String {
    let bars = separator
        .match_indices('|')
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    if bars.len() >= 2 {
        return safe_name(&separator[bars[0] + 1..bars[bars.len() - 1]]);
    }
    format!("locus_{ordinal}")
}

fn parse_loci(path: &Path) -> io::Result<Vec<ArmRecord>> {
    let reader = BufReader::new(File::open(path)?);
    let mut records = Vec::new();
    let mut pending = Vec::<(String, String)>::new();
    let mut ordinal = 0;
    let mut flush = |separator: Option<&str>, pending: &mut Vec<(String, String)>| {
        if pending.is_empty() {
            return;
        }
        ordinal += 1;
        let raw_locus = separator
            .map(|line| locus_id(line, ordinal))
            .unwrap_or_else(|| format!("locus_{ordinal}"));
        let locus = safe_name(&raw_locus);
        for (sample, sequence) in pending.drain(..) {
            if let Some((r1, r2)) = split_arms(&sequence) {
                records.push(ArmRecord {
                    locus: locus.clone(),
                    raw_locus: raw_locus.clone(),
                    sample: safe_name(&sample),
                    r1,
                    r2,
                });
            }
        }
    };
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed.starts_with('/') {
            flush(Some(trimmed), &mut pending);
            continue;
        }
        let fields = trimmed.split_whitespace().collect::<Vec<_>>();
        if fields.len() >= 2 {
            let sample = fields[0].split(".trimmed").next().unwrap_or(fields[0]);
            pending.push((sample.into(), fields[fields.len() - 1].into()));
        }
    }
    flush(None, &mut pending);
    if records.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no .loci records with an NNN arm separator",
        ));
    }
    let mut normalized_loci = BTreeMap::new();
    for record in &records {
        if let Some(previous) =
            normalized_loci.insert(record.locus.clone(), record.raw_locus.clone())
        {
            if previous != record.raw_locus {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "distinct locus names collide after safe-name normalization: {previous} and {}",
                        record.raw_locus
                    ),
                ));
            }
        }
    }
    let mut names = BTreeSet::new();
    for record in &records {
        if !names.insert((record.locus.clone(), record.sample.clone())) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "duplicate locus/sample after safe-name normalization: {}/{}",
                    record.locus, record.sample
                ),
            ));
        }
    }
    Ok(records)
}

fn write_fasta(path: &Path, records: &[(String, String)]) -> io::Result<()> {
    let mut out = File::create(path)?;
    for (id, sequence) in records {
        writeln!(out, ">{id}\n{sequence}")?;
    }
    Ok(())
}

fn reference(values: &HashMap<String, Vec<String>>) -> io::Result<()> {
    let loci = required_path(values, "--loci");
    let out = required_path(values, "--out");
    let arms = out.join("arms");
    if out.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("reference output already exists: {}", out.display()),
        ));
    }
    let records = parse_loci(&loci)?;
    fs::create_dir_all(&arms)?;
    let mut grouped = BTreeMap::<String, Vec<ArmRecord>>::new();
    for record in records {
        grouped
            .entry(record.locus.clone())
            .or_default()
            .push(record);
    }
    let mut manifest = File::create(out.join("locus_manifest.tsv"))?;
    writeln!(
        manifest,
        "locus\toriginal_samples\tr1_max_length\tr2_max_length"
    )?;
    for (locus, records) in grouped {
        let r1 = records
            .iter()
            .map(|record| (record.sample.clone(), record.r1.clone()))
            .collect::<Vec<_>>();
        let r2 = records
            .iter()
            .map(|record| (record.sample.clone(), record.r2.clone()))
            .collect::<Vec<_>>();
        let r1_max = r1
            .iter()
            .map(|(_, sequence)| sequence.len())
            .max()
            .unwrap_or(0);
        let r2_max = r2
            .iter()
            .map(|(_, sequence)| sequence.len())
            .max()
            .unwrap_or(0);
        write_fasta(&arms.join(format!("{locus}__R1.fasta")), &r1)?;
        write_fasta(&arms.join(format!("{locus}__R2.fasta")), &r2)?;
        writeln!(manifest, "{locus}\t{}\t{r1_max}\t{r2_max}", r1.len())?;
    }
    Ok(())
}

fn normalize_fasta_sequence(sequence: &str, path: &Path) -> io::Result<String> {
    let mut normalized = String::with_capacity(sequence.len());
    for base in sequence.bytes() {
        match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => normalized.push(base.to_ascii_uppercase() as char),
            b'U' => normalized.push('T'),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("FASTA contains a non-ACGTU base: {}", path.display()),
                ));
            }
        }
    }
    if normalized.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTA record has an empty sequence: {}", path.display()),
        ));
    }
    Ok(normalized)
}

fn read_fasta(path: &Path) -> io::Result<Vec<(String, String)>> {
    let mut result = Vec::new();
    let reader = BufReader::new(File::open(path)?);
    let mut name = None;
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = name.take() {
                result.push((previous, normalize_fasta_sequence(&sequence, path)?));
            }
            name = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or("unnamed")
                    .to_owned(),
            );
            sequence.clear();
        } else if !line.trim().is_empty() {
            if name.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("FASTA sequence occurs before a header: {}", path.display()),
                ));
            }
            sequence.push_str(line.trim());
        }
    }
    if let Some(previous) = name {
        result.push((previous, normalize_fasta_sequence(&sequence, path)?));
    }
    if result.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTA contains no records: {}", path.display()),
        ));
    }
    Ok(result)
}

fn validate_reference_loci(loci: &BTreeMap<String, ArmReferences>) -> io::Result<BTreeSet<String>> {
    let mut samples = BTreeSet::new();
    for (locus, (r1, r2)) in loci {
        if r1.is_empty() || r2.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference locus {locus} must contain both R1 and R2 arms"),
            ));
        }
        let r1_ids = r1.iter().map(|(id, _)| id).collect::<BTreeSet<_>>();
        let r2_ids = r2.iter().map(|(id, _)| id).collect::<BTreeSet<_>>();
        if r1_ids.len() != r1.len() || r2_ids.len() != r2.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference locus {locus} has duplicate sample identifiers"),
            ));
        }
        if r1_ids != r2_ids {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference locus {locus} has unmatched R1/R2 sample identifiers"),
            ));
        }
        samples.extend(r1_ids.into_iter().cloned());
    }
    Ok(samples)
}

fn longest_common_substring(left: &str, right: &str) -> usize {
    let (short, long) = if left.len() <= right.len() {
        (left.as_bytes(), right.as_bytes())
    } else {
        (right.as_bytes(), left.as_bytes())
    };
    let mut previous = vec![0usize; short.len() + 1];
    let mut best = 0;
    for &base in long {
        let mut current = vec![0usize; short.len() + 1];
        for (index, &other) in short.iter().enumerate() {
            if base == other {
                current[index + 1] = previous[index] + 1;
                best = best.max(current[index + 1]);
            }
        }
        previous = current;
    }
    best
}

fn recovered(
    path: &Path,
    references: &[(String, String)],
    breadth: f64,
) -> io::Result<Option<String>> {
    if !path.is_file() {
        return Ok(None);
    }
    let reference_length = references
        .iter()
        .map(|(_, sequence)| sequence.len())
        .max()
        .unwrap_or(0);
    if reference_length == 0 {
        return Ok(None);
    }
    let required = 15usize.max((reference_length as f64 * 0.30).ceil() as usize);
    let mut best = None;
    for (_, sequence) in read_fasta(path)? {
        if sequence.len() as f64 / (reference_length as f64) < breadth {
            continue;
        }
        let shared = references
            .iter()
            .map(|(_, reference)| longest_common_substring(&sequence, reference))
            .max()
            .unwrap_or(0);
        if shared >= required
            && best
                .as_ref()
                .is_none_or(|(best_shared, best_sequence): &(usize, String)| {
                    shared > *best_shared || (shared == *best_shared && sequence > *best_sequence)
                })
        {
            best = Some((shared, sequence));
        }
    }
    Ok(best.map(|(_, sequence)| sequence))
}

fn finalize(values: &HashMap<String, Vec<String>>) -> io::Result<()> {
    let reference_dir = required_path(values, "--reference");
    let recovery = required_path(values, "--recovery");
    let out = required_path(values, "--out");
    let samples = required_values(values, "--sample");
    let breadth = values
        .get("--min-arm-breadth")
        .and_then(|v| v.first())
        .map(String::as_str)
        .unwrap_or("0.80")
        .parse::<f64>()
        .ok()
        .filter(|v| (0.0..=1.0).contains(v))
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "--min-arm-breadth must be in [0,1]",
            )
        })?;
    let arms = reference_dir.join("arms");
    let mut loci = BTreeMap::<String, (Vec<(String, String)>, Vec<(String, String)>)>::new();
    for entry in fs::read_dir(&arms)? {
        let path = entry?.path();
        let Some(stem) = path.file_stem().and_then(|v| v.to_str()) else {
            continue;
        };
        let Some((locus, arm)) = stem.rsplit_once("__") else {
            continue;
        };
        let records = read_fasta(&path)?;
        match arm {
            "R1" => {
                loci.entry(locus.into())
                    .or_insert_with(|| (Vec::new(), Vec::new()))
                    .0 = records
            }
            "R2" => {
                loci.entry(locus.into())
                    .or_insert_with(|| (Vec::new(), Vec::new()))
                    .1 = records
            }
            _ => {}
        }
    }
    if loci.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reference has no RAD arm FASTA files",
        ));
    }
    let reference_samples = validate_reference_loci(&loci)?;
    let requested_samples = samples
        .iter()
        .map(|sample| safe_name(sample))
        .collect::<BTreeSet<_>>();
    if let Some(sample) = requested_samples.intersection(&reference_samples).next() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("recovery sample already exists in the RAD reference: {sample}"),
        ));
    }
    if out.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("matrix output already exists: {}", out.display()),
        ));
    }
    let phy_r1 = out.join("recovered_arms/R1");
    let phy_r2 = out.join("recovered_arms/R2");
    let strict_r1 = out.join("paired_arms/R1");
    let strict_r2 = out.join("paired_arms/R2");
    fs::create_dir_all(&phy_r1)?;
    fs::create_dir_all(&phy_r2)?;
    fs::create_dir_all(&strict_r1)?;
    fs::create_dir_all(&strict_r2)?;
    let mut report = File::create(out.join("rad_sample_locus.tsv"))?;
    writeln!(report, "sample\tlocus\tr1_status\tr2_status\tjoint_status")?;
    for (locus, (r1_ref, r2_ref)) in loci {
        let mut phy_left = r1_ref.clone();
        let mut phy_right = r2_ref.clone();
        let mut strict_left = r1_ref.clone();
        let mut strict_right = r2_ref.clone();
        let original = r1_ref
            .iter()
            .map(|(id, _)| id.clone())
            .collect::<BTreeSet<_>>();
        for sample in &samples {
            let base = recovery.join(sample).join("results");
            let left = recovered(&base.join(format!("{locus}__R1.fasta")), &r1_ref, breadth)?;
            let right = recovered(&base.join(format!("{locus}__R2.fasta")), &r2_ref, breadth)?;
            let left_state = if left.is_some() {
                "wgs_recovered"
            } else {
                "wgs_insufficient"
            };
            let right_state = if right.is_some() {
                "wgs_recovered"
            } else {
                "wgs_insufficient"
            };
            let joint = match (&left, &right) {
                (Some(_), Some(_)) => "rad_missing_wgs_recovered",
                (Some(_), None) | (None, Some(_)) => "partial_arm_recovery",
                _ => "unresolved",
            };
            writeln!(
                report,
                "{sample}\t{locus}\t{left_state}\t{right_state}\t{joint}"
            )?;
            if original.contains(sample) {
                continue;
            }
            if let Some(sequence) = left {
                phy_left.push((sample.clone(), sequence.clone()));
                if right.is_some() {
                    strict_left.push((sample.clone(), sequence));
                }
            }
            if let Some(sequence) = right {
                phy_right.push((sample.clone(), sequence.clone()));
                if left_state == "wgs_recovered" {
                    strict_right.push((sample.clone(), sequence));
                }
            }
        }
        write_fasta(&phy_r1.join(format!("{locus}.fasta")), &phy_left)?;
        write_fasta(&phy_r2.join(format!("{locus}.fasta")), &phy_right)?;
        write_fasta(&strict_r1.join(format!("{locus}.fasta")), &strict_left)?;
        write_fasta(&strict_r2.join(format!("{locus}.fasta")), &strict_right)?;
    }
    fs::write(out.join("README.txt"), "Unaligned RAD recovery matrices. R1 and R2 are independent observations; do not infer an intervening genomic sequence. paired_arms includes WGS samples only when both arms recover; recovered_arms retains supported individual arms. Use rad-validate before phylogenetic inference.\n")?;
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct LocalAlignment {
    score: i32,
    matches: usize,
    columns: usize,
    query_bases: usize,
    reference_bases: usize,
}

#[derive(Clone)]
struct BestAlignment {
    alignment: LocalAlignment,
    reference_id: String,
    reference_length: usize,
}

fn local_alignment(query: &str, reference: &str) -> LocalAlignment {
    // Smith-Waterman with deterministic tie breaking. RAD arms are short, so a
    // full local alignment is more reliable than seed-only comparison here.
    let query = query.as_bytes();
    let reference = reference.as_bytes();
    let width = reference.len() + 1;
    let mut score = vec![0i32; (query.len() + 1) * width];
    let mut direction = vec![0u8; score.len()]; // 1 diagonal; 2 up; 3 left
    let (mut best, mut best_i, mut best_j) = (0i32, 0usize, 0usize);
    for i in 1..=query.len() {
        for j in 1..=reference.len() {
            let index = i * width + j;
            let diagonal = score[(i - 1) * width + j - 1]
                + if query[i - 1] == reference[j - 1] {
                    2
                } else {
                    -3
                };
            let up = score[(i - 1) * width + j] - 4;
            let left = score[i * width + j - 1] - 4;
            let (value, step) = if diagonal > 0 && diagonal >= up && diagonal >= left {
                (diagonal, 1)
            } else if up > 0 && up >= left {
                (up, 2)
            } else if left > 0 {
                (left, 3)
            } else {
                (0, 0)
            };
            score[index] = value;
            direction[index] = step;
            if value > best {
                best = value;
                best_i = i;
                best_j = j;
            }
        }
    }
    let (mut i, mut j) = (best_i, best_j);
    let mut result = LocalAlignment {
        score: best,
        ..LocalAlignment::default()
    };
    while i > 0 && j > 0 {
        match direction[i * width + j] {
            0 => break,
            1 => {
                result.columns += 1;
                result.query_bases += 1;
                result.reference_bases += 1;
                if query[i - 1] == reference[j - 1] {
                    result.matches += 1;
                }
                i -= 1;
                j -= 1;
            }
            2 => {
                result.columns += 1;
                result.query_bases += 1;
                i -= 1;
            }
            3 => {
                result.columns += 1;
                result.reference_bases += 1;
                j -= 1;
            }
            _ => unreachable!(),
        }
    }
    result
}

fn best_alignment(query: &str, references: &[(String, String)]) -> Option<BestAlignment> {
    references
        .iter()
        .map(|(id, reference)| BestAlignment {
            alignment: local_alignment(query, reference),
            reference_id: id.clone(),
            reference_length: reference.len(),
        })
        .max_by(|left, right| {
            left.alignment
                .score
                .cmp(&right.alignment.score)
                .then_with(|| left.alignment.matches.cmp(&right.alignment.matches))
                .then_with(|| right.reference_id.cmp(&left.reference_id))
        })
}

fn decimal(values: &HashMap<String, Vec<String>>, flag: &str, default: f64) -> io::Result<f64> {
    let parsed = values
        .get(flag)
        .and_then(|items| items.first())
        .map(|item| item.parse::<f64>())
        .transpose()
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{flag} must be a number"),
            )
        })?
        .unwrap_or(default);
    if !(0.0..=1.0).contains(&parsed) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{flag} must be in [0, 1]"),
        ));
    }
    Ok(parsed)
}

fn validation_samples(
    recovery: &Path,
    values: &HashMap<String, Vec<String>>,
) -> io::Result<Vec<String>> {
    if let Some(samples) = values.get("--sample") {
        return Ok(samples.clone());
    }
    let mut samples = fs::read_dir(recovery)?
        .filter_map(Result::ok)
        .filter_map(|entry| {
            let path = entry.path();
            (path.is_dir() && path.join("results").is_dir())
                .then(|| entry.file_name().to_string_lossy().into_owned())
        })
        .collect::<Vec<_>>();
    samples.sort();
    if samples.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "recovery directory has no sample results/ directories",
        ));
    }
    Ok(samples)
}

fn validate(values: &HashMap<String, Vec<String>>) -> io::Result<()> {
    let reference_dir = required_path(values, "--reference");
    let recovery = required_path(values, "--recovery");
    let out = required_path(values, "--out");
    if out.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("validation output already exists: {}", out.display()),
        ));
    }
    let min_identity = decimal(values, "--min-identity", 0.90)?;
    let min_breadth = decimal(values, "--min-breadth", 0.80)?;
    let min_delta = decimal(values, "--min-delta", 0.05)?;
    let samples = validation_samples(&recovery, values)?;
    let arms = reference_dir.join("arms");
    let mut loci = BTreeMap::<String, (Vec<(String, String)>, Vec<(String, String)>)>::new();
    for entry in fs::read_dir(&arms)? {
        let path = entry?.path();
        let Some(stem) = path.file_stem().and_then(|value| value.to_str()) else {
            continue;
        };
        let Some((locus, arm)) = stem.rsplit_once("__") else {
            continue;
        };
        let records = read_fasta(&path)?;
        match arm {
            "R1" => {
                loci.entry(locus.into())
                    .or_insert_with(|| (Vec::new(), Vec::new()))
                    .0 = records
            }
            "R2" => {
                loci.entry(locus.into())
                    .or_insert_with(|| (Vec::new(), Vec::new()))
                    .1 = records
            }
            _ => {}
        }
    }
    if loci.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reference has no RAD arm FASTA files",
        ));
    }
    let reference_samples = validate_reference_loci(&loci)?;
    let requested_samples = samples
        .iter()
        .map(|sample| safe_name(sample))
        .collect::<BTreeSet<_>>();
    if let Some(sample) = requested_samples.intersection(&reference_samples).next() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("recovery sample already exists in the RAD reference: {sample}"),
        ));
    }
    for locus in loci.keys() {
        for sample in &samples {
            for arm in ["R1", "R2"] {
                let candidate = recovery
                    .join(sample)
                    .join("results")
                    .join(format!("{locus}__{arm}.fasta"));
                if candidate.is_file() {
                    read_fasta(&candidate)?;
                }
            }
        }
    }
    fs::create_dir_all(out.join("strict_arms/R1"))?;
    fs::create_dir_all(out.join("strict_arms/R2"))?;
    let mut report = File::create(out.join("rad_validation.tsv"))?;
    writeln!(report, "sample\tlocus\tarm\tstatus\tcandidate\tbest_reference\tidentity\tquery_breadth\treference_breadth\tbest_score\tforeign_score\tscore_delta")?;
    for (locus, (r1_refs, r2_refs)) in &loci {
        let mut strict_r1 = r1_refs.clone();
        let mut strict_r2 = r2_refs.clone();
        for sample in &samples {
            let mut accepted = [None, None];
            for (arm_index, (arm, own_refs)) in
                [("R1", r1_refs), ("R2", r2_refs)].into_iter().enumerate()
            {
                let candidate_path = recovery
                    .join(sample)
                    .join("results")
                    .join(format!("{locus}__{arm}.fasta"));
                let candidates = if candidate_path.is_file() {
                    read_fasta(&candidate_path)?
                } else {
                    Vec::new()
                };
                let mut choice: Option<(String, String, BestAlignment)> = None;
                for (candidate_id, candidate) in candidates {
                    let Some(best) = best_alignment(&candidate, own_refs) else {
                        continue;
                    };
                    if choice.as_ref().is_none_or(|(_, _, current)| {
                        best.alignment.score > current.alignment.score
                    }) {
                        choice = Some((candidate_id, candidate, best));
                    }
                }
                let (
                    status,
                    candidate_id,
                    best_ref,
                    identity,
                    query_breadth,
                    reference_breadth,
                    score,
                    foreign_score,
                    delta,
                    sequence,
                ) = if let Some((candidate_id, candidate, best)) = choice {
                    let foreign = loci
                        .iter()
                        .filter(|(other, _)| *other != locus)
                        .flat_map(|(_, refs)| {
                            if arm == "R1" {
                                refs.0.iter()
                            } else {
                                refs.1.iter()
                            }
                        })
                        .map(|(_, reference)| local_alignment(&candidate, reference).score)
                        .max()
                        .unwrap_or(0);
                    let identity = if best.alignment.columns == 0 {
                        0.0
                    } else {
                        best.alignment.matches as f64 / best.alignment.columns as f64
                    };
                    let query_breadth = if candidate.is_empty() {
                        0.0
                    } else {
                        best.alignment.query_bases as f64 / candidate.len() as f64
                    };
                    let reference_breadth = if best.reference_length == 0 {
                        0.0
                    } else {
                        best.alignment.reference_bases as f64 / best.reference_length as f64
                    };
                    let delta = if best.alignment.score <= 0 {
                        0.0
                    } else {
                        (best.alignment.score - foreign).max(0) as f64 / best.alignment.score as f64
                    };
                    let status = if query_breadth < min_breadth || reference_breadth < min_breadth {
                        "insufficient_coverage"
                    } else if identity < min_identity {
                        "low_identity"
                    } else if delta < min_delta {
                        "ambiguous_paralog"
                    } else {
                        "validated"
                    };
                    (
                        status,
                        candidate_id,
                        best.reference_id,
                        identity,
                        query_breadth,
                        reference_breadth,
                        best.alignment.score,
                        foreign,
                        delta,
                        candidate,
                    )
                } else {
                    (
                        "missing",
                        String::new(),
                        String::new(),
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0,
                        0.0,
                        String::new(),
                    )
                };
                writeln!(report, "{sample}\t{locus}\t{arm}\t{status}\t{candidate_id}\t{best_ref}\t{identity:.4}\t{query_breadth:.4}\t{reference_breadth:.4}\t{score}\t{foreign_score}\t{delta:.4}")?;
                if status == "validated" {
                    accepted[arm_index] = Some(sequence);
                }
            }
            if let [Some(left), Some(right)] = accepted {
                strict_r1.push((sample.clone(), left));
                strict_r2.push((sample.clone(), right));
            }
        }
        write_fasta(
            &out.join("strict_arms/R1").join(format!("{locus}.fasta")),
            &strict_r1,
        )?;
        write_fasta(
            &out.join("strict_arms/R2").join(format!("{locus}.fasta")),
            &strict_r2,
        )?;
    }
    fs::write(out.join("README.txt"), "Validated RAD arm matrix. A WGS sample enters strict only when both independent arms pass local alignment coverage, identity, and cross-locus score-separation checks. R1/R2 are not bridged.\n")?;
    Ok(())
}

const MAX_STACK_READS: usize = 128;

#[derive(Default)]
struct PairStack {
    r1: Vec<String>,
    r2: Vec<String>,
    total_pairs: usize,
    hashes: Vec<u64>,
}

fn stable_pair_hash(r1: &str, r2: &str) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for base in r1.bytes().chain([0xff]).chain(r2.bytes()) {
        hash ^= u64::from(base);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

impl PairStack {
    fn push(&mut self, r1: String, r2: String) {
        self.total_pairs += 1;
        let incoming = stable_pair_hash(&r1, &r2);
        if self.r1.len() < MAX_STACK_READS {
            self.r1.push(r1);
            self.r2.push(r2);
            self.hashes.push(incoming);
            return;
        }
        let (replace, maximum) = self
            .hashes
            .iter()
            .copied()
            .enumerate()
            .max_by_key(|(_, hash)| *hash)
            .expect("non-empty bounded stack");
        if incoming < maximum {
            self.r1[replace] = r1;
            self.r2[replace] = r2;
            self.hashes[replace] = incoming;
        }
    }
}

struct SampleConsensus {
    r1: String,
    r2: String,
    read_pairs: usize,
}

fn positive(
    values: &HashMap<String, Vec<String>>,
    flag: &str,
    default: usize,
) -> io::Result<usize> {
    let parsed = values
        .get(flag)
        .and_then(|items| items.first())
        .map(|item| item.parse::<usize>())
        .transpose()
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{flag} must be a positive integer"),
            )
        })?
        .unwrap_or(default);
    if parsed == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{flag} must be at least 1"),
        ));
    }
    Ok(parsed)
}

fn canonical_kmers(sequence: &[u8], k: usize) -> Vec<u64> {
    if !(1..=31).contains(&k) {
        return Vec::new();
    }
    let mask = (1u64 << (2 * k)) - 1;
    let high = 2 * (k - 1);
    let (mut forward, mut reverse, mut valid) = (0u64, 0u64, 0usize);
    let mut result = Vec::with_capacity(sequence.len().saturating_sub(k) + 1);
    for &base in sequence {
        let bits = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => {
                valid = 0;
                forward = 0;
                reverse = 0;
                continue;
            }
        };
        forward = ((forward << 2) | bits) & mask;
        reverse = (reverse >> 2) | ((3 - bits) << high);
        valid += 1;
        if valid >= k {
            result.push(forward.min(reverse));
        }
    }
    result
}

fn clean_rad_read(sequence: &[u8], overhang: &[u8], min_length: usize) -> Option<String> {
    let cleaned = sequence
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    if !overhang.is_empty() && !cleaned.starts_with(overhang) {
        return None;
    }
    let start = overhang.len();
    let trimmed = cleaned.get(start..)?;
    if trimmed.len() < min_length
        || !trimmed
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return None;
    }
    Some(String::from_utf8_lossy(trimmed).into_owned())
}

fn count_kmers(
    path: &Path,
    overhang: &[u8],
    k: usize,
    min_length: usize,
) -> io::Result<HashMap<u64, u32>> {
    let mut reader = FastxReader::open(path, FastxFormat::Fastq)?;
    let mut counts = HashMap::<u64, u32>::new();
    while let Some(record) = reader.next_record()? {
        if let Some(sequence) = clean_rad_read(&record.sequence, overhang, min_length) {
            for kmer in canonical_kmers(sequence.as_bytes(), k) {
                let count = counts.entry(kmer).or_default();
                *count = count.saturating_add(1);
            }
        }
    }
    Ok(counts)
}

fn solid_minimizer(
    sequence: &str,
    counts: &HashMap<u64, u32>,
    k: usize,
    min_count: u32,
) -> Option<u64> {
    canonical_kmers(sequence.as_bytes(), k)
        .into_iter()
        .filter(|kmer| counts.get(kmer).copied().unwrap_or(0) >= min_count)
        .min()
}

fn consensus(records: &[String], min_length: usize) -> Option<String> {
    if records.is_empty() {
        return None;
    }
    let length = records.iter().map(String::len).min()?;
    if length < min_length {
        return None;
    }
    let mut result = Vec::with_capacity(length);
    for index in 0..length {
        let mut count = [0usize; 4];
        for record in records {
            match record.as_bytes()[index] {
                b'A' => count[0] += 1,
                b'C' => count[1] += 1,
                b'G' => count[2] += 1,
                b'T' => count[3] += 1,
                _ => {}
            }
        }
        let (base, support) = count.iter().enumerate().max_by_key(|(_, value)| **value)?;
        // Retain ordinary heterozygosity, but reject stacks whose dominant path is not clear.
        if support * 2 < records.len() {
            return None;
        }
        result.push(*b"ACGT".get(base)?);
    }
    String::from_utf8(result).ok()
}

fn pair_id(header: &[u8]) -> &[u8] {
    let header = header.strip_prefix(b"@").unwrap_or(header);
    let end = header
        .iter()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(header.len());
    let id = &header[..end];
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
}

fn denovo(values: &HashMap<String, Vec<String>>) -> io::Result<()> {
    let out = required_path(values, "--out");
    let samples = required_values(values, "--sample");
    let r1_paths = required_values(values, "--read1");
    let r2_paths = required_values(values, "--read2");
    if samples.len() != r1_paths.len() || samples.len() != r2_paths.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "each --sample needs exactly one --read1 and --read2",
        ));
    }
    let k = positive(values, "--kmer", 31)?;
    if !(15..=31).contains(&k) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--kmer must be in [15, 31]",
        ));
    }
    let min_count = u32::try_from(positive(values, "--min-count", 3)?).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-count must fit in a 32-bit unsigned integer",
        )
    })?;
    let min_samples = positive(values, "--min-samples", 2)?;
    let min_length = positive(values, "--min-length", 60)?;
    let parse_overhang = |flag: &str| -> io::Result<Vec<u8>> {
        let overhang = values
            .get(flag)
            .and_then(|items| items.first())
            .map(|value| {
                value
                    .as_bytes()
                    .iter()
                    .map(|base| base.to_ascii_uppercase())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        if !overhang
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{flag} must contain only A/C/G/T"),
            ));
        }
        Ok(overhang)
    };
    let overhang_r1 = parse_overhang("--overhang")?;
    // R2 often starts in a size-selected insert rather than at the same cut site.
    // Only enforce a second end when the library design explicitly supplies one.
    let overhang_r2 = parse_overhang("--overhang-r2")?;
    let arms = out.join("arms");
    if out.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("reference output already exists: {}", out.display()),
        ));
    }
    let mut normalized_samples = BTreeSet::new();
    for sample in &samples {
        if !normalized_samples.insert(safe_name(sample)) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("duplicate sample after safe-name normalization: {sample}"),
            ));
        }
    }
    let mut loci = HashMap::<(u64, u64), BTreeMap<String, SampleConsensus>>::new();
    for ((sample, r1_path), r2_path) in samples.iter().zip(&r1_paths).zip(&r2_paths) {
        let counts1 = count_kmers(Path::new(r1_path), &overhang_r1, k, min_length)?;
        let counts2 = count_kmers(Path::new(r2_path), &overhang_r2, k, min_length)?;
        let mut sample_stacks = HashMap::<(u64, u64), PairStack>::new();
        let mut reader1 = FastxReader::open(Path::new(r1_path), FastxFormat::Fastq)?;
        let mut reader2 = FastxReader::open(Path::new(r2_path), FastxFormat::Fastq)?;
        loop {
            let left = reader1.next_record()?;
            let right = reader2.next_record()?;
            match (left, right) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("paired FASTQ length differs for sample {sample}"),
                    ))
                }
                (Some(left), Some(right)) => {
                    if pair_id(&left.header) != pair_id(&right.header) {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("paired FASTQ read identifiers differ for sample {sample}"),
                        ));
                    }
                    let Some(left) = clean_rad_read(&left.sequence, &overhang_r1, min_length)
                    else {
                        continue;
                    };
                    let Some(right) = clean_rad_read(&right.sequence, &overhang_r2, min_length)
                    else {
                        continue;
                    };
                    let (Some(seed1), Some(seed2)) = (
                        solid_minimizer(&left, &counts1, k, min_count),
                        solid_minimizer(&right, &counts2, k, min_count),
                    ) else {
                        continue;
                    };
                    sample_stacks
                        .entry((seed1, seed2))
                        .or_default()
                        .push(left, right);
                }
            }
        }
        for (seeds, stack) in sample_stacks {
            if stack.total_pairs < min_count as usize {
                continue;
            }
            let (Some(r1), Some(r2)) = (
                consensus(&stack.r1, min_length),
                consensus(&stack.r2, min_length),
            ) else {
                continue;
            };
            loci.entry(seeds).or_default().insert(
                safe_name(sample),
                SampleConsensus {
                    r1,
                    r2,
                    read_pairs: stack.total_pairs,
                },
            );
        }
    }
    fs::create_dir_all(&arms)?;
    let mut manifest = File::create(out.join("locus_manifest.tsv"))?;
    writeln!(
        manifest,
        "locus\tdenovo_status\tsamples\tr1_max_length\tr2_max_length\tr1_minimizer\tr2_minimizer"
    )?;
    let mut evidence = File::create(out.join("denovo_probe_evidence.tsv"))?;
    writeln!(
        evidence,
        "locus\tsample\tr1_reads\tr2_reads\tr1_status\tr2_status"
    )?;
    let mut ordinal = 0usize;
    let mut grouped = loci.into_iter().collect::<Vec<_>>();
    grouped.sort_by_key(|(seeds, _)| *seeds);
    for ((seed1, seed2), sample_stacks) in grouped {
        let usable = sample_stacks
            .into_iter()
            .map(|(sample, stack)| (sample, stack.read_pairs, stack.r1, stack.r2))
            .collect::<Vec<_>>();
        if usable.len() < min_samples {
            continue;
        }
        ordinal += 1;
        let locus = format!("denovo_{ordinal:06}");
        let r1 = usable
            .iter()
            .map(|(sample, _, seq, _)| (sample.clone(), seq.clone()))
            .collect::<Vec<_>>();
        let r2 = usable
            .iter()
            .map(|(sample, _, _, seq)| (sample.clone(), seq.clone()))
            .collect::<Vec<_>>();
        let r1_max = r1.iter().map(|(_, seq)| seq.len()).max().unwrap_or(0);
        let r2_max = r2.iter().map(|(_, seq)| seq.len()).max().unwrap_or(0);
        write_fasta(&arms.join(format!("{locus}__R1.fasta")), &r1)?;
        write_fasta(&arms.join(format!("{locus}__R2.fasta")), &r2)?;
        writeln!(
            manifest,
            "{locus}\tdenovo_candidate\t{}\t{r1_max}\t{r2_max}\t{seed1}\t{seed2}",
            usable.len()
        )?;
        for (sample, reads, _, _) in usable {
            writeln!(evidence, "{locus}\t{sample}\t{reads}\t{reads}\tpass\tpass")?;
        }
    }
    if ordinal == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no denovo RAD candidates passed the k-mer, depth, and sample-support filters",
        ));
    }
    Ok(())
}

fn main() {
    let args = env::args().skip(1).collect::<Vec<_>>();
    let Some(command) = args.first() else { usage() };
    let values = options(&args[1..]);
    let result = match command.as_str() {
        "reference" => reference(&values),
        "denovo" => denovo(&values),
        "validate" => validate(&values),
        "finalize" => finalize(&values),
        _ => usage(),
    };
    if let Err(error) = result {
        eprintln!("rad_workflow: {error}");
        std::process::exit(1);
    }
}
