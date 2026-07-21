use flate2::read::MultiGzDecoder;
use plotters::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

type Result<T> = std::result::Result<T, String>;
const EXCLUDED_FLAGS: u16 = 0x4 | 0x100 | 0x800;
const BASES: [u8; 6] = *b"-ACGNT";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    thresholds: Vec<f64>,
    out: PathBuf,
    prefix: Option<String>,
    min_depth: u32,
    fill: u8,
    max_del: usize,
    width: usize,
    save_mutations: bool,
}
impl Default for Args {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            thresholds: vec![0.25],
            out: PathBuf::new(),
            prefix: None,
            min_depth: 2,
            fill: b'-',
            max_del: 150,
            width: 0,
            save_mutations: true,
        }
    }
}

#[derive(Clone)]
struct Reference {
    name: String,
    counts: Vec<[u32; 6]>,
    coverage: Vec<u32>,
    insertions: BTreeMap<isize, Vec<String>>,
}

fn take(argv: &[String], i: &mut usize, flag: &str) -> Result<String> {
    *i += 1;
    argv.get(*i)
        .cloned()
        .ok_or_else(|| format!("{flag} requires a value"))
}
fn parse_args() -> Result<Args> {
    let argv: Vec<String> = env::args().skip(1).collect();
    let mut a = Args::default();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-i" | "--input" => a.input = take(&argv, &mut i, "-i")?.into(),
            "-c" | "--consensus-thresholds" => {
                a.thresholds = take(&argv, &mut i, "-c")?
                    .split(',')
                    .map(|x| {
                        x.parse::<f64>()
                            .map_err(|_| "invalid consensus threshold".to_string())
                    })
                    .collect::<Result<_>>()?;
            }
            "-o" | "--outfolder" => a.out = take(&argv, &mut i, "-o")?.into(),
            "-p" | "--prefix" => a.prefix = Some(take(&argv, &mut i, "-p")?),
            "-m" | "--min-depth" => {
                a.min_depth = take(&argv, &mut i, "-m")?
                    .parse()
                    .map_err(|_| "invalid minimum depth".to_string())?
            }
            "-f" | "--fill" => {
                let s = take(&argv, &mut i, "-f")?;
                if s.len() != 1 {
                    return Err("--fill must be one byte".into());
                }
                a.fill = s.as_bytes()[0];
            }
            "-d" | "--maxdel" => {
                a.max_del = take(&argv, &mut i, "-d")?
                    .parse()
                    .map_err(|_| "invalid maximum deletion".to_string())?
            }
            "-n" => {
                a.width = take(&argv, &mut i, "-n")?
                    .parse()
                    .map_err(|_| "invalid line width".to_string())?
            }
            "-s" => a.save_mutations = take(&argv, &mut i, "-s")? != "0",
            "-h" | "--help" => {
                println!("Usage: build_consensus -i INPUT.sam[.gz] -c THRESHOLD[,THRESHOLD] -o DIR [options]\n\nBuild read-supported IUPAC consensus FASTA from SAM alignments.");
                std::process::exit(0)
            }
            x => return Err(format!("unknown option {x}")),
        };
        i += 1;
    }
    if a.input.as_os_str().is_empty() {
        return Err("-i/--input is required".into());
    }
    if a.out.as_os_str().is_empty() {
        a.out = a.input.parent().unwrap_or(Path::new(".")).to_path_buf();
    }
    if a.thresholds.is_empty() || a.thresholds.iter().any(|x| !(0.0 < *x && *x <= 1.0)) {
        return Err("consensus thresholds must be in (0, 1]".into());
    }
    Ok(a)
}
fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).map_err(|e| e.to_string())?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}
fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'-' => Some(0),
        b'A' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        b'N' => Some(4),
        b'T' | b'U' => Some(5),
        _ => None,
    }
}
fn parse_headers(input: &mut dyn BufRead) -> Result<(Vec<Reference>, Vec<String>)> {
    let mut refs = Vec::new();
    let mut rows = Vec::new();
    let mut line = String::new();
    while input.read_line(&mut line).map_err(|e| e.to_string())? != 0 {
        let row = line.trim_end_matches(['\r', '\n']).to_string();
        line.clear();
        if row.starts_with("@SQ\t") {
            let mut name = None;
            let mut length = None;
            for field in row.split('\t').skip(1) {
                if let Some(x) = field.strip_prefix("SN:") {
                    name = Some(x.to_string())
                }
                if let Some(x) = field.strip_prefix("LN:") {
                    length = Some(
                        x.parse::<usize>()
                            .map_err(|_| "invalid @SQ length".to_string())?,
                    )
                }
            }
            let name = name.ok_or_else(|| "@SQ is missing SN".to_string())?;
            let length = length.ok_or_else(|| "@SQ is missing LN".to_string())?;
            refs.push(Reference {
                name,
                counts: vec![[0; 6]; length],
                coverage: vec![0; length],
                insertions: BTreeMap::new(),
            });
        } else if !row.starts_with('@') && !row.is_empty() {
            rows.push(row);
            break;
        }
    }
    Ok((refs, rows))
}
fn cigar_ops(cigar: &str) -> Result<Vec<(usize, u8)>> {
    let mut out = Vec::new();
    let mut n = 0usize;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n
                .checked_mul(10)
                .and_then(|x| x.checked_add((b - b'0') as usize))
                .ok_or_else(|| "CIGAR length overflow".to_string())?;
        } else {
            if n == 0
                || !matches!(
                    b,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'='
                )
            {
                return Err(format!("malformed CIGAR {cigar:?}"));
            };
            out.push((n, b));
            n = 0;
        }
    }
    if n != 0 || out.is_empty() {
        return Err(format!("malformed CIGAR {cigar:?}"));
    };
    Ok(out)
}
fn add_alignment(
    reference: &mut Reference,
    cigar: &str,
    sequence: &[u8],
    position: usize,
    max_del: usize,
) -> Result<()> {
    let ops = cigar_ops(cigar)?;
    let mut q = 0usize;
    let mut r = position;
    let mut max_gap = 0usize;
    let mut events: Vec<(usize, Option<u8>)> = Vec::new();
    let mut insertions = Vec::new();
    for (n, op) in ops {
        match op {
            b'M' | b'X' | b'=' => {
                if q + n > sequence.len() {
                    return Err("CIGAR consumes beyond SEQ".into());
                };
                for &b in &sequence[q..q + n] {
                    events.push((r, Some(b)));
                    r += 1;
                }
                q += n;
            }
            b'D' | b'N' => {
                max_gap = max_gap.max(n);
                for _ in 0..n {
                    events.push((r, None));
                    r += 1;
                }
            }
            b'I' => {
                if q + n > sequence.len() {
                    return Err("CIGAR consumes beyond SEQ".into());
                };
                insertions.push((
                    r as isize - 1,
                    String::from_utf8_lossy(&sequence[q..q + n]).to_ascii_uppercase(),
                ));
                q += n;
            }
            b'S' => {
                if q + n > sequence.len() {
                    return Err("CIGAR consumes beyond SEQ".into());
                };
                q += n;
            }
            b'H' | b'P' => {}
            _ => unreachable!(),
        }
    }
    if q != sequence.len() {
        return Err(format!(
            "CIGAR consumes {q} query bases, but SEQ has length {}",
            sequence.len()
        ));
    }
    for (pos, base) in events {
        if pos >= reference.counts.len() {
            return Err("alignment exceeds @SQ reference length".into());
        };
        if max_gap > max_del && base.is_none() {
            continue;
        };
        let b = base.unwrap_or(b'-');
        let idx = base_index(b).ok_or_else(|| format!("unexpected base {}", b as char))?;
        reference.counts[pos][idx] += 1;
    }
    for (pos, motif) in insertions {
        reference.insertions.entry(pos).or_default().push(motif);
    }
    Ok(())
}
fn parse_sam(path: &Path, max_del: usize) -> Result<(Vec<Reference>, u64, u64)> {
    let mut input = reader(path)?;
    let (mut refs, mut first) = parse_headers(&mut *input)?;
    let names: HashMap<String, usize> = refs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.clone(), i))
        .collect();
    let mut total = 0;
    let mut mapped = 0;
    loop {
        let row = if let Some(x) = first.pop() {
            x
        } else {
            let mut line = String::new();
            if input.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
                break;
            };
            line.trim_end_matches(['\r', '\n']).to_string()
        };
        if row.is_empty() || row.starts_with('@') {
            continue;
        };
        total += 1;
        let f: Vec<_> = row.split('\t').collect();
        if f.len() < 11 {
            return Err("malformed SAM record".into());
        };
        let flag: u16 = f[1].parse().map_err(|_| "invalid SAM flag".to_string())?;
        if flag & EXCLUDED_FLAGS != 0 || f[5] == "*" {
            continue;
        };
        let idx = *names
            .get(f[2])
            .ok_or_else(|| format!("SAM reference {} absent from @SQ", f[2]))?;
        let pos: usize = f[3]
            .parse::<usize>()
            .map_err(|_| "invalid SAM position".to_string())?
            .checked_sub(1)
            .ok_or_else(|| "SAM positions are 1-based".to_string())?;
        add_alignment(&mut refs[idx], f[5], f[9].as_bytes(), pos, max_del)?;
        mapped += 1;
    }
    Ok((refs, total, mapped))
}
fn code(chars: &[u8]) -> u8 {
    // Exact compatibility table for the legacy Python `amb` mapping.
    match chars {
        b"-" => b'-',
        b"A" => b'A',
        b"C" => b'C',
        b"G" => b'G',
        b"N" => b'N',
        b"T" => b'T',
        b"-A" => b'a',
        b"-C" => b'c',
        b"-G" => b'g',
        b"-N" => b'n',
        b"-T" => b't',
        b"AC" => b'M',
        b"AG" => b'R',
        b"AN" => b'a',
        b"AT" => b'W',
        b"CG" => b'S',
        b"CN" => b'c',
        b"CT" => b'Y',
        b"GN" => b'g',
        b"GT" => b'K',
        b"NT" => b't',
        b"-AC" => b'm',
        b"-AG" => b'r',
        b"-AN" => b'a',
        b"-AT" => b'w',
        b"-CG" => b's',
        b"-CN" => b'c',
        b"-CT" => b'y',
        b"-GN" => b'g',
        b"-GT" => b'k',
        b"-NT" => b't',
        b"ACG" => b'V',
        b"ACN" => b'm',
        b"ACT" => b'H',
        b"AGN" => b'r',
        b"AGT" => b'D',
        b"ANT" => b'w',
        b"CGN" => b's',
        b"CGT" => b'B',
        b"CNT" => b'y',
        b"GNT" => b'k',
        b"-ACG" => b'v',
        b"-ACN" => b'm',
        b"-ACT" => b'h',
        b"-AGN" => b'r',
        b"-AGT" => b'd',
        b"-ANT" => b'w',
        b"-CGN" => b's',
        b"-CGT" => b'b',
        b"-CNT" => b'y',
        b"-GNT" => b'k',
        b"ACGN" => b'v',
        b"ACGT" => b'N',
        b"ACNT" => b'h',
        b"AGNT" => b'd',
        b"CGNT" => b'b',
        b"-ACGN" => b'v',
        b"-ACGT" => b'N',
        b"-ACNT" => b'h',
        b"-AGNT" => b'd',
        b"-CGNT" => b'b',
        b"-ACGNT" => b'N',
        _ => b'N',
    }
}

fn called(counts: &[u32; 6], coverage: u32, threshold: f64) -> u8 {
    let mut groups: BTreeMap<u32, Vec<u8>> = BTreeMap::new();
    for (i, &n) in counts.iter().enumerate() {
        groups.entry(n).or_default().push(BASES[i]);
    }
    let mut chosen = Vec::new();
    let mut sum = 0u32;
    for (n, bases) in groups.into_iter().rev() {
        if (sum as f64) >= threshold * coverage as f64 {
            break;
        };
        sum += n.saturating_mul(bases.len() as u32);
        chosen.extend(bases);
    }
    chosen.sort_unstable();
    code(&chosen)
}
fn insertion_columns(reference: &Reference) -> BTreeMap<isize, Vec<[u32; 6]>> {
    let mut out = BTreeMap::new();
    for (&pos, motifs) in &reference.insertions {
        let max = motifs.iter().map(|x| x.len()).max().unwrap_or(0);
        let mut cols = vec![[0; 6]; max];
        for motif in motifs {
            for (i, b) in motif.bytes().enumerate() {
                if let Some(x) = base_index(b) {
                    cols[i][x] += 1;
                }
            }
        }
        let coverage = reference.coverage
            [pos.clamp(0, reference.coverage.len().saturating_sub(1) as isize) as usize];
        for col in &mut cols {
            let observed: u32 = col.iter().sum();
            col[0] = coverage.saturating_sub(observed);
        }
        out.insert(pos, cols);
    }
    out
}
fn build_sequence(reference: &Reference, threshold: f64, min_depth: u32, fill: u8) -> Vec<u8> {
    let ins = insertion_columns(reference);
    let mut out = Vec::new();
    let append = |out: &mut Vec<u8>, pos: isize| {
        if let Some(cols) = ins.get(&pos) {
            let cov = reference.coverage
                [pos.clamp(0, reference.coverage.len().saturating_sub(1) as isize) as usize];
            if cov >= min_depth {
                for col in cols {
                    let b = called(col, cov, threshold);
                    if b != b'-' {
                        out.push(b)
                    }
                }
            }
        }
    };
    append(&mut out, -1);
    for (pos, (&cov, counts)) in reference.coverage.iter().zip(&reference.counts).enumerate() {
        if cov >= min_depth {
            out.push(called(counts, cov, threshold));
            append(&mut out, pos as isize)
        } else {
            out.push(fill)
        }
    }
    out
}
fn python_round2(value: f64) -> String {
    let mut text = format!("{value:.2}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.push('0');
    }
    text
}
fn prefix(path: &Path) -> String {
    // Match Python os.path.splitext(basename): input.sam.gz -> input.sam.
    let file = path.file_name().unwrap_or_default().to_string_lossy();
    file.rsplit_once('.')
        .map(|x| x.0)
        .unwrap_or(&file)
        .to_string()
}
fn write_fasta(args: &Args, refs: &mut [Reference]) -> Result<()> {
    fs::create_dir_all(&args.out).map_err(|e| e.to_string())?;
    let name = args.prefix.clone().unwrap_or_else(|| prefix(&args.input));
    let mut out =
        File::create(args.out.join(format!("{name}.fasta"))).map_err(|e| e.to_string())?;
    for reference in refs.iter_mut() {
        for cov in &mut reference.coverage {
            *cov = 0
        }
        for (pos, counts) in reference.counts.iter().enumerate() {
            reference.coverage[pos] = counts.iter().sum()
        }
        if reference.coverage.iter().all(|&x| x == 0) {
            continue;
        }
        for &t in &args.thresholds {
            let seq = build_sequence(reference, t, args.min_depth, args.fill);
            if seq.iter().all(|&x| x == b'-') {
                continue;
            }
            let sum: u64 = reference.coverage.iter().map(|&x| x as u64).sum();
            let cov = sum as f64 / seq.len().max(1) as f64;
            let length = seq.iter().filter(|&&x| x != b'-').count();
            writeln!(
                out,
                ">{name}|c{} reference:{} coverage:{} length:{} consensus_threshold:{}%",
                (t * 100.0) as u32,
                reference.name,
                python_round2(cov),
                length,
                (t * 100.0) as u32
            )
            .map_err(|e| e.to_string())?;
            if args.width == 0 {
                out.write_all(&seq).map_err(|e| e.to_string())?;
                writeln!(out).map_err(|e| e.to_string())?
            } else {
                for line in seq.chunks(args.width) {
                    out.write_all(line).map_err(|e| e.to_string())?;
                    writeln!(out).map_err(|e| e.to_string())?
                }
            }
        }
    }
    Ok(())
}
fn write_mutation_plot(args: &Args, refs: &[Reference]) -> Result<()> {
    let mut values = Vec::new();
    for reference in refs {
        for (counts, &coverage) in reference.counts.iter().zip(&reference.coverage) {
            if coverage < args.min_depth {
                continue;
            }
            let normal = [counts[1], counts[2], counts[3], counts[5]];
            let total: u32 = normal.iter().sum();
            if total == 0 {
                continue;
            }
            let proportions: Vec<f64> = normal.iter().map(|&x| x as f64 / total as f64).collect();
            if proportions.iter().copied().fold(0.0, f64::max) < 0.9 {
                values.extend(proportions.into_iter().filter(|&x| x > 0.1 && x < 0.9));
            }
        }
    }
    if values.len() <= 4 {
        return Ok(());
    }
    let mut low = values.iter().copied().fold(f64::INFINITY, f64::min);
    let mut high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if (high - low).abs() < f64::EPSILON {
        low = (low - 0.05).max(0.0);
        high = (high + 0.05).min(1.0);
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / values.len() as f64;
    let bandwidth = (variance.sqrt() * 0.2).max(0.01);
    let points: Vec<(f64, f64)> = (0..200)
        .map(|i| {
            let x = low + (high - low) * i as f64 / 199.0;
            let density = values
                .iter()
                .map(|v| (-0.5 * ((x - v) / bandwidth).powi(2)).exp())
                .sum::<f64>()
                / (values.len() as f64 * bandwidth * (2.0 * std::f64::consts::PI).sqrt());
            (x, density)
        })
        .collect();
    let top = points.iter().map(|x| x.1).fold(0.0, f64::max).max(1e-9);
    let prefix = args
        .prefix
        .clone()
        .unwrap_or_else(|| prefix(&args.input))
        .replace("_tmp", "");
    let output = args.out.join(format!("{prefix}.png"));
    // Primitive bitmap drawing avoids fontconfig/FreeType dependencies.
    let area = BitMapBackend::new(&output, (800, 500)).into_drawing_area();
    area.fill(&WHITE).map_err(|e| e.to_string())?;
    let (left, right, top_y, bottom) = (40_i32, 780_i32, 20_i32, 470_i32);
    area.draw(&PathElement::new(
        vec![(left, bottom), (right, bottom)],
        BLACK,
    ))
    .map_err(|e| e.to_string())?;
    area.draw(&PathElement::new(
        vec![(left, bottom), (left, top_y)],
        BLACK,
    ))
    .map_err(|e| e.to_string())?;
    let curve: Vec<(i32, i32)> = points
        .iter()
        .map(|(x, y)| {
            let px = left as f64 + (x - low) / (high - low) * (right - left) as f64;
            let py = bottom as f64 - y / (top * 1.05) * (bottom - top_y) as f64;
            (px.round() as i32, py.round() as i32)
        })
        .collect();
    area.draw(&PathElement::new(curve, BLUE))
        .map_err(|e| e.to_string())?;
    area.present().map_err(|e| e.to_string())?;
    Ok(())
}
fn main() {
    match parse_args().and_then(|args| {
        let (input, total, mapped) = parse_sam(&args.input, args.max_del)?;
        println!(
            "A total of {total} reads were processed, out of which, {mapped} reads were mapped."
        );
        let mut refs = input;
        write_fasta(&args, &mut refs)?;
        if args.save_mutations {
            write_mutation_plot(&args, &refs)?;
        }
        Ok(())
    }) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cigar_padding_and_insertions() {
        let ops = cigar_ops("2M1I1P2M").unwrap();
        assert_eq!(ops.len(), 4);
        let mut r = Reference {
            name: "r".into(),
            counts: vec![[0; 6]; 4],
            coverage: vec![0; 4],
            insertions: BTreeMap::new(),
        };
        add_alignment(&mut r, "2M1I2M", b"ACTGT", 0, 150).unwrap();
        assert_eq!(r.insertions[&1], vec!["T"]);
        assert_eq!(r.counts[3][5], 1);
    }
    #[test]
    fn iupac_ties() {
        assert_eq!(code(b"AC"), b'M');
        assert_eq!(called(&[0, 1, 1, 0, 0, 0], 2, 0.75), b'M');
        assert_eq!(called(&[0, 1, 0, 0, 0, 0], 1, 0.75), b'A');
        assert_eq!(code(b"AN"), b'a');
        assert_eq!(code(b"ACN"), b'm');
        assert_eq!(code(b"-ACN"), b'm');
    }
}
