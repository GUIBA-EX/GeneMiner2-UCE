use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::env;
use std::ffi::OsStr;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::{Arc, Mutex};
use std::thread;

const RUN_LEN_CONST: f64 = 0.577_215_664_9_f64 / std::f64::consts::LN_2 - 1.5;
const THR_P95_2T: f64 = 1.96;
const THR_1E5_1T: f64 = 3.74;
const TOLERANCE: f64 = 1e-5;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FileType {
    Fasta,
    Fastq,
}

impl FileType {
    fn output_ext(self) -> &'static str {
        match self {
            FileType::Fasta => ".fasta",
            FileType::Fastq => ".fq",
        }
    }
}

#[derive(Clone)]
struct Record {
    title: String,
    seq: String,
    qual: String,
}

#[derive(Clone)]
struct Task {
    name: String,
    out_dir: PathBuf,
    ref_path: PathBuf,
    read_paths: Vec<PathBuf>,
    log_path: Option<PathBuf>,
    min_depth: i64,
    max_depth: i64,
    max_size: i64,
    copy_only: bool,
    keep_temporaries: bool,
    keep_linked_mates: bool,
    kmer_size: usize,
}

#[derive(Default)]
struct Args {
    se_dir: Option<PathBuf>,
    pe_dir: Option<PathBuf>,
    ref_dir: PathBuf,
    out_dir: PathBuf,
    log_file: Option<PathBuf>,
    min_depth: i64,
    max_depth: i64,
    max_size: i64,
    copy_only: bool,
    keep_temporaries: bool,
    keep_linked_mates: bool,
    use_gm2_format: bool,
    kmer_size: usize,
    processes: usize,
}

struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    fn next_record(&mut self) -> io::Result<Option<Record>> {
        let mut title = String::new();

        if self.reader.read_line(&mut title)? == 0 {
            return Ok(None);
        }

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();

        if self.reader.read_line(&mut seq)? == 0
            || self.reader.read_line(&mut plus)? == 0
            || self.reader.read_line(&mut qual)? == 0
        {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "truncated FASTQ record"));
        }

        let title = title.trim_end().strip_prefix('@').unwrap_or(title.trim_end()).to_string();
        Ok(Some(Record {
            title,
            seq: seq.trim_end().to_string(),
            qual: qual.trim_end().to_string(),
        }))
    }
}

struct FastaReader<R: BufRead> {
    reader: R,
    pending_title: Option<String>,
    done: bool,
}

impl<R: BufRead> FastaReader<R> {
    fn next_record(&mut self) -> io::Result<Option<Record>> {
        if self.done {
            return Ok(None);
        }

        let title = match self.pending_title.take() {
            Some(title) => title,
            None => {
                let mut line = String::new();
                loop {
                    line.clear();
                    if self.reader.read_line(&mut line)? == 0 {
                        self.done = true;
                        return Ok(None);
                    }

                    if let Some(title) = line.trim_end().strip_prefix('>') {
                        break title.to_string();
                    }
                }
            }
        };

        let mut seq = String::new();
        let mut line = String::new();

        loop {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                self.done = true;
                break;
            }

            if let Some(next_title) = line.trim_end().strip_prefix('>') {
                self.pending_title = Some(next_title.to_string());
                break;
            }

            seq.push_str(line.trim());
        }

        Ok(Some(Record {
            title,
            seq,
            qual: String::new(),
        }))
    }
}

struct Gm2Reader<R: Read> {
    reader: R,
    read_id: usize,
    suffix: String,
}

impl<R: Read> Gm2Reader<R> {
    fn next_record(&mut self) -> io::Result<Option<Record>> {
        let mut header = [0_u8; 6];
        let mut read_bytes = 0;

        while read_bytes < header.len() {
            let n = self.reader.read(&mut header[read_bytes..])?;
            if n == 0 {
                if read_bytes == 0 {
                    return Ok(None);
                }

                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated GM2 header"));
            }

            read_bytes += n;
        }

        let rec_len = ((header[0] as usize) << 16) | header[1] as usize;
        let has_phr = (header[2] & 0x80) != 0;
        let seq_len = (((header[2] & 0x7f) as usize) << 16) | header[3] as usize;

        if rec_len == 0 {
            return self.next_record();
        }

        let mut record = vec![0_u8; rec_len];
        self.reader.read_exact(&mut record)?;

        if seq_len == 0 {
            return self.next_record();
        }

        let (seq, qual) = parse_gm2_record(&record, has_phr, seq_len)?;
        self.read_id += 1;

        Ok(Some(Record {
            title: format!("read_{}{}", self.read_id, self.suffix),
            seq,
            qual,
        }))
    }
}

enum RecordReader {
    Fasta(FastaReader<BufReader<File>>),
    Fastq(FastqReader<BufReader<File>>),
    Gm2(Gm2Reader<BufReader<File>>),
}

impl RecordReader {
    fn from_path(path: &Path, file_type: FileType, gm2_format: bool, suffix: String) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        if gm2_format {
            return Ok(RecordReader::Gm2(Gm2Reader {
                reader,
                read_id: 0,
                suffix,
            }));
        }

        match file_type {
            FileType::Fasta => Ok(RecordReader::Fasta(FastaReader {
                reader,
                pending_title: None,
                done: false,
            })),
            FileType::Fastq => Ok(RecordReader::Fastq(FastqReader { reader })),
        }
    }

    fn next_record(&mut self) -> io::Result<Option<Record>> {
        match self {
            RecordReader::Fasta(reader) => reader.next_record(),
            RecordReader::Fastq(reader) => reader.next_record(),
            RecordReader::Gm2(reader) => reader.next_record(),
        }
    }
}

fn parse_gm2_record(record: &[u8], has_phr: bool, seq_len: usize) -> io::Result<(String, String)> {
    let mut i = 0_usize;
    let mut j = 0_usize;
    let mut last_chunk = 0_u8;
    let mut seq = vec![0_u8; seq_len];

    while i < seq_len {
        if j >= record.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated GM2 sequence"));
        }

        let delta = record[j] & 0x7f;
        let chunk = delta ^ last_chunk;
        let rep_num = ((chunk >> 5) & 3) as usize + 1;

        for k in 0..rep_num {
            if i + k < seq_len {
                seq[i + k] = (chunk & 31) + 64;
            }
        }

        i += rep_num;
        last_chunk = chunk;

        if record[j] >> 7 != 0 {
            let chunk = delta ^ last_chunk;
            let rep_num = ((chunk >> 5) & 3) as usize + 1;

            for k in 0..rep_num {
                if i + k < seq_len {
                    seq[i + k] = (chunk & 31) + 64;
                }
            }

            i += rep_num;
            last_chunk = chunk;
        }

        j += 1;
    }

    if !has_phr {
        return Ok((String::from_utf8_lossy(&seq).into_owned(), String::new()));
    }

    i = 0;
    let mut qual = vec![0_u8; seq_len];

    while i < seq_len {
        if j >= record.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated GM2 quality"));
        }

        let chunk = record[j];

        if chunk & 0x80 != 0 {
            if j + 1 >= record.len() {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated GM2 quality run"));
            }

            let chunk2 = record[j + 1];
            let rep_num = ((chunk & 0x7f) | (chunk2 & 0x80)) as usize + 1;

            for k in 0..rep_num {
                if i + k < seq_len {
                    qual[i + k] = chunk2 & 0x7f;
                }
            }

            i += rep_num;
            j += 2;
        } else {
            qual[i] = chunk;
            i += 1;
            j += 1;
        }
    }

    Ok((
        String::from_utf8_lossy(&seq).into_owned(),
        String::from_utf8_lossy(&qual).into_owned(),
    ))
}

fn main() {
    if let Err(error) = real_main() {
        eprintln!("Error: {error}");
        process::exit(1);
    }
}

fn real_main() -> Result<(), String> {
    let args = parse_args(env::args().skip(1).collect())?;
    let read_dict = get_read_dict(args.se_dir.as_deref(), args.pe_dir.as_deref(), args.use_gm2_format)?;
    let ref_dict = get_ref_dict(&args.ref_dir)?;
    fs::create_dir_all(args.out_dir.join("large_files")).map_err(|e| e.to_string())?;
    run(args, read_dict, ref_dict)
}

fn parse_args(argv: Vec<String>) -> Result<Args, String> {
    if argv.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        process::exit(0);
    }

    let mut args = Args {
        min_depth: 50,
        max_depth: 768,
        max_size: 6,
        kmer_size: 31,
        processes: 1,
        ..Args::default()
    };

    let mut i = 0_usize;

    while i < argv.len() {
        let key = &argv[i];

        match key.as_str() {
            "-qs" | "--se-dir" => args.se_dir = Some(next_path(&argv, &mut i, key)?),
            "-qd" | "--pe-dir" => args.pe_dir = Some(next_path(&argv, &mut i, key)?),
            "-r" | "--ref-dir" => args.ref_dir = next_path(&argv, &mut i, key)?,
            "-o" | "--out-dir" => args.out_dir = next_path(&argv, &mut i, key)?,
            "--log-file" => args.log_file = Some(next_path(&argv, &mut i, key)?),
            "--min-depth" => args.min_depth = next_parse(&argv, &mut i, key)?,
            "--max-depth" => args.max_depth = next_parse(&argv, &mut i, key)?,
            "--max-size" => args.max_size = next_parse(&argv, &mut i, key)?,
            "--copy-only" => args.copy_only = true,
            "--keep-temporaries" => args.keep_temporaries = true,
            "--keep-linked-mates" => args.keep_linked_mates = true,
            "--use-gm2-format" => args.use_gm2_format = true,
            "-kf" | "--kmer-size" => args.kmer_size = next_parse(&argv, &mut i, key)?,
            "-p" | "--processes" => args.processes = next_parse(&argv, &mut i, key)?,
            _ => return Err(format!("unknown argument {key}")),
        }

        i += 1;
    }

    if args.se_dir.is_some() == args.pe_dir.is_some() {
        return Err("exactly one of --se-dir or --pe-dir is required".to_string());
    }

    if args.ref_dir.as_os_str().is_empty() {
        return Err("--ref-dir is required".to_string());
    }

    if args.out_dir.as_os_str().is_empty() {
        return Err("--out-dir is required".to_string());
    }

    if args.kmer_size == 0 {
        return Err("--kmer-size must be positive".to_string());
    }

    if args.min_depth < 0 {
        return Err("--min-depth must be zero or positive".to_string());
    }

    if args.max_depth <= 0 {
        return Err("--max-depth must be positive".to_string());
    }

    if args.max_size <= 0 {
        return Err("--max-size must be positive".to_string());
    }

    if args.processes == 0 {
        return Err("--processes must be positive".to_string());
    }

    Ok(args)
}

fn next_path(argv: &[String], i: &mut usize, key: &str) -> Result<PathBuf, String> {
    *i += 1;
    argv.get(*i)
        .map(PathBuf::from)
        .ok_or_else(|| format!("{key} requires a value"))
}

fn next_parse<T: std::str::FromStr>(argv: &[String], i: &mut usize, key: &str) -> Result<T, String> {
    *i += 1;
    argv.get(*i)
        .ok_or_else(|| format!("{key} requires a value"))?
        .parse()
        .map_err(|_| format!("{key} has an invalid value"))
}

fn print_help() {
    println!(
        "usage: main_refilter_new -r REF_DIR -o OUT_DIR (-qs SE_DIR | -qd PE_DIR) [options]\n\
\n\
An improved NGS filtering gadget based on k-mers.\n\
\n\
options:\n\
  -qs, --se-dir DIR              Directory with single-read sequencing data\n\
  -qd, --pe-dir DIR              Directory with paired-end sequencing data\n\
  -r, --ref-dir DIR              Directory with reference sequences\n\
  -o, --out-dir DIR              Output directory\n\
  --log-file FILE                Log file\n\
  --min-depth INT                Min allowed coverage (default: 50)\n\
  --max-depth INT                Max allowed coverage (default: 768)\n\
  --max-size INT                 Max allowed size in million bases (default: 6)\n\
  --copy-only                    Interleave paired reads without filtering\n\
  --keep-temporaries             Keep temporary files\n\
  --keep-linked-mates            Keep both paired mates when either mate passes filtering\n\
  --use-gm2-format               Read reads from compressed binary format\n\
  -kf, --kmer-size INT           K-mer size (default: 31)\n\
  -p, --processes INT            Number of parallel worker threads (default: 1)"
    );
}

fn get_file_type(path: &Path) -> Option<FileType> {
    match path.extension().and_then(OsStr::to_str) {
        Some("fa" | "fas" | "fasta") => Some(FileType::Fasta),
        Some("fq" | "fastq" | "gm2") => Some(FileType::Fastq),
        _ => None,
    }
}

fn get_read_dict(
    se_dir: Option<&Path>,
    pe_dir: Option<&Path>,
    use_compressed: bool,
) -> Result<BTreeMap<String, Vec<PathBuf>>, String> {
    let mut read_dict = BTreeMap::new();

    if let Some(dir) = se_dir {
        if !dir.is_dir() {
            return Err("Argument --se-dir does not refer to a directory.".to_string());
        }

        for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
            let entry = entry.map_err(|e| e.to_string())?;
            let path = entry.path();

            if !path.is_file() || get_file_type(&path).is_none() {
                continue;
            }

            if !use_compressed && path.extension() == Some(OsStr::new("gm2")) {
                continue;
            }

            let basename = path.file_stem().and_then(OsStr::to_str).unwrap_or_default().to_string();
            insert_read_group(&mut read_dict, basename, vec![path])?;
        }
    }

    if let Some(dir) = pe_dir {
        if !dir.is_dir() {
            return Err("Argument --pe-dir does not refer to a directory.".to_string());
        }

        for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
            let entry = entry.map_err(|e| e.to_string())?;
            let path = entry.path();

            if !path.is_file() || get_file_type(&path).is_none() {
                continue;
            }

            if !use_compressed && path.extension() == Some(OsStr::new("gm2")) {
                continue;
            }

            let basename = path.file_stem().and_then(OsStr::to_str).unwrap_or_default();

            if !basename.ends_with("_1") {
                continue;
            }

            let gene_name = basename.trim_end_matches("_1").to_string();
            let ext = path.extension().and_then(OsStr::to_str).unwrap_or_default();
            let mate2 = dir.join(format!("{gene_name}_2.{ext}"));
            let paths = if mate2.is_file() { vec![path, mate2] } else { vec![path] };
            insert_read_group(&mut read_dict, gene_name, paths)?;
        }
    }

    Ok(read_dict)
}

fn insert_read_group(
    read_dict: &mut BTreeMap<String, Vec<PathBuf>>,
    name: String,
    paths: Vec<PathBuf>,
) -> Result<(), String> {
    if let Some(existing) = read_dict.get(&name) {
        let existing_is_gm2 = existing
            .first()
            .and_then(|p| p.extension())
            .is_some_and(|ext| ext == "gm2");
        let new_is_gm2 = paths
            .first()
            .and_then(|p| p.extension())
            .is_some_and(|ext| ext == "gm2");

        if existing_is_gm2 {
            return Ok(());
        }

        if !new_is_gm2 {
            return Err(format!("Duplicate read group name {name}."));
        }
    }

    read_dict.insert(name, paths);
    Ok(())
}

fn get_ref_dict(ref_dir: &Path) -> Result<BTreeMap<String, PathBuf>, String> {
    let mut ref_dict = BTreeMap::new();

    for entry in fs::read_dir(ref_dir).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();

        if !path.is_file() || get_file_type(&path).is_none() {
            continue;
        }

        let basename = path.file_stem().and_then(OsStr::to_str).unwrap_or_default().to_string();

        if ref_dict.insert(basename.clone(), path).is_some() {
            return Err(format!("Duplicate reference sequence name {basename}."));
        }
    }

    Ok(ref_dict)
}

fn run(
    args: Args,
    read_dict: BTreeMap<String, Vec<PathBuf>>,
    ref_dict: BTreeMap<String, PathBuf>,
) -> Result<(), String> {
    let mut tasks = VecDeque::new();

    for (name, ref_path) in ref_dict {
        let Some(read_paths) = read_dict.get(&name) else {
            print_log(args.log_file.as_deref(), &format!("No reads for gene {name}."))?;
            continue;
        };

        tasks.push_back(Task {
            name,
            out_dir: args.out_dir.clone(),
            ref_path,
            read_paths: read_paths.clone(),
            log_path: args.log_file.clone(),
            min_depth: args.min_depth,
            max_depth: args.max_depth,
            max_size: args.max_size,
            copy_only: args.copy_only,
            keep_temporaries: args.keep_temporaries,
            keep_linked_mates: args.keep_linked_mates,
            kmer_size: args.kmer_size,
        });
    }

    if tasks.is_empty() {
        print_log(args.log_file.as_deref(), "No genes with matching reads and references.")?;
        return Ok(());
    }

    let processes = args.processes.min(tasks.len());

    if processes > 1 {
        let queue = Arc::new(Mutex::new(tasks));
        let errors = Arc::new(Mutex::new(Vec::new()));

        thread::scope(|scope| {
            for _ in 0..processes {
                let queue = Arc::clone(&queue);
                let errors = Arc::clone(&errors);

                scope.spawn(move || loop {
                    let task = {
                        match queue.lock() {
                            Ok(mut queue) => queue.pop_front(),
                            Err(_) => {
                                if let Ok(mut errors) = errors.lock() {
                                    errors.push("task queue lock was poisoned".to_string());
                                }
                                None
                            }
                        }
                    };

                    let Some(task) = task else {
                        break;
                    };

                    if let Err(error) = filter_gene(task) {
                        if let Ok(mut errors) = errors.lock() {
                            errors.push(error);
                        }
                    }
                });
            }
        });

        let errors = errors.lock().map_err(|_| "error list poisoned".to_string())?;

        if let Some(error) = errors.first() {
            return Err(error.clone());
        }
    } else {
        for task in tasks {
            filter_gene(task)?;
        }
    }

    if !args.keep_temporaries {
        let _ = fs::remove_dir(args.out_dir.join("large_files"));
    }

    Ok(())
}

fn filter_gene(task: Task) -> Result<(), String> {
    let file_type = task
        .read_paths
        .first()
        .and_then(|path| get_file_type(path))
        .ok_or_else(|| format!("File '{}' has invalid file type.", task.read_paths[0].display()))?;

    if task.copy_only {
        print_log(task.log_path.as_deref(), &format!("Writing reads for gene {}.", task.name))?;
        copy_reads(&task.name, &task.out_dir, &task.read_paths, file_type)?;
        return Ok(());
    }

    print_log(task.log_path.as_deref(), &format!("Filtering gene {}.", task.name))?;
    let (ref_set, effective_len) = load_reference(&task.ref_path, task.kmer_size)?;

    if effective_len == 0.0 {
        print_log(task.log_path.as_deref(), &format!("Gene {} has no valid reference.", task.name))?;
        return Ok(());
    }

    let run_kmer = std::cmp::max(task.kmer_size / 2, task.kmer_size.saturating_sub(13)) | 1;
    let tmp_path = run_length_filter(
        &task.name,
        &task.out_dir,
        &ref_set,
        &task.read_paths,
        file_type,
        run_kmer,
        task.keep_linked_mates && task.read_paths.len() == 2,
    )?;

    kmer_filter(
        &task.name,
        &task.out_dir,
        task.log_path.as_deref(),
        &ref_set,
        effective_len,
        &tmp_path,
        file_type,
        task.kmer_size,
        task.min_depth,
        task.max_depth,
        task.max_size,
        task.keep_linked_mates && task.read_paths.len() == 2,
    )?;

    if !task.keep_temporaries {
        fs::remove_file(tmp_path).map_err(|e| e.to_string())?;
    }

    Ok(())
}

fn print_log(log_path: Option<&Path>, message: &str) -> Result<(), String> {
    if let Some(path) = log_path {
        let mut out = OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .map_err(|e| e.to_string())?;
        writeln!(out, "{message}").map_err(|e| e.to_string())?;
    }

    println!("{message}");
    Ok(())
}

fn load_reference(ref_path: &Path, kmer_size: usize) -> Result<(HashSet<String>, f64), String> {
    let mut reader = RecordReader::from_path(ref_path, FileType::Fasta, false, String::new())
        .map_err(|e| e.to_string())?;
    let mut ref_set = HashSet::new();

    while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
        if record.seq.len() >= kmer_size {
            ref_set.insert(record.seq);
        }
    }

    if ref_set.is_empty() {
        return Ok((ref_set, 0.0));
    }

    let max_len = ref_set.iter().map(|seq| seq.len()).max().unwrap_or(0) as f64;
    let effective_len = (max_len * ((ref_set.len() as f64).log10() + 1.0)).trunc();
    Ok((ref_set, effective_len))
}

fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(b'0'),
        b'C' | b'c' => Some(b'1'),
        b'G' | b'g' => Some(b'2'),
        b'T' | b't' | b'U' | b'u' => Some(b'3'),
        _ => None,
    }
}

fn translate_fwd(seq: &str) -> Vec<u8> {
    let mut encoded = Vec::with_capacity(seq.len());

    for base in seq.bytes() {
        let Some(base) = encode_base(base) else {
            // Never join the flanks around an ambiguity into a fake k-mer.
            return Vec::new();
        };
        encoded.push(base);
    }

    encoded
}

fn reference_runs(seq: &str) -> Vec<Vec<u8>> {
    let mut runs = Vec::new();
    let mut run = Vec::new();

    for base in seq.bytes() {
        if let Some(base) = encode_base(base) {
            run.push(base);
        } else if !run.is_empty() {
            runs.push(std::mem::take(&mut run));
        }
    }

    if !run.is_empty() {
        runs.push(run);
    }

    runs
}

fn reverse_complement(encoded: &[u8]) -> Vec<u8> {
    encoded
        .iter()
        .rev()
        .map(|base| b'3' - (*base - b'0'))
        .collect()
}

fn build_kmer_dict(ref_set: &HashSet<String>, kmer_size: usize) -> HashMap<Vec<u8>, u8> {
    let mut kmer_dict = HashMap::new();

    for seq in ref_set {
        for fwd in reference_runs(seq) {
            let rev = reverse_complement(&fwd);
            add_kmers(&mut kmer_dict, &fwd, kmer_size, 1);
            add_kmers(&mut kmer_dict, &rev, kmer_size, 2);
        }
    }

    kmer_dict
}

fn add_kmers(kmer_dict: &mut HashMap<Vec<u8>, u8>, seq: &[u8], kmer_size: usize, orient: u8) {
    if seq.len() < kmer_size {
        return;
    }

    for window in seq.windows(kmer_size) {
        let entry = kmer_dict.entry(window.to_vec()).or_insert(0);
        *entry |= orient;
    }
}

fn collect_runs_stats(read: &[u8], kmer_dict: &HashMap<Vec<u8>, u8>, kmer_size: usize) -> [usize; 13] {
    let mut results = [0_usize; 13];

    if read.len() < kmer_size {
        return results;
    }

    let mut curr_dir = 0_usize;
    let mut curr_len = 0_usize;
    results[12] = read.len() - kmer_size + 1;

    for window in read.windows(kmer_size) {
        let orient = *kmer_dict.get(window).unwrap_or(&0) as usize;

        if orient != curr_dir {
            if curr_len > results[curr_dir] {
                results[curr_dir] = curr_len;
            }

            results[curr_dir + 4] += 1;
            curr_dir = orient;
            curr_len = 0;
        }

        if curr_dir != 0 {
            curr_len += 1;
            results[curr_dir + 8] += 1;
        }
    }

    if curr_len > results[curr_dir] {
        results[curr_dir] = curr_len;
    }

    results[curr_dir + 4] += 1;
    results
}

fn filter_read(read: &[u8], kmer_dict: &HashMap<Vec<u8>, u8>, kmer_size: usize) -> bool {
    read.len() >= kmer_size && read.windows(kmer_size).any(|window| kmer_dict.contains_key(window))
}

fn is_close(a: f64, b: f64, abs_tol: f64) -> bool {
    let rel_tol = 1e-9;
    (a - b).abs() <= abs_tol.max(rel_tol * a.abs().max(b.abs()))
}

fn infer_orientation(stats: [usize; 13]) -> u8 {
    let fwd_l = stats[1] as f64;
    let rev_l = stats[2] as f64;
    let fwd_r = stats[5] as f64;
    let rev_r = stats[6] as f64;
    let fwd_n = stats[9] as f64;
    let rev_n = stats[10] as f64;
    let amb_n = stats[11] as f64;
    let tot_n = stats[12] as f64;

    if fwd_n <= 1.0 {
        if rev_n <= 1.0 {
            return if amb_n <= 1.0 { 0 } else { 3 };
        }

        return 2;
    } else if rev_n <= 1.0 {
        return 1;
    }

    let npr = 2.0 * fwd_n * rev_n;
    let nht = fwd_n + rev_n;
    let erc = npr / nht + 1.0;
    let vrn = npr * (npr - nht) / (nht * nht * (nht - 1.0));

    if (fwd_r + rev_r - erc) / vrn.sqrt() > -THR_1E5_1T {
        let ntt = if fwd_r > rev_r {
            fwd_n + fwd_r - rev_r
        } else {
            fwd_n + rev_r - fwd_r
        };
        let rex = fwd_n / ntt;

        if is_close(rex, 1.0, TOLERANCE)
            || (1.0 - rex) / (rex * (1.0 - rex) / ntt).sqrt() < THR_P95_2T
        {
            return 0;
        }
    }

    let erl = (tot_n.log2() + RUN_LEN_CONST).max(0.0) + 4.0;
    let mut orient = ((fwd_l > erl) as u8) + ((rev_l > erl) as u8) * 2;

    if orient != 3 {
        return orient;
    }

    let lpf = (1.0 / (1.0 - fwd_l + fwd_n)).ln().mul_add(1.0 / fwd_l, 0.0).exp();
    let lpr = (1.0 / (1.0 - rev_l + rev_n)).ln().mul_add(1.0 / rev_l, 0.0).exp();
    let fpz = is_close(lpf, 0.0, TOLERANCE);
    let rpz = is_close(lpr, 0.0, TOLERANCE);

    if fpz {
        orient = 2 - 2 * (rpz as u8);
    } else if rpz {
        orient = 1;
    } else if (is_close(lpf, 1.0, TOLERANCE) && is_close(lpr, 1.0, TOLERANCE))
        || (lpf - lpr).abs()
            / (lpf.powi(2) * (1.0 - lpf) / fwd_n + lpr.powi(2) * (1.0 - lpr) / rev_n).sqrt()
            < THR_P95_2T
    {
        orient = 0;
    }

    orient
}

fn make_readers(paths: &[PathBuf], file_type: FileType) -> Result<Vec<RecordReader>, String> {
    let gm2_format = paths
        .first()
        .and_then(|path| path.extension())
        .is_some_and(|ext| ext == "gm2");

    paths
        .iter()
        .enumerate()
        .map(|(i, path)| {
            let suffix = if gm2_format {
                format!("/{}", i + 1)
            } else {
                String::new()
            };
            RecordReader::from_path(path, file_type, gm2_format, suffix).map_err(|e| e.to_string())
        })
        .collect()
}

fn next_linked_read(readers: &mut [RecordReader]) -> Result<Option<Vec<Record>>, String> {
    let mut linked = Vec::with_capacity(readers.len());
    let mut ended = 0;

    for reader in readers.iter_mut() {
        match reader.next_record().map_err(|e| e.to_string())? {
            Some(record) => linked.push(record),
            None => ended += 1,
        }
    }

    if ended == 0 {
        Ok(Some(linked))
    } else if ended == readers.len() {
        Ok(None)
    } else {
        Err("paired input files contain different numbers of records".to_string())
    }
}

fn write_record<W: Write>(out: &mut W, record: &Record, file_type: FileType) -> Result<(), String> {
    match file_type {
        FileType::Fasta => writeln!(out, ">{}\n{}", record.title, record.seq),
        FileType::Fastq => writeln!(out, "@{}\n{}\n+\n{}", record.title, record.seq, record.qual),
    }
    .map_err(|e| e.to_string())
}

fn copy_reads(name: &str, out_dir: &Path, read_paths: &[PathBuf], file_type: FileType) -> Result<(), String> {
    let output_path = out_dir.join(format!("{}{}", name, file_type.output_ext()));
    let mut readers = make_readers(read_paths, file_type)?;
    let mut out = BufWriter::new(File::create(output_path).map_err(|e| e.to_string())?);

    while let Some(linked_reads) = next_linked_read(&mut readers)? {
        for record in linked_reads {
            write_record(&mut out, &record, file_type)?;
        }
    }

    Ok(())
}

fn run_length_filter(
    name: &str,
    out_dir: &Path,
    ref_set: &HashSet<String>,
    read_paths: &[PathBuf],
    file_type: FileType,
    kmer_size: usize,
    keep_linked_mates: bool,
) -> Result<PathBuf, String> {
    let output_path = out_dir.join("large_files").join(format!("{}{}", name, file_type.output_ext()));
    let kmer_dict = build_kmer_dict(ref_set, kmer_size);
    let mut readers = make_readers(read_paths, file_type)?;
    let mut out = BufWriter::new(File::create(&output_path).map_err(|e| e.to_string())?);

    while let Some(linked_reads) = next_linked_read(&mut readers)? {
        let mut orient = vec![0_u8; linked_reads.len()];

        for (i, record) in linked_reads.iter().enumerate() {
            let read = translate_fwd(&record.seq);
            orient[i] = infer_orientation(collect_runs_stats(&read, &kmer_dict, kmer_size));
        }

        if orient.len() == 2 && (1..=2).contains(&orient[0]) && orient[0] == orient[1] {
            continue;
        }

        if keep_linked_mates && linked_reads.len() == 2 && orient.iter().any(|value| *value != 0) {
            for record in &linked_reads {
                write_record(&mut out, record, file_type)?;
            }
        } else {
            for (record, value) in linked_reads.iter().zip(orient.iter()) {
                if *value != 0 {
                    write_record(&mut out, record, file_type)?;
                }
            }
        }
    }

    Ok(output_path)
}

fn next_temp_group(
    reader: &mut RecordReader,
    keep_linked_mates: bool,
) -> Result<Option<Vec<Record>>, String> {
    let Some(first) = reader.next_record().map_err(|e| e.to_string())? else {
        return Ok(None);
    };

    if !keep_linked_mates {
        return Ok(Some(vec![first]));
    }

    let second = reader
        .next_record()
        .map_err(|e| e.to_string())?
        .ok_or_else(|| "interleaved paired-read temporary file has an odd number of records".to_string())?;
    Ok(Some(vec![first, second]))
}

fn count_total_length_from_path(
    path: &Path,
    file_type: FileType,
    keep_linked_mates: bool,
    kmer_dict: Option<&HashMap<Vec<u8>, u8>>,
    kmer_size: usize,
) -> Result<usize, String> {
    let mut reader = RecordReader::from_path(path, file_type, false, String::new()).map_err(|e| e.to_string())?;
    let mut total = 0;

    while let Some(group) = next_temp_group(&mut reader, keep_linked_mates)? {
        let keep = match kmer_dict {
            Some(kmer_dict) => group.iter().any(|record| {
                filter_read(&translate_fwd(&record.seq), kmer_dict, kmer_size)
            }),
            None => true,
        };

        if keep {
            total += group.iter().map(|record| record.seq.len()).sum::<usize>();
        }
    }

    Ok(total)
}

#[allow(clippy::too_many_arguments)]
fn kmer_filter(
    name: &str,
    out_dir: &Path,
    log_path: Option<&Path>,
    ref_set: &HashSet<String>,
    ref_length: f64,
    temp_path: &Path,
    file_type: FileType,
    mut kmer_size: usize,
    min_depth: i64,
    max_depth: i64,
    max_size: i64,
    keep_linked_mates: bool,
) -> Result<(), String> {
    let output_path = out_dir.join(format!("{}{}", name, file_type.output_ext()));
    let mut total_length = count_total_length_from_path(temp_path, file_type, keep_linked_mates, None, kmer_size)?;
    let mut coverage = total_length as f64 / ref_length;
    let mut too_deep = coverage > max_depth as f64;
    let mut too_large = total_length / 1_000_000 > max_size as usize;

    if !too_deep && !too_large {
        fs::copy(temp_path, output_path).map_err(|e| e.to_string())?;
        return Ok(());
    }

    let min_depth = (min_depth as f64).min(max_depth as f64 / 4.0);
    let initial_kmer_size = kmer_size;

    while kmer_size < 64 && (too_deep || too_large) {
        let last_kmer_size = kmer_size;
        let last_length = total_length;

        if coverage > 8.0 * max_depth as f64 || total_length / 1_000_000 > (6 * max_size) as usize {
            kmer_size += 6;
        } else {
            kmer_size += 2;
        }

        print_log(log_path, &format!("K-mer size for {name}: {kmer_size}"))?;
        let kmer_dict = build_kmer_dict(ref_set, kmer_size);
        total_length = count_total_length_from_path(temp_path, file_type, keep_linked_mates, Some(&kmer_dict), kmer_size)?;
        coverage = total_length as f64 / ref_length;
        too_deep = coverage > max_depth as f64;
        too_large = total_length / 1_000_000 > max_size as usize;

        if coverage < min_depth {
            kmer_size = last_kmer_size;
            total_length = last_length;
            too_large = total_length / 1_000_000 > max_size as usize;
            break;
        }
    }

    if kmer_size == initial_kmer_size && !too_large {
        fs::copy(temp_path, output_path).map_err(|e| e.to_string())?;
        return Ok(());
    }

    let kmer_dict = build_kmer_dict(ref_set, kmer_size);
    let interval = std::cmp::max((total_length as f64 / 1e6 / max_size as f64).trunc() as usize, 2);
    let mut i = 0_usize;
    let mut out = BufWriter::new(File::create(output_path).map_err(|e| e.to_string())?);

    let mut reader = RecordReader::from_path(temp_path, file_type, false, String::new()).map_err(|e| e.to_string())?;
    while let Some(linked_reads) = next_temp_group(&mut reader, keep_linked_mates)? {
        if keep_linked_mates {
            if linked_reads
                .iter()
                .any(|record| filter_read(&translate_fwd(&record.seq), &kmer_dict, kmer_size))
            {
                i += 1;

                if too_large && !i.is_multiple_of(interval) {
                    continue;
                }

                for record in &linked_reads {
                    write_record(&mut out, record, file_type)?;
                }
            }
        }
        else if let Some(record) = linked_reads.first() {
            if filter_read(&translate_fwd(&record.seq), &kmer_dict, kmer_size) {
                i += 1;

                if too_large && !i.is_multiple_of(interval) {
                    continue;
                }

                write_record(&mut out, record, file_type)?;
            }
        }
    }

    Ok(())
}
