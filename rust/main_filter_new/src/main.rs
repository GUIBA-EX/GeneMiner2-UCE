use ahash::AHashMap;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::env;
use std::ffi::{c_char, c_int, c_void, CString};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

const CACHE_MAGIC: &[u8; 4] = b"GM2K";
const CACHE_VERSION: u16 = 3;
const FLAG_CANONICAL_KMERS: u16 = 1;
const MEBIBYTE_READS: u64 = 1_048_576;
const TOTAL_BUFFER_BUDGET: usize = 64 * 1024 * 1024;
const MIN_FILE_BUDGET: usize = 128 * 1024;

type AppResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FileKind {
    Fasta,
    Fastq,
}

#[derive(Clone, Debug)]
struct Args {
    reference: PathBuf,
    q1: Vec<PathBuf>,
    q2: Vec<PathBuf>,
    kmer: usize,
    step: usize,
    max_read_blocks: u64,
    output: PathBuf,
    out_subdir: String,
    dictionary: Option<PathBuf>,
    get_reverse: bool,
    use_composition_pattern: bool,
    mode: u8,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            reference: PathBuf::new(),
            q1: Vec::new(),
            q2: Vec::new(),
            kmer: 31,
            step: 3,
            max_read_blocks: 0,
            output: PathBuf::new(),
            out_subdir: "filtered".to_string(),
            dictionary: None,
            get_reverse: false,
            use_composition_pattern: false,
            mode: 0,
        }
    }
}

// 把能用的开关一次说明白，用户别靠猜参数整活儿。
fn print_help() {
    println!(
        "MainFilterNew (Rust implementation)\n\
         Usage: MainFilterNew -r REF -o OUT [options]\n\n\
         -r PATH       Reference FASTA file or directory\n\
         -q1 FILE...   Read 1 / single-end FASTA, FASTQ, or gzipped FASTQ\n\
         -q2 FILE...   Read 2 files\n\
         -kf INT       K-mer length (default: 31)\n\
         -s INT        Read scanning step (default: 3)\n\
         -m_reads INT  Maximum 2^20 read records per input file; 0 disables\n\
         -o PATH       Output directory\n\
         -subdir NAME  Filtered-read subdirectory (default: filtered)\n\
         -lkd PATH     Load or write a reusable k-mer dictionary\n\
         -gr           Index reverse-complement reference k-mers\n\
         -lb           Accept the legacy composition-pattern option\n\
         -m INT        Mode 0..5 (default: 0)\n\
         --version     Print version"
    );
}

// 取紧跟在开关后头的值；少了就当场吱声，别往下带病跑。
fn take_value(argv: &[String], index: &mut usize, option: &str) -> AppResult<String> {
    *index += 1;
    argv.get(*index)
        .cloned()
        .ok_or_else(|| format!("option {option} requires an argument"))
}

// 数字参数统一在这儿验，免得各处转来转去整乱套。
fn parse_usize(value: String, option: &str) -> AppResult<usize> {
    value
        .parse::<usize>()
        .map_err(|_| format!("cannot convert {value:?} for {option} to an integer"))
}

fn parse_args(argv: Vec<String>) -> AppResult<Args> {
    // 先把参数归拢明白，后面干活儿就不抓瞎了。
    if argv.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        process::exit(0);
    }
    if argv.iter().any(|arg| arg == "--version") {
        println!("MainFilterNew 0.4.0");
        process::exit(0);
    }

    let mut parsed = Args::default();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-r" => parsed.reference = PathBuf::from(take_value(&argv, &mut i, "-r")?),
            "-o" => parsed.output = PathBuf::from(take_value(&argv, &mut i, "-o")?),
            "-kf" => parsed.kmer = parse_usize(take_value(&argv, &mut i, "-kf")?, "-kf")?,
            "-s" => parsed.step = parse_usize(take_value(&argv, &mut i, "-s")?, "-s")?,
            "-m_reads" => {
                parsed.max_read_blocks = take_value(&argv, &mut i, "-m_reads")?
                    .parse::<u64>()
                    .map_err(|_| "-m_reads requires a non-negative integer".to_string())?;
            }
            "-subdir" => parsed.out_subdir = take_value(&argv, &mut i, "-subdir")?,
            "-lkd" => parsed.dictionary = Some(PathBuf::from(take_value(&argv, &mut i, "-lkd")?)),
            "-m" => {
                parsed.mode = take_value(&argv, &mut i, "-m")?
                    .parse::<u8>()
                    .map_err(|_| "-m requires an integer from 0 to 5".to_string())?;
            }
            "-gr" => parsed.get_reverse = true,
            "-lb" => parsed.use_composition_pattern = true,
            "-q1" | "-q2" => {
                let is_q1 = argv[i] == "-q1";
                let mut paths = Vec::new();
                while i + 1 < argv.len() && !argv[i + 1].starts_with('-') {
                    i += 1;
                    paths.push(PathBuf::from(&argv[i]));
                }
                if paths.is_empty() {
                    return Err(format!(
                        "option {} requires at least one file",
                        if is_q1 { "-q1" } else { "-q2" }
                    ));
                }
                if is_q1 {
                    parsed.q1.extend(paths);
                } else {
                    parsed.q2.extend(paths);
                }
            }
            option if option.starts_with('-') => return Err(format!("unknown option {option}")),
            value => return Err(format!("unexpected argument {value}")),
        }
        i += 1;
    }

    if parsed.reference.as_os_str().is_empty() {
        return Err("reference sequences must be supplied".to_string());
    }
    if parsed.output.as_os_str().is_empty() {
        return Err("output directory must be supplied".to_string());
    }
    if parsed.kmer < 16 {
        return Err("k-mer size must be at least 16".to_string());
    }
    if parsed.step == 0 {
        return Err("read scanning step must be at least 1".to_string());
    }
    if parsed.mode > 5 {
        return Err("mode must be between 0 and 5".to_string());
    }
    if !parsed.q2.is_empty() && parsed.q1.len() != parsed.q2.len() {
        return Err("-q1 and -q2 must contain the same number of files".to_string());
    }
    if parsed.mode != 2 && parsed.q1.is_empty() {
        return Err("at least one sequencing file is required".to_string());
    }
    Ok(parsed)
}

#[link(name = "z")]
extern "C" {
    fn gzopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn gzread(file: *mut c_void, buffer: *mut c_void, length: u32) -> c_int;
    fn gzclose(file: *mut c_void) -> c_int;
}

struct GzipReader {
    handle: *mut c_void,
}

impl GzipReader {
    fn open(path: &Path) -> io::Result<Self> {
        let path = CString::new(path.to_string_lossy().as_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
        let mode = CString::new("rb").expect("static string contains no NUL");
        let handle = unsafe { gzopen(path.as_ptr(), mode.as_ptr()) };
        if handle.is_null() {
            Err(io::Error::new(
                io::ErrorKind::NotFound,
                "cannot open gzip file",
            ))
        } else {
            Ok(Self { handle })
        }
    }
}

impl Read for GzipReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        let size = buffer.len().min(c_int::MAX as usize) as u32;
        let result = unsafe { gzread(self.handle, buffer.as_mut_ptr().cast(), size) };
        if result < 0 {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gzip decompression failed",
            ))
        } else {
            Ok(result as usize)
        }
    }
}

impl Drop for GzipReader {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { gzclose(self.handle) };
        }
    }
}

// 看文件头认 gzip，比光瞅后缀靠谱点儿。
fn is_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
}

fn sequence_extension(path: &Path) -> Option<String> {
    let base = if is_gzip(path) {
        path.file_stem().map(PathBuf::from)?
    } else {
        path.to_path_buf()
    };
    base.extension()
        .and_then(|value| value.to_str())
        .map(|value| value.to_ascii_lowercase())
}

// 先认清是 FASTA 还是 FASTQ，后面读法才不能串台。
fn detect_kind(path: &Path) -> AppResult<FileKind> {
    match sequence_extension(path).as_deref() {
        Some("fa" | "fas" | "fasta") => Ok(FileKind::Fasta),
        Some("fq" | "fastq") => Ok(FileKind::Fastq),
        _ if is_gzip(path) => Ok(FileKind::Fastq),
        _ => Err(format!(
            "unrecognized sequence file type: {}",
            path.display()
        )),
    }
}

// gzip 和普通文件都包成同一种读口，后头不用分两套流程。
fn open_input(path: &Path) -> io::Result<BufReader<Box<dyn Read>>> {
    let input: Box<dyn Read> = if is_gzip(path) {
        Box::new(GzipReader::open(path)?)
    } else {
        Box::new(File::open(path)?)
    };
    Ok(BufReader::new(input))
}

#[derive(Clone, Debug)]
struct Record {
    lines: Option<Vec<String>>,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
}

struct SequenceReader {
    input: BufReader<Box<dyn Read>>,
    kind: FileKind,
    pending_header: Option<String>,
    finished: bool,
    keep_text_lines: bool,
}

impl SequenceReader {
    fn open(path: &Path, kind: FileKind, keep_text_lines: bool) -> AppResult<Self> {
        Ok(Self {
            input: open_input(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?,
            kind,
            pending_header: None,
            finished: false,
            keep_text_lines,
        })
    }

    fn read_line(&mut self) -> AppResult<Option<String>> {
        let mut line = String::new();
        if self.input.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
            return Ok(None);
        }
        Ok(Some(line.trim().to_string()))
    }

    fn next_record(&mut self) -> AppResult<Option<Record>> {
        match self.kind {
            FileKind::Fasta => self.next_fasta(),
            FileKind::Fastq => self.next_fastq(),
        }
    }

    fn next_fasta(&mut self) -> AppResult<Option<Record>> {
        if self.finished {
            return Ok(None);
        }
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                match self.read_line()? {
                    Some(line) if line.starts_with('>') => break line,
                    Some(line) if line.is_empty() => continue,
                    Some(_) => return Err("FASTA sequence encountered before a header".to_string()),
                    None => {
                        self.finished = true;
                        return Ok(None);
                    }
                }
            }
        };

        let mut sequence = Vec::new();
        loop {
            match self.read_line()? {
                Some(line) if line.starts_with('>') => {
                    self.pending_header = Some(line);
                    break;
                }
                Some(line) => sequence.extend(line.bytes().map(|base| base.to_ascii_uppercase())),
                None => {
                    self.finished = true;
                    break;
                }
            }
        }
        Ok(Some(Record {
            lines: self
                .keep_text_lines
                .then(|| vec![header, String::from_utf8_lossy(&sequence).into_owned()]),
            sequence,
            quality: None,
        }))
    }

    fn next_fastq(&mut self) -> AppResult<Option<Record>> {
        let header = loop {
            match self.read_line()? {
                Some(line) if line.is_empty() => continue,
                Some(line) => break line,
                None => return Ok(None),
            }
        };
        let sequence_line = self
            .read_line()?
            .ok_or_else(|| "truncated FASTQ sequence".to_string())?;
        let plus = self
            .read_line()?
            .ok_or_else(|| "truncated FASTQ plus line".to_string())?;
        let quality_line = self
            .read_line()?
            .ok_or_else(|| "truncated FASTQ quality".to_string())?;
        if !header.starts_with('@') || !plus.starts_with('+') {
            return Err("malformed FASTQ record".to_string());
        }
        let sequence: Vec<u8> = sequence_line
            .bytes()
            .map(|base| base.to_ascii_uppercase())
            .collect();
        let quality = quality_line.into_bytes();
        if quality.len() != sequence.len() {
            return Err("FASTQ sequence and quality lengths differ".to_string());
        }
        Ok(Some(Record {
            lines: self.keep_text_lines.then(|| {
                vec![
                    header,
                    String::from_utf8_lossy(&sequence).into_owned(),
                    plus,
                    String::from_utf8_lossy(&quality).into_owned(),
                ]
            }),
            sequence,
            quality: Some(quality),
        }))
    }
}

const INVALID_BASE_CODE: u8 = u8::MAX;

const fn build_base_code_table() -> [u8; 256] {
    let mut table = [INVALID_BASE_CODE; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
}

const BASE_CODE_TABLE: [u8; 256] = build_base_code_table();

#[inline(always)]
// 碱基压成 2-bit；碰上 N 就返回空，窗口从这儿断开。
fn base_code(base: u8) -> Option<u8> {
    match BASE_CODE_TABLE[base as usize] {
        INVALID_BASE_CODE => None,
        code => Some(code),
    }
}

#[derive(Clone, Copy, Debug)]
enum ReferenceHits {
    One(u32),
    // Used only while building: each k-mer accumulates its locus IDs once.
    Pending(u32),
    // Final/read-only representation: a slice in packed_hits.
    Packed { offset: u32, len: u32 },
}

#[derive(Clone, Debug)]
enum KmerStore {
    Short(AHashMap<u64, ReferenceHits>),
    Long(AHashMap<Vec<u8>, ReferenceHits>),
}

impl KmerStore {
    fn new(k: usize) -> Self {
        if k <= 32 {
            Self::Short(AHashMap::new())
        } else {
            Self::Long(AHashMap::new())
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Short(map) => map.len(),
            Self::Long(map) => map.len(),
        }
    }
    fn reserve(&mut self, additional: usize) {
        match self {
            Self::Short(map) => map.reserve(additional),
            Self::Long(map) => map.reserve(additional),
        }
    }
}

#[derive(Clone, Debug)]
struct KmerIndex {
    k: usize,
    reference_names: Vec<String>,
    store: KmerStore,
    pending_hits: Vec<Vec<u32>>,
    packed_hits: Vec<u32>,
}

impl KmerIndex {
    fn new(k: usize, reference_names: Vec<String>) -> Self {
        Self {
            k,
            reference_names,
            store: KmerStore::new(k),
            pending_hits: Vec::new(),
            packed_hits: Vec::new(),
        }
    }
    fn len(&self) -> usize {
        self.store.len()
    }

    fn add_reference_hit(
        hit: &mut ReferenceHits,
        reference: u32,
        pending_hits: &mut Vec<Vec<u32>>,
    ) {
        match *hit {
            ReferenceHits::One(previous) if previous != reference => {
                let index = pending_hits.len() as u32;
                pending_hits.push(vec![previous, reference]);
                *hit = ReferenceHits::Pending(index);
            }
            ReferenceHits::Pending(index) => {
                let hits = &mut pending_hits[index as usize];
                if hits.last().copied() != Some(reference) {
                    hits.push(reference);
                }
            }
            ReferenceHits::One(_) => {}
            ReferenceHits::Packed { .. } => unreachable!("cannot extend finalized reference hits"),
        }
    }

    fn insert_short(&mut self, key: u64, reference: u32) {
        let KmerStore::Short(map) = &mut self.store else {
            unreachable!("short k-mer inserted into long-kmer index");
        };
        match map.entry(key) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(ReferenceHits::One(reference));
            }
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                Self::add_reference_hit(entry.get_mut(), reference, &mut self.pending_hits)
            }
        }
    }

    fn insert_long(&mut self, key: &[u8], reference: u32) {
        let KmerStore::Long(map) = &mut self.store else {
            unreachable!("long k-mer inserted into short-kmer index");
        };
        if let Some(hit) = map.get_mut(key) {
            Self::add_reference_hit(hit, reference, &mut self.pending_hits);
        } else {
            map.insert(key.to_vec(), ReferenceHits::One(reference));
        }
    }

    fn store_loaded_hits(&mut self, hits: Vec<u32>) -> AppResult<ReferenceHits> {
        if hits.is_empty() {
            return Err("dictionary k-mer has no reference hits".to_string());
        }
        if hits.len() == 1 {
            return Ok(ReferenceHits::One(hits[0]));
        }
        let offset = u32::try_from(self.packed_hits.len())
            .map_err(|_| "packed reference hits exceed u32 range".to_string())?;
        let len = u32::try_from(hits.len())
            .map_err(|_| "reference hit list exceeds u32 range".to_string())?;
        offset
            .checked_add(len)
            .ok_or_else(|| "packed reference hit range exceeds u32 range".to_string())?;
        self.packed_hits.extend(hits);
        Ok(ReferenceHits::Packed { offset, len })
    }

    fn pack_pending_hit(
        hit: &mut ReferenceHits,
        pending_hits: &[Vec<u32>],
        packed_hits: &mut Vec<u32>,
    ) -> AppResult<()> {
        if let ReferenceHits::Pending(index) = *hit {
            let values = &pending_hits[index as usize];
            let offset = u32::try_from(packed_hits.len())
                .map_err(|_| "packed reference hits exceed u32 range".to_string())?;
            let len = u32::try_from(values.len())
                .map_err(|_| "reference hit list exceeds u32 range".to_string())?;
            offset
                .checked_add(len)
                .ok_or_else(|| "packed reference hit range exceeds u32 range".to_string())?;
            packed_hits.extend_from_slice(values);
            *hit = ReferenceHits::Packed { offset, len };
        }
        Ok(())
    }

    fn finalize_hits(&mut self) -> AppResult<()> {
        match &mut self.store {
            KmerStore::Short(map) => {
                for hit in map.values_mut() {
                    Self::pack_pending_hit(hit, &self.pending_hits, &mut self.packed_hits)?;
                }
            }
            KmerStore::Long(map) => {
                for hit in map.values_mut() {
                    Self::pack_pending_hit(hit, &self.pending_hits, &mut self.packed_hits)?;
                }
            }
        }
        self.pending_hits.clear();
        self.pending_hits.shrink_to_fit();
        Ok(())
    }

    fn hit_slice<'a>(&'a self, hits: &'a ReferenceHits) -> &'a [u32] {
        match hits {
            ReferenceHits::One(reference) => std::slice::from_ref(reference),
            ReferenceHits::Packed { offset, len } => {
                let end = offset
                    .checked_add(*len)
                    .expect("validated packed reference-hit range");
                &self.packed_hits[*offset as usize..end as usize]
            }
            ReferenceHits::Pending(_) => {
                unreachable!("reference hits must be finalized before querying")
            }
        }
    }

    fn short_mask(&self) -> u64 {
        if self.k == 32 {
            u64::MAX
        } else {
            (1_u64 << (2 * self.k)) - 1
        }
    }

    fn add_reference_sequence(&mut self, sequence: &[u8], reference: u32) {
        // 参考序列滚着编码进索引，短 k-mer 不整 Vec，省内存老鼻子了。
        if sequence.len() < self.k {
            return;
        }
        if self.k <= 32 {
            let mask = self.short_mask();
            let reverse_shift = 2 * (self.k - 1);
            let mut forward = 0_u64;
            let mut reverse = 0_u64;
            let mut valid = 0_usize;
            for &base in sequence {
                if let Some(code) = base_code(base) {
                    forward = ((forward << 2) | code as u64) & mask;
                    reverse = (reverse >> 2) | (((3 - code) as u64) << reverse_shift);
                    valid += 1;
                    if valid >= self.k {
                        self.insert_short(forward.min(reverse), reference);
                    }
                } else {
                    forward = 0;
                    reverse = 0;
                    valid = 0;
                }
            }
        } else {
            let mut forward = Vec::with_capacity(self.k);
            let mut reverse = Vec::with_capacity(self.k);
            for start in 0..=sequence.len() - self.k {
                if long_kmer_into(&sequence[start..start + self.k], &mut forward, &mut reverse) {
                    self.insert_long(
                        if forward <= reverse {
                            &forward
                        } else {
                            &reverse
                        },
                        reference,
                    );
                }
            }
        }
    }

    fn collect_hits(&self, sequence: &[u8], step: usize, collector: &mut HitCollector) {
        // Canonical k-mers preserve bidirectional matching while one lookup per window suffices.
        if sequence.len() < self.k {
            return;
        }
        let tail = sequence.len() - self.k;
        if self.k <= 32 {
            let mask = self.short_mask();
            let reverse_shift = 2 * (self.k - 1);
            let mut forward = 0_u64;
            let mut reverse = 0_u64;
            let mut valid = 0_usize;
            let mut next_probe = 0_usize;
            for (end, &base) in sequence.iter().enumerate() {
                if let Some(code) = base_code(base) {
                    forward = ((forward << 2) | code as u64) & mask;
                    reverse = (reverse >> 2) | (((3 - code) as u64) << reverse_shift);
                    valid += 1;
                } else {
                    forward = 0;
                    reverse = 0;
                    valid = 0;
                }
                if end + 1 < self.k {
                    continue;
                }
                let start = end + 1 - self.k;
                let sampled = start == next_probe;
                if sampled {
                    next_probe = next_probe.saturating_add(step);
                }
                if valid >= self.k && (sampled || start == tail) {
                    self.collect_short(forward.min(reverse), collector);
                }
            }
        } else {
            let mut forward = Vec::with_capacity(self.k);
            let mut reverse = Vec::with_capacity(self.k);
            for start in (0..=tail)
                .step_by(step)
                .chain((!tail.is_multiple_of(step)).then_some(tail))
            {
                if long_kmer_into(&sequence[start..start + self.k], &mut forward, &mut reverse) {
                    self.collect_long(
                        if forward <= reverse {
                            &forward
                        } else {
                            &reverse
                        },
                        collector,
                    );
                }
            }
        }
    }

    fn collect_short(&self, key: u64, collector: &mut HitCollector) {
        let KmerStore::Short(map) = &self.store else {
            unreachable!("short lookup in long-kmer index");
        };
        if let Some(hits) = map.get(&key) {
            for &reference in self.hit_slice(hits) {
                collector.mark(reference as usize);
            }
        }
    }

    fn collect_long(&self, key: &[u8], collector: &mut HitCollector) {
        let KmerStore::Long(map) = &self.store else {
            unreachable!("long lookup in short-kmer index");
        };
        if let Some(hits) = map.get(key) {
            for &reference in self.hit_slice(hits) {
                collector.mark(reference as usize);
            }
        }
    }
}

// k 太长塞不进 u64 时，老老实实造正反链字符串键。
fn long_kmer_into(sequence: &[u8], forward: &mut Vec<u8>, reverse: &mut Vec<u8>) -> bool {
    forward.clear();
    reverse.clear();
    if forward.capacity() < sequence.len() {
        forward.reserve(sequence.len() - forward.capacity());
    }
    if reverse.capacity() < sequence.len() {
        reverse.reserve(sequence.len() - reverse.capacity());
    }
    for &base in sequence {
        let Some(code) = base_code(base) else {
            forward.clear();
            reverse.clear();
            return false;
        };
        forward.push(b"ACGT"[code as usize]);
    }
    for &base in sequence.iter().rev() {
        let Some(code) = base_code(base) else {
            unreachable!("validated above");
        };
        reverse.push(b"ACGT"[(3 - code) as usize]);
    }
    true
}

struct HitCollector {
    generation: u32,
    seen: Vec<u32>,
    hits: Vec<u32>,
}

impl HitCollector {
    fn new(reference_count: usize) -> Self {
        Self {
            generation: 0,
            seen: vec![0; reference_count],
            hits: Vec::new(),
        }
    }

    fn begin(&mut self) {
        // 换个批次号就当清空了，省得每条 read 都把整张表擦一遍。
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.seen.fill(0);
            self.generation = 1;
        }
        self.hits.clear();
    }

    fn mark(&mut self, reference: usize) {
        if self.seen[reference] != self.generation {
            self.seen[reference] = self.generation;
            self.hits.push(reference as u32);
        }
    }
}

// 文件名就是 locus 名，去掉 .gz 和序列后缀，别把扩展名带进结果。
fn reference_basename(path: &Path) -> AppResult<String> {
    let without_gzip = if is_gzip(path) {
        path.file_stem()
            .map(PathBuf::from)
            .ok_or_else(|| "invalid reference name".to_string())?
    } else {
        path.to_path_buf()
    };
    without_gzip
        .file_stem()
        .and_then(|value| value.to_str())
        .map(str::to_string)
        .ok_or_else(|| format!("invalid reference name: {}", path.display()))
}

// 参考既能给一个文件，也能给一摞文件；这儿统一排好队。
fn reference_paths(reference: &Path) -> AppResult<Vec<PathBuf>> {
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if path.is_file() && detect_kind(&path).ok() == Some(FileKind::Fasta) {
                paths.push(path);
            }
        }
    } else if reference.is_file() {
        paths.push(reference.to_path_buf());
    }
    paths.sort();
    if paths.is_empty() {
        return Err("no reference FASTA file found".to_string());
    }
    Ok(paths)
}

fn reference_content_hash(reference: &Path) -> AppResult<[u8; 32]> {
    let paths = reference_paths(reference)?;
    let mut hasher = Sha256::new();
    hasher.update(b"GM2-MainFilter-reference-v1\0");
    for path in paths {
        let name = reference_basename(&path)?;
        hasher.update((name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
        let mut reader = SequenceReader::open(&path, FileKind::Fasta, false)?;
        while let Some(record) = reader.next_record()? {
            hasher.update((record.sequence.len() as u64).to_le_bytes());
            hasher.update(&record.sequence);
        }
    }
    Ok(hasher.finalize().into())
}

fn build_index(reference: &Path, k: usize) -> AppResult<(KmerIndex, [u8; 32])> {
    // 首次构建时把内容哈希和索引合在一遍参考读取里，避免额外 I/O。
    let paths = reference_paths(reference)?;
    let names: Vec<String> = paths
        .iter()
        .map(|path| reference_basename(path))
        .collect::<AppResult<_>>()?;
    let mut unique = HashSet::new();
    for name in &names {
        if !unique.insert(name.clone()) {
            return Err(format!("duplicate reference locus name: {name}"));
        }
    }
    let mut hasher = Sha256::new();
    hasher.update(b"GM2-MainFilter-reference-v1\0");
    let mut index = KmerIndex::new(k, names);
    for (reference_id, path) in paths.iter().enumerate() {
        let name = &index.reference_names[reference_id];
        hasher.update((name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
        let mut reader = SequenceReader::open(path, FileKind::Fasta, false)?;
        while let Some(record) = reader.next_record()? {
            hasher.update((record.sequence.len() as u64).to_le_bytes());
            hasher.update(&record.sequence);
            index.add_reference_sequence(&record.sequence, reference_id as u32);
        }
    }
    index.finalize_hits()?;
    Ok((index, hasher.finalize().into()))
}

fn write_u16(out: &mut impl Write, value: u16) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}
fn write_u32(out: &mut impl Write, value: u32) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}
fn write_u64(out: &mut impl Write, value: u64) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}

fn write_dictionary(index: &KmerIndex, reference_hash: &[u8; 32], path: &Path) -> AppResult<()> {
    // 先写临时文件再改名，半道断了也不至于把老缓存整坏喽。
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let temporary = path.with_extension(format!("tmp.{}", process::id()));
    let mut out = File::create(&temporary).map_err(|e| e.to_string())?;
    out.write_all(CACHE_MAGIC).map_err(|e| e.to_string())?;
    write_u16(&mut out, CACHE_VERSION).map_err(|e| e.to_string())?;
    write_u16(&mut out, FLAG_CANONICAL_KMERS).map_err(|e| e.to_string())?;
    write_u32(&mut out, index.k as u32).map_err(|e| e.to_string())?;
    write_u32(&mut out, index.reference_names.len() as u32).map_err(|e| e.to_string())?;
    write_u64(&mut out, index.len() as u64).map_err(|e| e.to_string())?;
    out.write_all(reference_hash).map_err(|e| e.to_string())?;
    for name in &index.reference_names {
        write_u32(&mut out, name.len() as u32).map_err(|e| e.to_string())?;
        out.write_all(name.as_bytes()).map_err(|e| e.to_string())?;
    }
    match &index.store {
        KmerStore::Short(map) => {
            for (value, hits) in map {
                out.write_all(&[0]).map_err(|e| e.to_string())?;
                write_u64(&mut out, *value).map_err(|e| e.to_string())?;
                let hits = index.hit_slice(hits);
                write_u32(&mut out, hits.len() as u32).map_err(|e| e.to_string())?;
                for &hit in hits {
                    write_u32(&mut out, hit).map_err(|e| e.to_string())?;
                }
            }
        }
        KmerStore::Long(map) => {
            for (value, hits) in map {
                out.write_all(&[1]).map_err(|e| e.to_string())?;
                write_u32(&mut out, value.len() as u32).map_err(|e| e.to_string())?;
                out.write_all(value).map_err(|e| e.to_string())?;
                let hits = index.hit_slice(hits);
                write_u32(&mut out, hits.len() as u32).map_err(|e| e.to_string())?;
                for &hit in hits {
                    write_u32(&mut out, hit).map_err(|e| e.to_string())?;
                }
            }
        }
    }
    out.flush().map_err(|e| e.to_string())?;
    fs::rename(&temporary, path).map_err(|e| e.to_string())?;
    Ok(())
}

fn read_array<const N: usize>(input: &mut impl Read) -> AppResult<[u8; N]> {
    let mut bytes = [0_u8; N];
    input
        .read_exact(&mut bytes)
        .map_err(|_| "truncated k-mer dictionary".to_string())?;
    Ok(bytes)
}

fn read_u8(input: &mut impl Read) -> AppResult<u8> {
    Ok(read_array::<1>(input)?[0])
}

fn read_u16(input: &mut impl Read) -> AppResult<u16> {
    Ok(u16::from_le_bytes(read_array(input)?))
}

fn read_u32(input: &mut impl Read) -> AppResult<u32> {
    Ok(u32::from_le_bytes(read_array(input)?))
}

fn read_u64(input: &mut impl Read) -> AppResult<u64> {
    Ok(u64::from_le_bytes(read_array(input)?))
}

fn read_vec(input: &mut impl Read, length: usize) -> AppResult<Vec<u8>> {
    let mut bytes = vec![0_u8; length];
    input
        .read_exact(&mut bytes)
        .map_err(|_| "truncated k-mer dictionary".to_string())?;
    Ok(bytes)
}

// 优先尝试捡现成字典；格式和参数对不上就拒绝，不能瞎凑合。
fn load_dictionary(
    path: &Path,
    requested_k: usize,
    expected_reference_hash: &[u8; 32],
) -> AppResult<KmerIndex> {
    let file = File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut input = BufReader::new(file);
    let prefix = read_array::<4>(&mut input)?;
    if &prefix != CACHE_MAGIC {
        return Err(
            "legacy k-mer dictionary lacks a reference-content hash; rebuilding".to_string(),
        );
    }
    load_v3_dictionary(&mut input, requested_k, expected_reference_hash)
}

// 读取 Rust 缓存，核对 reference 内容、k 与 canonical 链策略，避免静默拿错字典。
fn load_v3_dictionary(
    input: &mut impl BufRead,
    requested_k: usize,
    expected_reference_hash: &[u8; 32],
) -> AppResult<KmerIndex> {
    if read_u16(input)? != CACHE_VERSION {
        return Err("unsupported k-mer dictionary version".to_string());
    }
    let flags = read_u16(input)?;
    if flags & FLAG_CANONICAL_KMERS == 0 {
        return Err("dictionary does not use canonical k-mers".to_string());
    }
    let k = read_u32(input)? as usize;
    if k != requested_k {
        return Err(format!(
            "requested {requested_k}-mer but dictionary contains {k}-mer"
        ));
    }
    let reference_count = read_u32(input)? as usize;
    let entry_count = usize::try_from(read_u64(input)?)
        .map_err(|_| "dictionary entry count exceeds this platform".to_string())?;
    let found_reference_hash = read_array::<32>(input)?;
    if &found_reference_hash != expected_reference_hash {
        return Err("dictionary reference-content hash does not match; rebuilding".to_string());
    }
    let mut names = Vec::with_capacity(reference_count);
    for _ in 0..reference_count {
        let length = read_u32(input)? as usize;
        names.push(
            String::from_utf8(read_vec(input, length)?)
                .map_err(|_| "reference name is not UTF-8".to_string())?,
        );
    }
    let mut index = KmerIndex::new(k, names);
    index.store.reserve(entry_count);
    for _ in 0..entry_count {
        let key_type = read_u8(input)?;
        let short_key = match key_type {
            0 if k <= 32 => Some(read_u64(input)?),
            1 => {
                let length = read_u32(input)? as usize;
                if k <= 32 || length != k {
                    return Err("dictionary k-mer type or length does not match k".to_string());
                }
                let key = read_vec(input, length)?;
                let hit_count = read_u32(input)? as usize;
                let hits = read_dictionary_hits(input, hit_count, reference_count)?;
                let value = index.store_loaded_hits(hits)?;
                let KmerStore::Long(map) = &mut index.store else {
                    unreachable!();
                };
                if map.insert(key, value).is_some() {
                    return Err("duplicate k-mer in dictionary".to_string());
                }
                None
            }
            _ => return Err("unknown k-mer key type".to_string()),
        };
        if let Some(key) = short_key {
            let hit_count = read_u32(input)? as usize;
            let hits = read_dictionary_hits(input, hit_count, reference_count)?;
            let value = index.store_loaded_hits(hits)?;
            let KmerStore::Short(map) = &mut index.store else {
                unreachable!();
            };
            if map.insert(key, value).is_some() {
                return Err("duplicate k-mer in dictionary".to_string());
            }
        }
    }
    let mut trailing = [0_u8; 1];
    if input.read(&mut trailing).map_err(|e| e.to_string())? != 0 {
        return Err("unexpected trailing dictionary data".to_string());
    }
    Ok(index)
}

fn read_dictionary_hits(
    input: &mut impl Read,
    hit_count: usize,
    reference_count: usize,
) -> AppResult<Vec<u32>> {
    let mut hits = Vec::with_capacity(hit_count);
    for _ in 0..hit_count {
        let hit = read_u32(input)?;
        if hit as usize >= reference_count {
            return Err("dictionary reference id is out of range".to_string());
        }
        if hits.last().copied().is_some_and(|last| hit <= last) {
            return Err("dictionary reference hits are not strictly increasing".to_string());
        }
        hits.push(hit);
    }
    Ok(hits)
}

struct OutputManager {
    paths: Vec<PathBuf>,
    buffers: Vec<Vec<u8>>,
    total_buffered: usize,
    file_budget: usize,
    paired_paths: bool,
    binary: bool,
    encode_scratch: Vec<u8>,
}

impl OutputManager {
    fn new(
        output: &Path,
        subdir: &str,
        names: &[String],
        kind: FileKind,
        mode: u8,
    ) -> AppResult<Self> {
        let directory = output.join(subdir);
        fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
        let extension = if kind == FileKind::Fasta {
            ".fasta"
        } else {
            ".fq"
        };
        let mut paths = Vec::new();
        match mode {
            0 => {
                for name in names {
                    paths.push(directory.join(format!("{name}{extension}")));
                }
            }
            1 => {
                paths.push(directory.join("all_1.fq"));
                paths.push(directory.join("all_2.fq"));
            }
            3 => {}
            4 => {
                for name in names {
                    paths.push(directory.join(format!("{name}_1{extension}")));
                    paths.push(directory.join(format!("{name}_2{extension}")));
                }
            }
            5 => {
                for name in names {
                    paths.push(directory.join(format!("{name}_1.gm2")));
                    paths.push(directory.join(format!("{name}_2.gm2")));
                }
            }
            _ => return Err("output manager does not support this mode".to_string()),
        }
        for path in &paths {
            if path.exists() {
                fs::remove_file(path).map_err(|e| e.to_string())?;
            }
        }
        let file_budget = if paths.is_empty() {
            MIN_FILE_BUDGET
        } else {
            (TOTAL_BUFFER_BUDGET / paths.len()).max(MIN_FILE_BUDGET)
        };
        let buffer_count = paths.len();
        Ok(Self {
            paths,
            buffers: vec![Vec::new(); buffer_count],
            total_buffered: 0,
            file_budget,
            paired_paths: mode == 4 || mode == 5,
            binary: mode == 5,
            encode_scratch: Vec::with_capacity(4096),
        })
    }

    fn write_encoded_record(
        &mut self,
        reference: usize,
        encoded: &[u8],
        second_mate: bool,
    ) -> AppResult<()> {
        let key = if self.paired_paths {
            2 * reference + usize::from(second_mate)
        } else {
            reference
        };
        self.write_encoded_key(key, encoded)
    }

    fn write_key(&mut self, key: usize, record: &Record) -> AppResult<()> {
        if self.paths.is_empty() {
            return Ok(());
        }
        encode_record_into(record, self.binary, &mut self.encode_scratch)?;
        let encoded = std::mem::take(&mut self.encode_scratch);
        let result = self.write_encoded_key(key, &encoded);
        self.encode_scratch = encoded;
        result
    }

    fn write_encoded_key(&mut self, key: usize, encoded: &[u8]) -> AppResult<()> {
        if self.paths.is_empty() {
            return Ok(());
        }
        self.total_buffered += encoded.len();
        self.buffers[key].extend_from_slice(encoded);
        if self.buffers[key].len() >= self.file_budget {
            self.flush_key(key)?;
        } else if self.total_buffered >= TOTAL_BUFFER_BUDGET {
            self.flush_largest()?;
        }
        Ok(())
    }

    fn flush_key(&mut self, key: usize) -> AppResult<()> {
        if self.buffers[key].is_empty() {
            return Ok(());
        }
        let mut out = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.paths[key])
            .map_err(|e| e.to_string())?;
        out.write_all(&self.buffers[key])
            .map_err(|e| e.to_string())?;
        self.total_buffered -= self.buffers[key].len();
        self.buffers[key].clear();
        Ok(())
    }

    fn flush_largest(&mut self) -> AppResult<()> {
        if let Some((key, _)) = self
            .buffers
            .iter()
            .enumerate()
            .max_by_key(|(_, buffer)| buffer.len())
        {
            self.flush_key(key)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> AppResult<()> {
        for key in 0..self.paths.len() {
            self.flush_key(key)?;
        }
        Ok(())
    }
}

// 文本输出保留原样，方便人眼复查。
fn encode_text_into(record: &Record, output: &mut Vec<u8>) {
    output.clear();
    for line in record
        .lines
        .as_ref()
        .expect("text output requires record lines")
    {
        output.extend_from_slice(line.as_bytes());
        output.push(b'\n');
    }
}

fn encode_record_into(record: &Record, binary: bool, output: &mut Vec<u8>) -> AppResult<()> {
    if binary {
        encode_gm2_into(record, output)
    } else {
        encode_text_into(record, output);
        Ok(())
    }
}

fn gm2_base_value(base: u8) -> u8 {
    let upper = base.to_ascii_uppercase();
    if upper.is_ascii_uppercase() {
        upper - 64
    } else {
        b'N' - 64
    }
}

fn append_gm2_sequence_chunk(output: &mut Vec<u8>, last_chunk: &mut u8, chunk: u8) {
    let delta = chunk ^ *last_chunk;
    if output.len() > 6 && output.last().copied() == Some(delta) {
        if let Some(last) = output.last_mut() {
            *last |= 0x80;
        }
    } else {
        output.push(delta);
    }
    *last_chunk = chunk;
}

// GM2 把序列压紧，磁盘慢的时候少搬点儿字节就挺顶用。
fn encode_gm2_into(record: &Record, output: &mut Vec<u8>) -> AppResult<()> {
    output.clear();
    if record.sequence.is_empty() {
        return Ok(());
    }
    if record.sequence.len() > 0x7f_ffff {
        return Err("GM2 sequence is too long".to_string());
    }
    let capacity = record
        .sequence
        .len()
        .checked_mul(2)
        .and_then(|size| size.checked_add(6))
        .ok_or_else(|| "GM2 record size overflow".to_string())?;
    if output.capacity() < capacity {
        output.reserve(capacity - output.capacity());
    }
    output.resize(6, 0);

    let mut last_chunk = 0_u8;
    let mut last_value = gm2_base_value(record.sequence[0]);
    let mut duplicate_count = 0_u8;
    for &base in &record.sequence[1..] {
        let value = gm2_base_value(base);
        if value != last_value || duplicate_count == 3 {
            append_gm2_sequence_chunk(output, &mut last_chunk, (duplicate_count << 5) | last_value);
            last_value = value;
            duplicate_count = 0;
        } else {
            duplicate_count += 1;
        }
    }
    append_gm2_sequence_chunk(output, &mut last_chunk, (duplicate_count << 5) | last_value);

    let has_quality = record
        .quality
        .as_ref()
        .is_some_and(|quality| quality.len() >= record.sequence.len());
    if let Some(quality) = record.quality.as_ref().filter(|_| has_quality) {
        let mut start = 0;
        while start < record.sequence.len() {
            let value = quality[start] & 0x7f;
            let mut length = 1_usize;
            while start + length < record.sequence.len()
                && quality[start + length] & 0x7f == value
                && length < 256
            {
                length += 1;
            }
            if length == 1 {
                output.push(value);
            } else {
                let duplicate_count = length - 1;
                output.push(0x80 | (duplicate_count as u8 & 0x7f));
                output.push(value | if duplicate_count & 0x80 != 0 { 0x80 } else { 0 });
            }
            start += length;
        }
    }
    let payload_length = output.len() - 6;
    if payload_length > 0xff_ffff {
        return Err("GM2 payload is too long".to_string());
    }
    output[..6].copy_from_slice(&[
        ((payload_length >> 16) & 0xff) as u8,
        ((payload_length >> 8) & 0xff) as u8,
        (payload_length & 0xff) as u8,
        (((record.sequence.len() >> 16) & 0x7f) as u8) | if has_quality { 0x80 } else { 0 },
        ((record.sequence.len() >> 8) & 0xff) as u8,
        (record.sequence.len() & 0xff) as u8,
    ]);
    Ok(())
}

#[cfg(test)]
fn encode_gm2(record: &Record) -> AppResult<Vec<u8>> {
    let mut output = Vec::new();
    encode_gm2_into(record, &mut output)?;
    Ok(output)
}

struct Logger {
    file: BufWriter<File>,
    last_flush: Instant,
}

impl Logger {
    fn new(output: &Path) -> AppResult<Self> {
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(output.join("log.txt"))
            .map_err(|e| e.to_string())?;
        Ok(Self {
            file: BufWriter::with_capacity(64 * 1024, file),
            last_flush: Instant::now(),
        })
    }
    fn log(&mut self, message: &str) {
        println!("{message}");
        let _ = writeln!(self.file, "{message}");
        if self.last_flush.elapsed().as_millis() >= 500 {
            let _ = self.file.flush();
            self.last_flush = Instant::now();
        }
    }
}
impl Drop for Logger {
    fn drop(&mut self) {
        let _ = self.file.flush();
    }
}

// 每个 locus 最后留个计数单子，方便看过滤到底捞着多少 reads。
fn write_read_counts(output: &Path, names: &[String], counts: &[u64]) -> AppResult<()> {
    let mut out =
        File::create(output.join("ref_reads_count_dict.txt")).map_err(|e| e.to_string())?;
    for (name, count) in names.iter().zip(counts) {
        if *count > 0 {
            writeln!(out, "{name},{count}").map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn run(args: Args) -> AppResult<()> {
    // 索引只建一次，随后按输入 read 逐条分到命中的 locus。
    fs::create_dir_all(&args.output).map_err(|e| e.to_string())?;
    let mut logger = Logger::new(&args.output)?;
    logger.log("Getting information from references...");
    let started = Instant::now();
    let cached = args
        .dictionary
        .as_ref()
        .filter(|path| path.exists())
        .map(|path| {
            let reference_hash = reference_content_hash(&args.reference)?;
            load_dictionary(path, args.kmer, &reference_hash).map(|index| (index, reference_hash))
        });
    let (index, reference_hash) = match cached {
        Some(Ok(value)) => {
            logger.log(&format!(
                "Loaded k-mer dictionary with {} entries.",
                value.0.len()
            ));
            value
        }
        Some(Err(reason)) => {
            logger.log(&format!("Ignoring reusable dictionary: {reason}"));
            let (index, hash) = build_index(&args.reference, args.kmer)?;
            logger.log(&format!(
                "Built k-mer dictionary with {} entries.",
                index.len()
            ));
            if let Some(path) = &args.dictionary {
                write_dictionary(&index, &hash, path)?;
            }
            (index, hash)
        }
        None => {
            let (index, hash) = build_index(&args.reference, args.kmer)?;
            logger.log(&format!(
                "Built k-mer dictionary with {} entries.",
                index.len()
            ));
            if let Some(path) = &args.dictionary {
                write_dictionary(&index, &hash, path)?;
            }
            (index, hash)
        }
    };
    let _ = reference_hash;
    logger.log(&format!(
        "Dictionary stage took {:.3} seconds.",
        started.elapsed().as_secs_f64()
    ));
    if args.mode == 2 {
        return Ok(());
    }
    if args.use_composition_pattern {
        logger
            .log("Note: -lb is accepted for compatibility; exact rolling k-mer scanning is used.");
    }

    let kind = detect_kind(&args.q1[0])?;
    for path in args.q1.iter().chain(args.q2.iter()) {
        if detect_kind(path)? != kind {
            return Err("all read files must use the same FASTA/FASTQ format".to_string());
        }
    }
    // GM2 和只扫描模式不留文本行，少分配几套没用的 String。
    let keep_text_lines = matches!(args.mode, 0 | 1 | 4);
    let mut output = OutputManager::new(
        &args.output,
        &args.out_subdir,
        &index.reference_names,
        kind,
        args.mode,
    )?;
    let mut counts = vec![0_u64; index.reference_names.len()];
    let mut collector = HitCollector::new(index.reference_names.len());
    let mut encoded1 = Vec::with_capacity(4096);
    let mut encoded2 = Vec::with_capacity(4096);
    if args.get_reverse {
        logger.log(
            "Note: -gr is retained for compatibility; canonical k-mers already match both strands.",
        );
    }
    let filter_started = Instant::now();

    for file_number in 0..args.q1.len() {
        let mut reader1 = SequenceReader::open(&args.q1[file_number], kind, keep_text_lines)?;
        let mut reader2 = if args.q2.is_empty() {
            None
        } else {
            Some(SequenceReader::open(
                &args.q2[file_number],
                kind,
                keep_text_lines,
            )?)
        };
        let mut read_count = 0_u64;
        let max_reads = args.max_read_blocks.saturating_mul(MEBIBYTE_READS);
        let mut stopped_at_limit = false;
        while let Some(record1) = reader1.next_record()? {
            let record2 = match reader2.as_mut() {
                Some(reader) => Some(reader.next_record()?.ok_or_else(|| {
                    "paired input files contain different numbers of records".to_string()
                })?),
                None => None,
            };
            collector.begin();
            index.collect_hits(&record1.sequence, args.step, &mut collector);
            if let Some(record) = &record2 {
                index.collect_hits(&record.sequence, args.step, &mut collector);
            }
            if !collector.hits.is_empty() {
                if args.mode == 1 {
                    output.write_key(0, &record1)?;
                    if let Some(record) = &record2 {
                        output.write_key(1, record)?;
                    }
                }
                let write_per_reference = matches!(args.mode, 0 | 4 | 5);
                if write_per_reference {
                    encode_record_into(&record1, args.mode == 5, &mut encoded1)?;
                    if let Some(record) = &record2 {
                        encode_record_into(record, args.mode == 5, &mut encoded2)?;
                    }
                }
                for &reference in &collector.hits {
                    let reference = reference as usize;
                    counts[reference] += if record2.is_some() { 2 } else { 1 };
                    if write_per_reference {
                        output.write_encoded_record(reference, &encoded1, false)?;
                    }
                    if write_per_reference && record2.is_some() {
                        output.write_encoded_record(reference, &encoded2, true)?;
                    }
                }
            }
            read_count += 1;
            if read_count.is_multiple_of(MEBIBYTE_READS) {
                logger.log(&format!(
                    "Handled {} Mi reads from {}.",
                    read_count / MEBIBYTE_READS,
                    args.q1[file_number].display()
                ));
            }
            if args.max_read_blocks > 0 && read_count >= max_reads {
                stopped_at_limit = true;
                break;
            }
        }
        if !stopped_at_limit {
            if let Some(reader) = reader2.as_mut() {
                if reader.next_record()?.is_some() {
                    return Err(
                        "paired input files contain different numbers of records".to_string()
                    );
                }
            }
        }
    }
    output.flush()?;
    write_read_counts(&args.output, &index.reference_names, &counts)?;
    logger.log(&format!(
        "Filtering took {:.3} seconds.",
        filter_started.elapsed().as_secs_f64()
    ));
    Ok(())
}

fn main() {
    let argv: Vec<String> = env::args().skip(1).collect();
    let args = match parse_args(argv) {
        Ok(args) => args,
        Err(error) => {
            eprintln!("Invalid argument: {error}");
            process::exit(2);
        }
    };
    if let Err(error) = run(args) {
        eprintln!("Error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn normalized_kmer(sequence: &[u8]) -> Option<Vec<u8>> {
        sequence
            .iter()
            .map(|&base| base_code(base).map(|code| b"ACGT"[code as usize]))
            .collect()
    }

    fn reverse_complement(sequence: &[u8]) -> Option<Vec<u8>> {
        sequence
            .iter()
            .rev()
            .map(|&base| base_code(base).map(|code| b"ACGT"[(3 - code) as usize]))
            .collect()
    }

    fn oracle_hits(
        references: &[Vec<u8>],
        read: &[u8],
        k: usize,
        step: usize,
        reverse_indexed: bool,
    ) -> Vec<u32> {
        if read.len() < k {
            return Vec::new();
        }
        let reference_kmers: Vec<HashSet<Vec<u8>>> = references
            .iter()
            .map(|reference| {
                let mut kmers = HashSet::new();
                for window in reference.windows(k) {
                    if let Some(forward) = normalized_kmer(window) {
                        if reverse_indexed {
                            kmers.insert(reverse_complement(&forward).unwrap());
                        }
                        kmers.insert(forward);
                    }
                }
                kmers
            })
            .collect();
        let tail = read.len() - k;
        let mut starts: Vec<usize> = (0..=tail).step_by(step).collect();
        if starts.last().copied() != Some(tail) {
            starts.push(tail);
        }
        let mut hits = Vec::new();
        for (reference, kmers) in reference_kmers.iter().enumerate() {
            let matched = starts.iter().any(|&start| {
                let window = &read[start..start + k];
                normalized_kmer(window).is_some_and(|forward| {
                    kmers.contains(&forward)
                        || (!reverse_indexed
                            && reverse_complement(&forward)
                                .is_some_and(|reverse| kmers.contains(&reverse)))
                })
            });
            if matched {
                hits.push(reference as u32);
            }
        }
        hits
    }

    fn pseudo_random_sequence(seed: &mut u64, length: usize) -> Vec<u8> {
        let alphabet = b"ACGTACGTACGTACGTACGTACGTACGTACGTUN";
        (0..length)
            .map(|_| {
                *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                alphabet[(*seed as usize) % alphabet.len()]
            })
            .collect()
    }

    #[test]
    fn ambiguous_bases_split_reference_kmers() {
        let mut index = KmerIndex::new(16, vec!["locus".to_string()]);
        index.add_reference_sequence(b"AAAAAAAAAAAAAAAANCCCCCCCCCCCCCCCC", 0);
        assert_eq!(index.len(), 2);
        let KmerStore::Short(map) = &index.store else {
            unreachable!();
        };
        assert!(!map.contains_key(&0b01));
    }

    #[test]
    fn reverse_complement_lookup_matches_original_modes() {
        let mut index = KmerIndex::new(16, vec!["locus".to_string()]);
        index.add_reference_sequence(b"AAAACCCCGGGGTTTA", 0);
        let mut hits = HitCollector::new(1);
        hits.begin();
        index.collect_hits(b"TAAACCCCGGGGTTTT", 3, &mut hits);
        assert_eq!(hits.hits, vec![0]);
    }

    #[test]
    fn gm2_header_uses_all_six_bytes() {
        let record = Record {
            lines: None,
            sequence: vec![b'A'; 300],
            quality: Some(vec![b'I'; 300]),
        };
        let encoded = encode_gm2(&record).unwrap();
        let payload =
            ((encoded[0] as usize) << 16) | ((encoded[1] as usize) << 8) | encoded[2] as usize;
        let sequence = (((encoded[3] & 0x7f) as usize) << 16)
            | ((encoded[4] as usize) << 8)
            | encoded[5] as usize;
        assert_eq!(encoded.len(), payload + 6);
        assert_eq!(sequence, 300);
        assert_ne!(encoded[3] & 0x80, 0);
    }

    #[test]
    fn randomized_lookup_matches_naive_oracle_across_boundaries() {
        let mut seed = 0x4d59_5df4_d0f3_3173;
        for &k in &[16, 31, 32, 33] {
            for &step in &[1, 3, k + 2] {
                for reverse_indexed in [false, true] {
                    for iteration in 0..8 {
                        let references: Vec<Vec<u8>> = (0..3)
                            .map(|_| pseudo_random_sequence(&mut seed, 96))
                            .collect();
                        let mut read = if iteration % 3 == 0 {
                            pseudo_random_sequence(&mut seed, 64)
                        } else {
                            references[iteration % references.len()][12..76].to_vec()
                        };
                        if iteration % 2 == 1 {
                            read = reverse_complement(&read).unwrap_or(read);
                        }
                        if iteration % 4 == 2 {
                            read[3] = b'N';
                        }

                        let mut index = KmerIndex::new(
                            k,
                            vec!["a".to_string(), "b".to_string(), "c".to_string()],
                        );
                        for (reference, sequence) in references.iter().enumerate() {
                            index.add_reference_sequence(sequence, reference as u32);
                        }
                        let mut collector = HitCollector::new(references.len());
                        collector.begin();
                        index.collect_hits(&read, step, &mut collector);
                        collector.hits.sort_unstable();

                        assert_eq!(
                            collector.hits,
                            oracle_hits(&references, &read, k, step, reverse_indexed),
                            "k={k}, step={step}, reverse_indexed={reverse_indexed}, iteration={iteration}"
                        );
                    }
                }
            }
        }
    }
}
