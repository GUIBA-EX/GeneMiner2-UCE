use ahash::{AHashMap, AHashSet};
use flate2::read::MultiGzDecoder;
use gm2_tools::fastx::{FastxFormat, FastxReader};
use rayon::prelude::*;
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BTreeSet, BinaryHeap, HashSet},
    env,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    process::Command,
};
type R<T> = Result<T, String>;
#[derive(Clone)]
struct S {
    taxon: String,
    id: String,
    r1: PathBuf,
    r2: Option<PathBuf>,
}
#[derive(Default)]
struct A {
    samples: PathBuf,
    out: PathBuf,
    stage: String,
    k: usize,
    min: u32,
    cap: usize,
    threads: usize,
    ledger: Option<PathBuf>,
    library: Option<PathBuf>,
    mainfilter: String,
    annotation_min_fragment: usize,
    annotation_max_fragment: usize,
    annotation_min_support: usize,
    annotation_min_identity: f64,
    annotation_min_coverage: f64,
    annotation_min_delta: f64,
    assemble_min_kmer_count: u32,
    assemble_branch_ratio: f64,
    assemble_max_fragments: usize,
    quantify_pairs: usize,
    bootstrap_replicates: usize,
    estimate_genome_fraction: bool,
}
fn val(v: &[String], i: &mut usize, n: &str) -> R<String> {
    *i += 1;
    v.get(*i)
        .cloned()
        .ok_or_else(|| format!("{n} needs a value"))
}
fn args() -> R<A> {
    let v: Vec<_> = env::args().skip(1).collect();
    let mut a = A {
        stage: "all".into(),
        k: 25,
        min: 8,
        cap: 10_000,
        threads: 1,
        annotation_min_fragment: 80,
        annotation_max_fragment: 800,
        annotation_min_support: 5,
        annotation_min_identity: 0.80,
        annotation_min_coverage: 0.60,
        annotation_min_delta: 0.10,
        assemble_min_kmer_count: 3,
        assemble_branch_ratio: 1.5,
        assemble_max_fragments: 3,
        quantify_pairs: 0,
        bootstrap_replicates: 200,
        ..A::default()
    };
    let mut i = 0;
    while i < v.len() {
        match v[i].as_str(){"--samples"=>a.samples=val(&v,&mut i,"--samples")?.into(),"--output"=>a.out=val(&v,&mut i,"--output")?.into(),"--stage"=>a.stage=val(&v,&mut i,"--stage")?,"--kmer"=>a.k=val(&v,&mut i,"--kmer")?.parse().map_err(|_|"bad --kmer")?,"--min-kmer-count"=>a.min=val(&v,&mut i,"--min-kmer-count")?.parse().map_err(|_|"bad --min-kmer-count")?,"--catalog-pairs"=>a.cap=val(&v,&mut i,"--catalog-pairs")?.parse().map_err(|_|"bad --catalog-pairs")?,"--threads"=>a.threads=val(&v,&mut i,"--threads")?.parse().map_err(|_|"bad --threads")?,"--read-ledger"=>a.ledger=Some(val(&v,&mut i,"--read-ledger")?.into()),"--te-library"=>a.library=Some(val(&v,&mut i,"--te-library")?.into()),"--mainfilter"=>a.mainfilter=val(&v,&mut i,"--mainfilter")?,"--annotation-min-fragment"=>a.annotation_min_fragment=val(&v,&mut i,"--annotation-min-fragment")?.parse().map_err(|_|"bad --annotation-min-fragment")?,"--annotation-max-fragment"=>a.annotation_max_fragment=val(&v,&mut i,"--annotation-max-fragment")?.parse().map_err(|_|"bad --annotation-max-fragment")?,"--annotation-min-support"=>a.annotation_min_support=val(&v,&mut i,"--annotation-min-support")?.parse().map_err(|_|"bad --annotation-min-support")?,"--annotation-min-identity"=>a.annotation_min_identity=val(&v,&mut i,"--annotation-min-identity")?.parse().map_err(|_|"bad --annotation-min-identity")?,"--annotation-min-coverage"=>a.annotation_min_coverage=val(&v,&mut i,"--annotation-min-coverage")?.parse().map_err(|_|"bad --annotation-min-coverage")?,"--annotation-min-delta"=>a.annotation_min_delta=val(&v,&mut i,"--annotation-min-delta")?.parse().map_err(|_|"bad --annotation-min-delta")?,"--assemble-min-kmer-count"=>a.assemble_min_kmer_count=val(&v,&mut i,"--assemble-min-kmer-count")?.parse().map_err(|_|"bad --assemble-min-kmer-count")?,"--assemble-branch-ratio"=>a.assemble_branch_ratio=val(&v,&mut i,"--assemble-branch-ratio")?.parse().map_err(|_|"bad --assemble-branch-ratio")?,"--assemble-max-fragments"=>a.assemble_max_fragments=val(&v,&mut i,"--assemble-max-fragments")?.parse().map_err(|_|"bad --assemble-max-fragments")?,"--quantify-pairs"=>a.quantify_pairs=val(&v,&mut i,"--quantify-pairs")?.parse().map_err(|_|"bad --quantify-pairs")?,"--bootstrap-replicates"=>a.bootstrap_replicates=val(&v,&mut i,"--bootstrap-replicates")?.parse().map_err(|_|"bad --bootstrap-replicates")?,"--estimate-genome-fraction"=>a.estimate_genome_fraction=true,"-h"|"--help"=>return Err("Usage: main_repeat --samples taxa.tsv --output DIR --mainfilter PATH [--stage all|discover|curate|annotate|quantify]".into()),x=>return Err(format!("unknown option {x}"))};
        i += 1
    }
    if a.samples.as_os_str().is_empty() || a.out.as_os_str().is_empty() {
        return Err("--samples and --output are required".into());
    }
    if !(16..=31).contains(&a.k) {
        return Err("--kmer must be 16..31".into());
    }
    if !matches!(
        a.stage.as_str(),
        "all" | "discover" | "curate" | "annotate" | "interspersed" | "quantify" | "compare"
    ) {
        return Err(
            "--stage must be all, discover, curate, annotate, interspersed, quantify, or compare"
                .into(),
        );
    }
    if a.threads == 0
        || a.cap == 0
        || a.min == 0
        || a.annotation_min_fragment < a.k
        || a.annotation_max_fragment < a.annotation_min_fragment
        || a.annotation_min_support == 0
        || !(0.0..=1.0).contains(&a.annotation_min_identity)
        || !(0.0..=1.0).contains(&a.annotation_min_coverage)
        || a.annotation_min_delta < 0.0
        || a.assemble_min_kmer_count == 0
        || a.assemble_branch_ratio < 1.0
        || a.bootstrap_replicates == 0
        || a.assemble_max_fragments == 0
        || a.assemble_max_fragments > 8
    {
        return Err("invalid annotation or discovery parameters".into());
    }
    Ok(a)
}
fn open(p: &Path) -> R<Box<dyn Read>> {
    let f = File::open(p).map_err(|e| e.to_string())?;
    if p.extension().is_some_and(|x| x == "gz") {
        Ok(Box::new(MultiGzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}
fn samples(p: &Path) -> R<Vec<S>> {
    let r = BufReader::new(open(p)?);
    let mut x = Vec::new();
    let mut first_data_row = true;
    for (n, l) in r.lines().enumerate() {
        let l = l.map_err(|e| e.to_string())?;
        if l.trim().is_empty() || l.starts_with('#') {
            continue;
        }
        let c: Vec<_> = l.split_whitespace().collect();
        if first_data_row
            && c.first()
                .is_some_and(|x| x.eq_ignore_ascii_case("taxon_id"))
        {
            first_data_row = false;
            continue;
        }
        first_data_row = false;
        let s = match c.len() {
            4 => S {
                taxon: c[0].into(),
                id: c[1].into(),
                r1: c[2].into(),
                r2: Some(c[3].into()),
            },
            3 => S {
                taxon: c[0].into(),
                id: c[1].into(),
                r1: c[2].into(),
                r2: None,
            },
            2 => S {
                taxon: c[0].into(),
                id: c[0].into(),
                r1: c[1].into(),
                r2: None,
            },
            _ => return Err(format!("bad samples line {}", n + 1)),
        };
        if s.id.is_empty()
            || Path::new(&s.id)
                .file_name()
                .is_none_or(|name| name != std::ffi::OsStr::new(&s.id))
            || s.id == "."
            || s.id == ".."
        {
            return Err(format!("unsafe sample_id on line {}", n + 1));
        }
        x.push(s)
    }
    if x.is_empty() {
        Err("no samples".into())
    } else {
        let mut ids = HashSet::new();
        for s in &x {
            if !ids.insert(&s.id) {
                return Err(format!("duplicate sample_id '{}'", s.id));
            }
        }
        Ok(x)
    }
}
fn id(h: &str) -> String {
    h.trim_start_matches('@')
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .into()
}
fn pair_hash(value: &str) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    for byte in value.bytes() {
        h = (h ^ byte as u64).wrapping_mul(0x100000001b3);
    }
    h
}
fn ledger(p: Option<&Path>) -> R<AHashMap<String, AHashSet<u64>>> {
    let mut o = AHashMap::new();
    let Some(p) = p else { return Ok(o) };
    for l in BufReader::new(open(p)?).lines() {
        let l = l.map_err(|e| e.to_string())?;
        let c: Vec<_> = l.split_whitespace().collect();
        if c.len() >= 2 && c[0] != "sample_id" {
            o.entry(c[0].into())
                .or_insert_with(AHashSet::new)
                .insert(pair_hash(c[1]));
        }
    }
    Ok(o)
}
#[derive(Clone)]
struct Rec {
    id: String,
    seq: Vec<u8>,
    qual: String,
}
struct PairReader {
    r1: FastxReader,
    r2: Option<FastxReader>,
}
fn read_fastq(r: &mut FastxReader) -> R<Option<Rec>> {
    let Some(record) = r.next_record().map_err(|e| e.to_string())? else {
        return Ok(None);
    };
    Ok(Some(Rec {
        id: id(&String::from_utf8_lossy(&record.header)),
        seq: record.sequence,
        qual: String::from_utf8_lossy(&record.quality).into_owned(),
    }))
}
impl PairReader {
    fn new(r1: &Path, r2: Option<&Path>) -> R<Self> {
        Ok(Self {
            r1: FastxReader::open(r1, FastxFormat::Fastq).map_err(|e| e.to_string())?,
            r2: r2
                .map(|path| FastxReader::open(path, FastxFormat::Fastq).map_err(|e| e.to_string()))
                .transpose()?,
        })
    }
    fn next_pair(&mut self) -> R<Option<(Rec, Option<Rec>)>> {
        let one = read_fastq(&mut self.r1)?;
        match one {
            None => {
                if let Some(r) = self.r2.as_mut() {
                    if read_fastq(r)?.is_some() {
                        return Err("R2 has more records than R1".into());
                    }
                }
                Ok(None)
            }
            Some(one) => {
                let two = match self.r2.as_mut() {
                    Some(r) => {
                        let two = read_fastq(r)?.ok_or_else(|| "R2 ended before R1".to_string())?;
                        if two.id != one.id {
                            return Err("paired FASTQ read IDs differ".into());
                        }
                        Some(two)
                    }
                    None => None,
                };
                Ok(Some((one, two)))
            }
        }
    }
}
fn enc(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
fn kmers(s: &[u8], k: usize, mut f: impl FnMut(u64)) {
    let mut x = 0;
    let mut y = 0;
    let mut n = 0;
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    for &b in s {
        if let Some(z) = enc(b) {
            x = ((x << 2) | z) & mask;
            y = (y >> 2) | ((3 - z) << (2 * (k - 1)));
            n += 1;
            if n >= k {
                f(x.min(y))
            }
        } else {
            x = 0;
            y = 0;
            n = 0
        }
    }
}
fn seq(mut x: u64, k: usize) -> String {
    let mut b = vec![b'A'; k];
    for i in (0..k).rev() {
        b[i] = b"ACGT"[(x & 3) as usize];
        x >>= 2
    }
    String::from_utf8(b).unwrap()
}
fn period(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    (1..=b.len() / 2)
        .find(|&p| (0..b.len()).filter(|&i| b[i] == b[i % p]).count() * 100 >= 85 * b.len())
}
#[derive(Eq, PartialEq)]
struct Ranked {
    rank: u64,
    r1: Vec<u8>,
    r2: Option<Vec<u8>>,
}
impl Ord for Ranked {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.rank.cmp(&o.rank)
    }
}
impl PartialOrd for Ranked {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
fn count_one(s: &S, cap: usize, k: usize, deny: &AHashSet<u64>) -> R<AHashMap<u64, u32>> {
    let mut reader = PairReader::new(&s.r1, s.r2.as_deref())?;
    let mut heap = BinaryHeap::new();
    while let Some((r1, r2)) = reader.next_pair()? {
        let rank = pair_hash(&r1.id);
        if deny.contains(&rank) {
            continue;
        }
        let entry = Ranked {
            rank,
            r1: r1.seq,
            r2: r2.map(|x| x.seq),
        };
        if heap.len() < cap {
            heap.push(entry)
        } else if rank < heap.peek().unwrap().rank {
            heap.pop();
            heap.push(entry)
        }
    }
    let mut m = AHashMap::new();
    for x in heap {
        kmers(&x.r1, k, |z| *m.entry(z).or_insert(0) += 1);
        if let Some(y) = x.r2 {
            kmers(&y, k, |z| *m.entry(z).or_insert(0) += 1)
        }
    }
    Ok(m)
}
struct Dsu {
    p: Vec<usize>,
}
impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            p: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.p[x] != x {
            let r = self.find(self.p[x]);
            self.p[x] = r
        }
        self.p[x]
    }
    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.p[b] = a
        }
    }
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}
fn canonical_motif(m: &str) -> String {
    let rc_m = rc(m);
    (0..m.len())
        .flat_map(|i| {
            [
                format!("{}{}", &m[i..], &m[..i]),
                format!("{}{}", &rc_m[i..], &rc_m[..i]),
            ]
        })
        .min()
        .unwrap_or_default()
}
fn file_hash(path: &Path) -> R<u64> {
    let mut r = open(path)?;
    let mut buf = [0_u8; 65536];
    let mut h = 0xcbf29ce484222325_u64;
    loop {
        let n = r.read(&mut buf).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        for &b in &buf[..n] {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3)
        }
    }
    Ok(h)
}
fn mix_hash(mut state: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        state = (state ^ byte as u64).wrapping_mul(0x100000001b3);
    }
    state
}
fn manifest_fields(path: &Path) -> R<BTreeMap<String, String>> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    let header: Vec<_> = lines
        .next()
        .ok_or_else(|| "malformed manifest".to_string())?
        .split('\t')
        .collect();
    let row: Vec<_> = lines
        .next()
        .ok_or_else(|| "malformed manifest".to_string())?
        .split('\t')
        .collect();
    if header.is_empty() || header.len() != row.len() || lines.next().is_some() {
        return Err("malformed manifest".into());
    }
    let mut fields = BTreeMap::new();
    for (key, value) in header.into_iter().zip(row) {
        if fields.insert(key.into(), value.into()).is_some() {
            return Err("malformed manifest".into());
        }
    }
    Ok(fields)
}
fn manifest_value<'a>(fields: &'a BTreeMap<String, String>, key: &str) -> R<&'a str> {
    fields
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| format!("malformed manifest: missing {key}"))
}
fn manifest_hash(fields: &BTreeMap<String, String>, key: &str) -> R<u64> {
    u64::from_str_radix(manifest_value(fields, key)?, 16)
        .map_err(|_| format!("malformed manifest: invalid {key}"))
}
fn manifest_number<T: std::str::FromStr>(fields: &BTreeMap<String, String>, key: &str) -> R<T> {
    manifest_value(fields, key)?
        .parse()
        .map_err(|_| format!("malformed manifest: invalid {key}"))
}
fn tree_hash(root: &Path) -> R<u64> {
    fn visit(path: &Path, files: &mut Vec<PathBuf>) -> R<()> {
        for entry in fs::read_dir(path).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if path.is_dir() {
                visit(&path, files)?;
            } else if path.is_file() {
                files.push(path);
            }
        }
        Ok(())
    }
    let mut files = Vec::new();
    visit(root, &mut files)?;
    files.sort();
    let mut state = 0xcbf29ce484222325_u64;
    for path in files {
        let name = path.strip_prefix(root).map_err(|e| e.to_string())?;
        state = mix_hash(state, name.to_string_lossy().as_bytes());
        state = mix_hash(state, &file_hash(&path)?.to_le_bytes());
    }
    Ok(state)
}
fn discovery_artifact_hash(root: &Path) -> R<u64> {
    let mut state = 0xcbf29ce484222325_u64;
    for name in ["atomic_catalog.tsv", "atomic_seeds.fasta"] {
        state = mix_hash(state, name.as_bytes());
        state = mix_hash(state, &file_hash(&root.join(name))?.to_le_bytes());
    }
    state = mix_hash(state, &tree_hash(&root.join("seeds"))?.to_le_bytes());
    Ok(state)
}
fn curated_artifact_hash(root: &Path) -> R<u64> {
    let mut state = 0xcbf29ce484222325_u64;
    for name in [
        "candidate_samples.tsv",
        "atomic_repeat_signal.tsv",
        "equivalence_map.tsv",
        "curated_catalog.tsv",
        "curation_evidence.tsv",
        "repeat_linkage.tsv",
        "topology.tsv",
    ] {
        state = mix_hash(state, name.as_bytes());
        state = mix_hash(state, &file_hash(&root.join(name))?.to_le_bytes());
    }
    for dir in ["library", "candidate_recruit"] {
        state = mix_hash(state, dir.as_bytes());
        state = mix_hash(state, &tree_hash(&root.join(dir))?.to_le_bytes());
    }
    Ok(state)
}
fn validate_discover(a: &A) -> R<()> {
    let root = a.out.join("01_discover");
    let fields = manifest_fields(&root.join("manifest.tsv")).map_err(|_| {
        "missing or malformed discovery manifest; rerun --stage discover".to_string()
    })?;
    let ledger_hash = a.ledger.as_deref().map(file_hash).transpose()?.unwrap_or(0);
    let library_hash = a
        .library
        .as_deref()
        .map(file_hash)
        .transpose()?
        .unwrap_or(0);
    if manifest_value(&fields, "schema_version")? != "1"
        || manifest_value(&fields, "stage")? != "discover"
        || manifest_number::<usize>(&fields, "kmer")? != a.k
        || manifest_number::<u32>(&fields, "min_kmer_count")? != a.min
        || manifest_number::<usize>(&fields, "catalog_pairs")? != a.cap
        || manifest_hash(&fields, "sample_manifest_hash")? != file_hash(&a.samples)?
        || manifest_hash(&fields, "ledger_hash")? != ledger_hash
        || manifest_hash(&fields, "te_library_hash")? != library_hash
        || manifest_hash(&fields, "artifact_hash")? != discovery_artifact_hash(&root)?
    {
        return Err(
            "discovery output does not match current inputs or artifacts; rerun --stage discover"
                .into(),
        );
    }
    Ok(())
}
fn discover(a: &A, ss: &[S], led: &AHashMap<String, AHashSet<u64>>) -> R<()> {
    let labels: AHashMap<u64, String> = AHashMap::new();
    let final_d = a.out.join("01_discover");
    let d = a.out.join(format!(".discover.tmp.{}", std::process::id()));
    if d.exists() {
        fs::remove_dir_all(&d).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(d.join("seeds")).map_err(|e| e.to_string())?;
    let mut by = BTreeMap::new();
    for s in ss {
        *by.entry(s.taxon.clone()).or_insert(0usize) += 1
    }
    let maps: Vec<_> = ss
        .par_iter()
        .map(|s| {
            count_one(
                s,
                (a.cap / by[&s.taxon]).max(1),
                a.k,
                led.get(&s.id).unwrap_or(&AHashSet::new()),
            )
        })
        .collect::<Result<_, _>>()?;
    let mut all = AHashMap::new();
    for m in maps {
        for (x, n) in m {
            *all.entry(x).or_insert(0) += n
        }
    }
    let mut v: Vec<_> = all.into_iter().filter(|(_, n)| *n >= a.min).collect();
    v.sort_by_key(|x| Reverse(x.1));
    v.truncate(2000);
    let raws: Vec<String> = v.iter().map(|x| seq(x.0, a.k)).collect();
    let mut dsu = Dsu::new(raws.len());
    let mut pref: AHashMap<String, Vec<usize>> = AHashMap::new();
    let mut suff: AHashMap<String, Vec<usize>> = AHashMap::new();
    let mut motif: AHashMap<String, usize> = AHashMap::new();
    for (i, raw) in raws.iter().enumerate() {
        for z in [raw.clone(), rc(raw)] {
            let pre = z[..a.k - 1].to_string();
            let suf = z[1..].to_string();
            pref.entry(pre).or_default().push(i);
            suff.entry(suf).or_default().push(i);
        }
        if let Some(p) = period(raw) {
            let key = canonical_motif(&raw[..p]);
            if let Some(&j) = motif.get(&key) {
                dsu.join(i, j)
            } else {
                motif.insert(key, i);
            }
        }
    }
    for (key, left) in &suff {
        if let Some(right) = pref.get(key) {
            let left_unique: AHashSet<_> = left.iter().copied().collect();
            let right_unique: AHashSet<_> = right.iter().copied().collect();
            if left_unique.len() == 1 && right_unique.len() == 1 {
                dsu.join(
                    *left_unique.iter().next().unwrap(),
                    *right_unique.iter().next().unwrap(),
                );
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..raws.len() {
        let r = dsu.find(i);
        groups.entry(r).or_default().push(i)
    }
    let mut meta =
        BufWriter::new(File::create(d.join("atomic_catalog.tsv")).map_err(|e| e.to_string())?);
    writeln!(
        meta,
        "repeat_id\tclass\tseed_count\tkmer_count\ttentative_annotation"
    )
    .unwrap();
    let mut fa =
        BufWriter::new(File::create(d.join("atomic_seeds.fasta")).map_err(|e| e.to_string())?);
    for (g, members) in groups.values().enumerate() {
        let rid = format!("R{:05}", g + 1);
        let is_tandem = members.iter().any(|&i| period(&raws[i]).is_some());
        let class = if is_tandem {
            "tandem-like"
        } else {
            "dispersed-repeat"
        };
        let mut f = BufWriter::new(
            File::create(d.join("seeds").join(format!("{rid}.fasta")))
                .map_err(|e| e.to_string())?,
        );
        let mut total = 0;
        let mut hits: AHashMap<String, usize> = AHashMap::new();
        for &i in members {
            let raw = &raws[i];
            let seed = if let Some(p) = period(raw) {
                raw[..p].repeat((3 * a.k).div_ceil(p))
            } else {
                raw.clone()
            };
            writeln!(f, ">{rid}_{i}\n{seed}").unwrap();
            writeln!(fa, ">{rid}_{i}\n{seed}").unwrap();
            total += v[i].1;
            if let Some(x) = labels.get(&v[i].0) {
                *hits.entry(x.clone()).or_insert(0) += 1
            }
        }
        let annotation = hits
            .into_iter()
            .max_by_key(|x| x.1)
            .filter(|x| x.1 * 2 >= members.len())
            .map(|x| x.0)
            .unwrap_or_else(|| ".".into());
        writeln!(
            meta,
            "{rid}\t{class}\t{}\t{total}\t{annotation}",
            members.len()
        )
        .unwrap();
    }
    drop(meta);
    drop(fa);
    let artifact_hash = discovery_artifact_hash(&d)?;
    let mut manifest =
        BufWriter::new(File::create(d.join("manifest.tsv")).map_err(|e| e.to_string())?);
    writeln!(manifest, "schema_version\tstage\tkmer\tmin_kmer_count\tcatalog_pairs\tsample_count\tsample_manifest_hash\tledger_hash\tte_library_hash\tartifact_hash").unwrap();
    writeln!(
        manifest,
        "1\tdiscover\t{}\t{}\t{}\t{}\t{:016x}\t{:016x}\t{:016x}\t{artifact_hash:016x}",
        a.k,
        a.min,
        a.cap,
        ss.len(),
        file_hash(&a.samples)?,
        a.ledger.as_deref().map(file_hash).transpose()?.unwrap_or(0),
        a.library
            .as_deref()
            .map(file_hash)
            .transpose()?
            .unwrap_or(0),
    )
    .unwrap();
    let backup = a
        .out
        .join(format!(".discover.previous.{}", std::process::id()));
    if backup.exists() {
        fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
    }
    if final_d.exists() {
        fs::rename(&final_d, &backup).map_err(|e| e.to_string())?;
    }
    if let Err(error) = fs::rename(&d, &final_d) {
        if backup.exists() {
            let _ = fs::rename(&backup, &final_d);
        }
        return Err(error.to_string());
    }
    if backup.exists() {
        fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
    }
    Ok(())
}
type SeedIndex = AHashMap<u64, Vec<usize>>;
type Adjacency = AHashMap<Vec<u8>, [u32; 4]>;

fn seed_index_dir(a: &A, dir: &Path) -> R<(Vec<String>, SeedIndex)> {
    let mut files: Vec<_> = fs::read_dir(dir)
        .map_err(|e| e.to_string())?
        .filter_map(|x| x.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fasta"))
        .collect();
    files.sort();
    let mut names = Vec::new();
    let mut idx = AHashMap::new();
    for p in files {
        let id = p.file_stem().unwrap().to_string_lossy().to_string();
        let n = names.len();
        let mut sequence = String::new();
        for line in BufReader::new(open(&p)?).lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with('>') {
                if !sequence.is_empty() {
                    kmers(sequence.as_bytes(), a.k, |x| {
                        idx.entry(x).or_insert_with(Vec::new).push(n)
                    });
                    sequence.clear();
                }
            } else {
                sequence.push_str(line.trim());
            }
        }
        if !sequence.is_empty() {
            kmers(sequence.as_bytes(), a.k, |x| {
                idx.entry(x).or_insert_with(Vec::new).push(n)
            });
        }
        names.push(id);
    }
    Ok((names, idx))
}
fn eligible(s: &S, deny: &AHashSet<u64>, out: &Path) -> R<(PathBuf, Option<PathBuf>, u64)> {
    if deny.is_empty() {
        let mut reader = PairReader::new(&s.r1, s.r2.as_deref())?;
        let mut n = 0;
        while reader.next_pair()?.is_some() {
            n += 1;
        }
        return Ok((s.r1.clone(), s.r2.clone(), n));
    }
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let p1 = out.join(format!("{}_1.fq", s.id));
    let p2 = s.r2.as_ref().map(|_| out.join(format!("{}_2.fq", s.id)));
    let mut w1 = BufWriter::new(File::create(&p1).map_err(|e| e.to_string())?);
    let mut w2 = match &p2 {
        Some(p) => Some(BufWriter::new(File::create(p).map_err(|e| e.to_string())?)),
        None => None,
    };
    let mut reader = PairReader::new(&s.r1, s.r2.as_deref())?;
    let mut n = 0;
    while let Some((one, two)) = reader.next_pair()? {
        if deny.contains(&pair_hash(&one.id)) {
            continue;
        }
        writeln!(
            w1,
            "@{}\n{}\n+\n{}",
            one.id,
            String::from_utf8_lossy(&one.seq),
            one.qual
        )
        .unwrap();
        if let (Some(w), Some(two)) = (w2.as_mut(), two) {
            writeln!(
                w,
                "@{}\n{}\n+\n{}",
                two.id,
                String::from_utf8_lossy(&two.seq),
                two.qual
            )
            .unwrap();
        }
        n += 1;
    }
    Ok((p1, p2, n))
}
#[derive(Default)]
struct Audit {
    specific: Vec<u64>,
    ambiguous: Vec<u64>,
    mapped_bases: Vec<u64>,
    covered_kmers: Vec<AHashSet<u64>>,
    identity_sum: Vec<f64>,
    identity_n: Vec<u64>,
    // Percent-divergence bins 0..=20, derived only from anchored reads.
    identity_bins: Vec<[u64; 21]>,
    links: BTreeMap<(usize, usize), (u64, u64)>,
}
fn inverted_repeat_identity(sequence: &[u8]) -> f64 {
    if sequence.len() < 200 {
        return 0.0;
    }
    let reversed = rc(&String::from_utf8_lossy(sequence));
    diagonal_similarity(sequence, reversed.as_bytes(), 0).0
}
fn groups_for(seq: &[u8], k: usize, index: &AHashMap<u64, Vec<usize>>) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    kmers(seq, k, |x| {
        if let Some(v) = index.get(&x) {
            out.extend(v.iter().copied());
        }
    });
    out
}
fn transitions_for(
    seq: &[u8],
    k: usize,
    index: &AHashMap<u64, Vec<usize>>,
) -> BTreeSet<(usize, usize)> {
    let mut edges = BTreeSet::new();
    let mut last: Option<usize> = None;
    kmers(seq, k, |x| {
        let Some(v) = index.get(&x) else {
            last = None;
            return;
        };
        let groups: BTreeSet<_> = v.iter().copied().collect();
        if groups.len() != 1 {
            last = None;
            return;
        }
        let current = *groups.iter().next().unwrap();
        if let Some(previous) = last {
            if previous != current {
                edges.insert((previous.min(current), previous.max(current)));
            }
        }
        last = Some(current);
    });
    edges
}
fn audit(
    path1: &Path,
    path2: Option<&Path>,
    k: usize,
    index: &AHashMap<u64, Vec<usize>>,
    count: usize,
) -> R<Audit> {
    let mut out = Audit {
        specific: vec![0; count],
        ambiguous: vec![0; count],
        mapped_bases: vec![0; count],
        covered_kmers: (0..count).map(|_| AHashSet::new()).collect(),
        identity_sum: vec![0.0; count],
        identity_n: vec![0; count],
        identity_bins: vec![[0; 21]; count],
        ..Audit::default()
    };
    if !path1.exists() {
        return Ok(out);
    }
    let mut reader = PairReader::new(path1, path2)?;
    while let Some((one, two)) = reader.next_pair()? {
        let left = groups_for(&one.seq, k, index);
        let right = two
            .as_ref()
            .map(|x| groups_for(&x.seq, k, index))
            .unwrap_or_default();
        let all: BTreeSet<_> = left.union(&right).copied().collect();
        if all.len() == 1 {
            let group = *all.iter().next().unwrap();
            out.specific[group] += 1;
            let left_maps = left.contains(&group);
            let right_maps = right.contains(&group);
            if left_maps {
                out.mapped_bases[group] += one.seq.len() as u64;
            }
            if right_maps {
                out.mapped_bases[group] += two.as_ref().map(|x| x.seq.len() as u64).unwrap_or(0);
            }
            if left_maps {
                kmers(&one.seq, k, |x| {
                    if index
                        .get(&x)
                        .is_some_and(|v| v.iter().all(|&owner| owner == group))
                    {
                        out.covered_kmers[group].insert(x);
                    }
                });
            }
            if right_maps {
                if let Some(mate) = two.as_ref() {
                    kmers(&mate.seq, k, |x| {
                        if index
                            .get(&x)
                            .is_some_and(|v| v.iter().all(|&owner| owner == group))
                        {
                            out.covered_kmers[group].insert(x);
                        }
                    });
                }
            }
        } else if all.len() > 1 {
            for i in all {
                out.ambiguous[i] += 1;
            }
        }
        if left.len() == 1 && right.len() == 1 {
            let a = *left.iter().next().unwrap();
            let b = *right.iter().next().unwrap();
            if a != b {
                out.links.entry((a.min(b), a.max(b))).or_default().0 += 1;
            }
        }
        let mut edges = transitions_for(&one.seq, k, index);
        if let Some(two) = two.as_ref() {
            edges.extend(transitions_for(&two.seq, k, index));
        }
        for edge in edges {
            out.links.entry(edge).or_default().1 += 1;
        }
    }
    Ok(out)
}
struct FragmentRecord {
    group: usize,
    seq: Vec<u8>,
}
struct FragmentMapIndex {
    records: Vec<FragmentRecord>,
    positions: AHashMap<u64, Vec<(usize, usize)>>,
}
fn fragment_map_index(dir: &Path, groups: &[String], k: usize) -> R<FragmentMapIndex> {
    let mut records: Vec<FragmentRecord> = Vec::new();
    let mut positions: AHashMap<u64, Vec<(usize, usize)>> = AHashMap::new();
    for (group, name) in groups.iter().enumerate() {
        for (_, seq) in fasta_records(&dir.join(format!("{name}.fasta")))? {
            if seq.len() < k {
                continue;
            }
            let record = records.len();
            records.push(FragmentRecord { group, seq });
            let mut pos = 0;
            kmers(&records[record].seq, k, |x| {
                let hits = positions.entry(x).or_default();
                // Cap positions per fragment: alternative paths retain anchors.
                if hits.iter().filter(|&&(other, _)| other == record).count() < 8 {
                    hits.push((record, pos));
                }
                pos += 1;
            });
        }
    }
    Ok(FragmentMapIndex { records, positions })
}
fn anchored_alignment_one(
    seq: &[u8],
    group: usize,
    map: &FragmentMapIndex,
    k: usize,
) -> Option<(usize, f64)> {
    let mut offsets: AHashMap<(usize, isize), u32> = AHashMap::new();
    let mut read_pos = 0;
    kmers(seq, k, |x| {
        if let Some(hits) = map.positions.get(&x) {
            for &(record, ref_pos) in hits {
                if map.records[record].group == group {
                    *offsets
                        .entry((record, ref_pos as isize - read_pos as isize))
                        .or_default() += 1;
                }
            }
        }
        read_pos += 1;
    });
    let mut candidates: Vec<_> = offsets
        .into_iter()
        .map(|((record, offset), support)| (record, offset, support))
        .collect();
    candidates.sort_by(|left, right| {
        right
            .2
            .cmp(&left.2)
            .then_with(|| left.0.cmp(&right.0))
            .then_with(|| left.1.cmp(&right.1))
    });
    candidates
        .into_iter()
        .take(32)
        .filter_map(|(record, offset, support)| {
            let (identity, coverage) = diagonal_similarity(seq, &map.records[record].seq, offset);
            let aligned = (coverage * seq.len() as f64).round() as usize;
            (aligned >= k && identity >= 0.80)
                .then_some((aligned, identity, support, record, offset))
        })
        .max_by(|left, right| {
            left.2
                .cmp(&right.2)
                .then_with(|| left.1.partial_cmp(&right.1).unwrap())
                .then_with(|| left.0.cmp(&right.0))
                .then_with(|| right.3.cmp(&left.3))
                .then_with(|| right.4.cmp(&left.4))
        })
        .map(|x| (x.0, x.1))
}
fn anchored_alignment(
    seq: &[u8],
    group: usize,
    map: &FragmentMapIndex,
    k: usize,
) -> Option<(usize, f64)> {
    let forward = anchored_alignment_one(seq, group, map, k);
    let reverse =
        anchored_alignment_one(rc(&String::from_utf8_lossy(seq)).as_bytes(), group, map, k);
    match (forward, reverse) {
        (Some(a), Some(b)) => Some(if a.1 >= b.1 { a } else { b }),
        (Some(hit), None) | (None, Some(hit)) => Some(hit),
        (None, None) => None,
    }
}
fn selected_pair_hashes(
    sample: &S,
    deny: &AHashSet<u64>,
    limit: usize,
) -> R<Option<AHashSet<u64>>> {
    if limit == 0 {
        return Ok(None);
    }
    let mut selected = BinaryHeap::new();
    let mut reader = PairReader::new(&sample.r1, sample.r2.as_deref())?;
    while let Some((one, _)) = reader.next_pair()? {
        let hash = pair_hash(&one.id);
        if deny.contains(&hash) {
            continue;
        }
        if selected.len() < limit {
            selected.push(hash);
        } else if selected.peek().is_some_and(|largest| hash < *largest) {
            selected.pop();
            selected.push(hash);
        }
    }
    Ok(Some(selected.into_iter().collect()))
}

fn next_random(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / ((1_u64 << 53) as f64)
}

fn bootstrap_rpm_ci(specific: u64, effective: u64, replicates: usize, seed: u64) -> (f64, f64) {
    if effective == 0 {
        return (0.0, 0.0);
    }
    let p = specific as f64 / effective as f64;
    let sigma = (p * (1.0 - p) / effective as f64).sqrt();
    let mut state = seed.max(1);
    let mut values = Vec::with_capacity(replicates);
    for _ in 0..replicates {
        let u1 = next_random(&mut state).max(f64::MIN_POSITIVE);
        let u2 = next_random(&mut state);
        let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
        values.push((p + sigma * z).clamp(0.0, 1.0) * 1e6);
    }
    values.sort_by(|left, right| left.partial_cmp(right).unwrap());
    let low = values[((replicates - 1) * 25) / 1000];
    let high = values[((replicates - 1) * 975) / 1000];
    (low, high)
}

fn audit_sample(
    sample: &S,
    deny: &AHashSet<u64>,
    k: usize,
    index: &AHashMap<u64, Vec<usize>>,
    fragments: &FragmentMapIndex,
    count: usize,
    selected: Option<&AHashSet<u64>>,
) -> R<(Audit, u64)> {
    let mut out = Audit {
        specific: vec![0; count],
        ambiguous: vec![0; count],
        mapped_bases: vec![0; count],
        covered_kmers: (0..count).map(|_| AHashSet::new()).collect(),
        identity_sum: vec![0.0; count],
        identity_n: vec![0; count],
        identity_bins: vec![[0; 21]; count],
        ..Audit::default()
    };
    let mut reader = PairReader::new(&sample.r1, sample.r2.as_deref())?;
    let mut effective = 0_u64;
    while let Some((one, two)) = reader.next_pair()? {
        let hash = pair_hash(&one.id);
        if deny.contains(&hash) || selected.is_some_and(|ids| !ids.contains(&hash)) {
            continue;
        }
        effective += 1;
        let left = groups_for(&one.seq, k, index);
        let right = two
            .as_ref()
            .map(|x| groups_for(&x.seq, k, index))
            .unwrap_or_default();
        let all: BTreeSet<_> = left.union(&right).copied().collect();
        if all.len() == 1 {
            let group = *all.iter().next().unwrap();
            let left_alignment = left
                .contains(&group)
                .then(|| anchored_alignment(&one.seq, group, fragments, k))
                .flatten();
            let right_alignment = right
                .contains(&group)
                .then(|| {
                    two.as_ref()
                        .and_then(|mate| anchored_alignment(&mate.seq, group, fragments, k))
                })
                .flatten();
            if left_alignment.is_some() || right_alignment.is_some() {
                out.specific[group] += 1;
            }
            for (read, alignment) in [
                (&one.seq, left_alignment),
                (
                    two.as_ref().map(|x| &x.seq).unwrap_or(&Vec::new()),
                    right_alignment,
                ),
            ] {
                if let Some((aligned, identity)) = alignment {
                    out.mapped_bases[group] += aligned as u64;
                    out.identity_sum[group] += identity;
                    out.identity_n[group] += 1;
                    let bin = ((1.0 - identity).max(0.0) * 100.0).round() as usize;
                    out.identity_bins[group][bin.min(20)] += 1;
                    kmers(read, k, |x| {
                        if index
                            .get(&x)
                            .is_some_and(|v| v.iter().all(|&owner| owner == group))
                        {
                            out.covered_kmers[group].insert(x);
                        }
                    });
                }
            }
        } else if all.len() > 1 {
            for group in all {
                out.ambiguous[group] += 1;
            }
        }
    }
    Ok((out, effective))
}
fn fragment_kmer_counts(
    dir: &Path,
    k: usize,
    index: &SeedIndex,
    groups: &[String],
) -> R<AHashMap<String, usize>> {
    let mut out = AHashMap::new();
    for (group, id) in groups.iter().enumerate() {
        let mut seen = AHashSet::new();
        for (_, sequence) in fasta_records(&dir.join(format!("{id}.fasta")))? {
            kmers(&sequence, k, |x| {
                if index
                    .get(&x)
                    .is_some_and(|owners| owners.iter().all(|&owner| owner == group))
                {
                    seen.insert(x);
                }
            });
        }
        out.insert(id.clone(), seen.len());
    }
    Ok(out)
}

fn fragment_lengths(dir: &Path) -> R<AHashMap<String, usize>> {
    let mut out = AHashMap::new();
    for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
        let path = entry.map_err(|e| e.to_string())?.path();
        if path.extension().is_some_and(|x| x == "fasta") {
            let id = path.file_stem().unwrap().to_string_lossy().to_string();
            let length = fasta_records(&path)?
                .iter()
                .map(|x| x.1.len())
                .max()
                .unwrap_or(0);
            out.insert(id, length);
        }
    }
    Ok(out)
}
fn state(n: u64, specific: u64, ambiguous: u64) -> &'static str {
    if n < 100 {
        "NA"
    } else if specific >= 3 && specific as f64 / (specific + ambiguous).max(1) as f64 >= 0.70 {
        "PRESENT"
    } else if specific == 0 {
        "ABSENT"
    } else {
        "NA"
    }
}
#[derive(Default, Clone)]
struct AtomicMeta {
    class: String,
    seed_count: u64,
    kmer_count: u64,
    annotation: String,
}
fn read_atomic_meta(a: &A) -> R<BTreeMap<String, AtomicMeta>> {
    let mut out = BTreeMap::new();
    for (n, line) in fs::read_to_string(a.out.join("01_discover/atomic_catalog.tsv"))
        .map_err(|e| e.to_string())?
        .lines()
        .enumerate()
    {
        if n == 0 {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() >= 5 {
            out.insert(
                f[0].into(),
                AtomicMeta {
                    class: f[1].into(),
                    seed_count: f[2].parse().unwrap_or(0),
                    kmer_count: f[3].parse().unwrap_or(0),
                    annotation: f[4].into(),
                },
            );
        }
    }
    Ok(out)
}
fn fasta_key(path: &Path) -> R<String> {
    let mut records = BTreeSet::new();
    let mut sequence = String::new();
    for line in BufReader::new(open(path)?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('>') {
            if !sequence.is_empty() {
                let clean: String = sequence
                    .chars()
                    .filter(|x| matches!(x, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
                    .map(|x| x.to_ascii_uppercase())
                    .collect();
                let reverse = rc(&clean);
                records.insert(if clean <= reverse { clean } else { reverse });
                sequence.clear();
            }
        } else {
            sequence.push_str(line.trim());
        }
    }
    if !sequence.is_empty() {
        let clean: String = sequence
            .chars()
            .filter(|x| matches!(x, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
            .map(|x| x.to_ascii_uppercase())
            .collect();
        let reverse = rc(&clean);
        records.insert(clean.min(reverse));
    }
    Ok(records.into_iter().collect::<Vec<_>>().join("|"))
}
fn entropy2(path: &Path) -> R<f64> {
    let key = fasta_key(path)?;
    let b = key.as_bytes();
    if b.len() < 2 {
        return Ok(0.0);
    }
    let mut counts = [0_u64; 16];
    let mut total = 0_u64;
    for pair in b.windows(2) {
        if let (Some(a), Some(c)) = (enc(pair[0]), enc(pair[1])) {
            counts[(a * 4 + c) as usize] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return Ok(0.0);
    }
    Ok(counts
        .iter()
        .filter(|&&n| n > 0)
        .map(|&n| {
            let p = n as f64 / total as f64;
            -p * p.log2()
        })
        .sum())
}
fn commit(temp: &Path, final_d: &Path, label: &str) -> R<()> {
    let backup = final_d.with_file_name(format!(".{label}.previous.{}", std::process::id()));
    if backup.exists() {
        fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
    }
    if final_d.exists() {
        fs::rename(final_d, &backup).map_err(|e| e.to_string())?;
    }
    if let Err(e) = fs::rename(temp, final_d) {
        if backup.exists() {
            let _ = fs::rename(&backup, final_d);
        }
        return Err(e.to_string());
    }
    if backup.exists() {
        fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
    }
    Ok(())
}
fn curate(a: &A, ss: &[S], led: &AHashMap<String, AHashSet<u64>>) -> R<()> {
    if a.mainfilter.is_empty() {
        return Err("--mainfilter is required for curate".into());
    }
    validate_discover(a)?;
    let discover = a.out.join("01_discover");
    let (atomic, index) = seed_index_dir(a, &discover.join("seeds"))?;
    if atomic.is_empty() {
        return Err("discovery has no atomic seeds".into());
    }
    let final_d = a.out.join("02_curate");
    let temp = a.out.join(format!(".curate.tmp.{}", std::process::id()));
    if temp.exists() {
        fs::remove_dir_all(&temp).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(temp.join("candidate_recruit")).map_err(|e| e.to_string())?;
    fs::create_dir_all(temp.join("library")).map_err(|e| e.to_string())?;
    let mut signal = BufWriter::new(
        File::create(temp.join("atomic_repeat_signal.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(signal, "sample_id\ttaxon_id\tatomic_id\teffective_pairs\tspecific_pairs\tambiguous_pairs\tsignal_rpm\tstate").unwrap();
    let mut samples_out = BufWriter::new(
        File::create(temp.join("candidate_samples.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(samples_out, "sample_id\ttaxon_id\teffective_pairs\tpaired").unwrap();
    let mut all_links: BTreeMap<(usize, usize), (u64, u64)> = BTreeMap::new();
    let mut support = vec![(0_u64, 0_u64, 0_u64); atomic.len()];
    for s in ss {
        let (q1, q2, effective) = eligible(
            s,
            led.get(&s.id).unwrap_or(&AHashSet::new()),
            &temp.join("eligible"),
        )?;
        let rdir = temp.join("candidate_recruit").join(&s.id);
        fs::create_dir_all(&rdir).map_err(|e| e.to_string())?;
        let mut command = Command::new(&a.mainfilter);
        command.args([
            "-r",
            discover.join("seeds").to_str().unwrap(),
            "-q1",
            q1.to_str().unwrap(),
            "-o",
            rdir.to_str().unwrap(),
            "-kf",
            &a.k.to_string(),
            "-s",
            "1",
            "-gr",
            "-m",
            "1",
        ]);
        if let Some(q) = &q2 {
            command.args(["-q2", q.to_str().unwrap()]);
        }
        if !command.status().map_err(|e| e.to_string())?.success() {
            return Err(format!("MainFilter failed for {}", s.id));
        }
        let result = audit(
            &rdir.join("filtered/all_1.fq"),
            q2.as_ref()
                .map(|_| rdir.join("filtered/all_2.fq"))
                .as_deref(),
            a.k,
            &index,
            atomic.len(),
        )?;
        for i in 0..atomic.len() {
            support[i].0 += result.specific[i];
            support[i].1 += result.ambiguous[i];
            if result.specific[i] > 0 {
                support[i].2 += 1;
            }
            let rpm = if effective == 0 {
                0.0
            } else {
                1e6 * result.specific[i] as f64 / effective as f64
            };
            writeln!(
                signal,
                "{}\t{}\t{}\t{}\t{}\t{}\t{rpm:.6}\t{}",
                s.id,
                s.taxon,
                atomic[i],
                effective,
                result.specific[i],
                result.ambiguous[i],
                state(effective, result.specific[i], result.ambiguous[i])
            )
            .unwrap();
        }
        for (edge, counts) in result.links {
            let item = all_links.entry(edge).or_default();
            item.0 += counts.0;
            item.1 += counts.1;
        }
        writeln!(
            samples_out,
            "{}\t{}\t{}\t{}",
            s.id,
            s.taxon,
            effective,
            if q2.is_some() { 1 } else { 0 }
        )
        .unwrap();
    }
    let meta = read_atomic_meta(a)?;
    let mut key_to_eq: BTreeMap<String, usize> = BTreeMap::new();
    let mut eq_members: Vec<Vec<usize>> = Vec::new();
    let mut eq_for = vec![0_usize; atomic.len()];
    for (i, id) in atomic.iter().enumerate() {
        let key = fasta_key(&discover.join("seeds").join(format!("{id}.fasta")))?;
        let eq = match key_to_eq.get(&key) {
            Some(&existing) => existing,
            None => {
                let next = eq_members.len();
                key_to_eq.insert(key, next);
                eq_members.push(Vec::new());
                next
            }
        };
        eq_members[eq].push(i);
        eq_for[i] = eq;
    }
    let mut map =
        BufWriter::new(File::create(temp.join("equivalence_map.tsv")).map_err(|e| e.to_string())?);
    let mut catalog =
        BufWriter::new(File::create(temp.join("curated_catalog.tsv")).map_err(|e| e.to_string())?);
    let mut topology =
        BufWriter::new(File::create(temp.join("topology.tsv")).map_err(|e| e.to_string())?);
    let mut evidence = BufWriter::new(
        File::create(temp.join("curation_evidence.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        map,
        "atomic_id\tequivalence_id\tmerge_evidence\tsequence_key_hash"
    )
    .unwrap();
    writeln!(
        catalog,
        "equivalence_id\tatomic_count\tclass\tcuration_status"
    )
    .unwrap();
    writeln!(topology, "equivalence_id\tclass\tevidence").unwrap();
    writeln!(evidence, "equivalence_id\tatomic_ids\tseed_count\tdiscovery_kmer_count\ttentative_annotation\tspecific_pairs_total\tambiguous_pairs_total\tspecific_fraction\tsamples_with_specific_support\tpe_link_degree\tread_transition_degree\tsequence_complexity\tconfidence\tdecision\treason").unwrap();
    for (eq_i, members) in eq_members.iter().enumerate() {
        let eq = format!("EQ{:05}", eq_i + 1);
        let representative = members[0];
        let id = &atomic[representative];
        fs::copy(
            discover.join("seeds").join(format!("{id}.fasta")),
            temp.join("library").join(format!("{eq}.fasta")),
        )
        .map_err(|e| e.to_string())?;
        let first = meta.get(id).cloned().unwrap_or_default();
        let key = fasta_key(&discover.join("seeds").join(format!("{id}.fasta")))?;
        let specific = members.iter().map(|&i| support[i].0).max().unwrap_or(0);
        let ambiguous = members.iter().map(|&i| support[i].1).max().unwrap_or(0);
        let sample_support = members.iter().map(|&i| support[i].2).max().unwrap_or(0);
        let complexity = entropy2(&discover.join("seeds").join(format!("{id}.fasta")))?;
        let mut bridge_degree = 0_u64;
        let mut transition_degree = 0_u64;
        for ((left, right), (bridges, transitions)) in &all_links {
            if eq_for[*left] == eq_i && eq_for[*right] != eq_i
                || eq_for[*right] == eq_i && eq_for[*left] != eq_i
            {
                bridge_degree += *bridges;
                transition_degree += *transitions;
            }
        }
        let fraction = specific as f64 / (specific + ambiguous).max(1) as f64;
        let (confidence, decision, reason) = if complexity < 1.2 {
            (
                "low_complexity",
                "flagged_low_complexity",
                "2mer_entropy_below_1.2",
            )
        } else if specific >= 10 && fraction >= 0.70 {
            ("high", "retained", "specific_support")
        } else if specific > 0 {
            (
                "provisional",
                "retained_low_support",
                "limited_specific_support",
            )
        } else if ambiguous > 0 {
            ("ambiguous", "flagged_ambiguous", "no_specific_support")
        } else {
            ("unresolved", "retained_low_support", "no_candidate_support")
        };
        let atomic_ids = members
            .iter()
            .map(|&i| atomic[i].clone())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            catalog,
            "{eq}\t{}\t{}\t{decision}",
            members.len(),
            first.class
        )
        .unwrap();
        writeln!(topology, "{eq}\t{}\tseed_motif_only", first.class).unwrap();
        writeln!(evidence, "{eq}\t{atomic_ids}\t{}\t{}\t{}\t{specific}\t{ambiguous}\t{fraction:.6}\t{sample_support}\t{bridge_degree}\t{transition_degree}\t{complexity:.6}\t{confidence}\t{decision}\t{reason}", first.seed_count, first.kmer_count, first.annotation).unwrap();
        for &i in members {
            let evidence = if members.len() > 1 {
                "exact_sequence_or_rc"
            } else {
                "identity_only"
            };
            writeln!(
                map,
                "{}\t{eq}\t{evidence}\t{:016x}",
                atomic[i],
                pair_hash(&key)
            )
            .unwrap();
        }
    }
    let mut linkage =
        BufWriter::new(File::create(temp.join("repeat_linkage.tsv")).map_err(|e| e.to_string())?);
    writeln!(
        linkage,
        "equivalence_a\tequivalence_b\tpe_bridges\tread_transitions\tsupport_score\trelation\tmerge_allowed"
    )
    .unwrap();
    let mut equivalence_links: BTreeMap<(usize, usize), (u64, u64)> = BTreeMap::new();
    for ((a_id, b_id), (bridges, transitions)) in all_links {
        let left = eq_for[a_id];
        let right = eq_for[b_id];
        if left != right {
            let item = equivalence_links
                .entry((left.min(right), left.max(right)))
                .or_default();
            item.0 += bridges;
            item.1 += transitions;
        }
    }
    for ((left, right), (bridges, transitions)) in equivalence_links {
        if bridges > 0 || transitions > 0 {
            let score = (bridges as f64 + 1.0).ln() + (transitions as f64 + 1.0).ln();
            writeln!(
                linkage,
                "EQ{:05}\tEQ{:05}\t{bridges}\t{transitions}\t{score:.6}\tlinked_not_merged\tfalse",
                left + 1,
                right + 1
            )
            .unwrap();
        }
    }
    drop(signal);
    drop(samples_out);
    drop(map);
    drop(catalog);
    drop(topology);
    drop(evidence);
    drop(linkage);
    let artifact_hash = curated_artifact_hash(&temp)?;
    let mut manifest =
        BufWriter::new(File::create(temp.join("manifest.tsv")).map_err(|e| e.to_string())?);
    writeln!(manifest, "schema_version\tstage\tdiscover_manifest_hash\tkmer\tsample_manifest_hash\tledger_hash\tartifact_hash").unwrap();
    writeln!(
        manifest,
        "1\tcurate\t{:016x}\t{}\t{:016x}\t{:016x}\t{artifact_hash:016x}",
        file_hash(&discover.join("manifest.tsv"))?,
        a.k,
        file_hash(&a.samples)?,
        a.ledger.as_deref().map(file_hash).transpose()?.unwrap_or(0),
    )
    .unwrap();
    commit(&temp, &final_d, "curate")
}
fn validate_curate(a: &A) -> R<()> {
    validate_discover(a)?;
    let root = a.out.join("02_curate");
    let fields = manifest_fields(&root.join("manifest.tsv"))
        .map_err(|_| "missing or malformed curated library; rerun --stage curate".to_string())?;
    let ledger_hash = a.ledger.as_deref().map(file_hash).transpose()?.unwrap_or(0);
    if manifest_value(&fields, "schema_version")? != "1"
        || manifest_value(&fields, "stage")? != "curate"
        || manifest_hash(&fields, "discover_manifest_hash")?
            != file_hash(&a.out.join("01_discover/manifest.tsv"))?
        || manifest_number::<usize>(&fields, "kmer")? != a.k
        || manifest_hash(&fields, "sample_manifest_hash")? != file_hash(&a.samples)?
        || manifest_hash(&fields, "ledger_hash")? != ledger_hash
        || manifest_hash(&fields, "artifact_hash")? != curated_artifact_hash(&root)?
    {
        return Err(
            "curated output does not match current inputs or artifacts; rerun --stage curate"
                .into(),
        );
    }
    Ok(())
}
fn curation_status(a: &A) -> R<AHashMap<String, (String, String)>> {
    let text = fs::read_to_string(a.out.join("02_curate/curation_evidence.tsv"))
        .map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| "malformed curation evidence".to_string())?;
    let positions: BTreeMap<_, _> = header
        .split('\t')
        .enumerate()
        .map(|(i, x)| (x, i))
        .collect();
    let id = *positions
        .get("equivalence_id")
        .ok_or_else(|| "malformed curation evidence".to_string())?;
    let confidence = *positions
        .get("confidence")
        .ok_or_else(|| "malformed curation evidence".to_string())?;
    let decision = *positions
        .get("decision")
        .ok_or_else(|| "malformed curation evidence".to_string())?;
    let mut out = AHashMap::new();
    for line in lines {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() > decision {
            out.insert(f[id].into(), (f[confidence].into(), f[decision].into()));
        }
    }
    Ok(out)
}

#[derive(Clone)]
struct LibraryRecord {
    id: String,
    class: String,
    seq: Vec<u8>,
}
struct LibraryIndex {
    records: Vec<LibraryRecord>,
    positions: AHashMap<u64, Vec<(usize, usize)>>,
}
fn fasta_records(path: &Path) -> R<Vec<(String, Vec<u8>)>> {
    let mut out = Vec::new();
    let mut header = String::new();
    let mut sequence = Vec::new();
    for line in BufReader::new(open(path)?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(h) = line.strip_prefix('>') {
            if !header.is_empty() && !sequence.is_empty() {
                out.push((header.clone(), sequence.clone()));
            }
            header = h.split_whitespace().next().unwrap_or("unknown").into();
            sequence.clear();
        } else {
            sequence.extend(
                line.trim()
                    .bytes()
                    .filter(|b| enc(*b).is_some())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if !header.is_empty() && !sequence.is_empty() {
        out.push((header, sequence));
    }
    Ok(out)
}
fn library_index(path: Option<&Path>, k: usize) -> R<Option<LibraryIndex>> {
    let Some(path) = path else { return Ok(None) };
    let mut records = Vec::new();
    for (header, seq) in fasta_records(path)? {
        let mut parts = header.splitn(2, '#');
        let id = parts.next().unwrap_or("unknown").into();
        let class = parts.next().unwrap_or("unknown").into();
        if seq.len() >= k {
            records.push(LibraryRecord { id, class, seq });
        }
    }
    let mut positions: AHashMap<u64, Vec<(usize, usize)>> = AHashMap::new();
    for (record_i, record) in records.iter().enumerate() {
        let mut pos = 0;
        kmers(&record.seq, k, |x| {
            let hits = positions.entry(x).or_default();
            if hits.len() < 64 {
                hits.push((record_i, pos));
            }
            pos += 1;
        });
    }
    Ok(Some(LibraryIndex { records, positions }))
}
fn oriented_reads(reads: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(reads.len() * 2);
    for read in reads {
        out.push(read.clone());
        let text = String::from_utf8_lossy(read);
        out.push(rc(&text).into_bytes());
    }
    out
}
fn local_adjacency(reads: &[Vec<u8>], k: usize) -> (Adjacency, Adjacency) {
    let mut right = AHashMap::new();
    let mut left = AHashMap::new();
    for read in oriented_reads(reads) {
        if read.len() < k {
            continue;
        }
        for i in 0..=read.len() - k {
            if let Some(base) = enc(read[i + k - 1]) {
                right.entry(read[i..i + k - 1].to_vec()).or_insert([0; 4])[base as usize] += 1;
            }
            if let Some(base) = enc(read[i]) {
                left.entry(read[i + 1..i + k].to_vec()).or_insert([0; 4])[base as usize] += 1;
            }
        }
    }
    (right, left)
}
fn extend_fragments(
    seed: &[u8],
    reads: &[Vec<u8>],
    k: usize,
    max_len: usize,
    min_count: u32,
    branch_ratio: f64,
    max_fragments: usize,
) -> Vec<Vec<u8>> {
    // Construct the local DBG adjacency once. Keep cumulative edge support so a
    // retained bubble path is selected by evidence, never lexical sequence order.
    let (right_edges, left_edges) = local_adjacency(reads, k);
    let mut paths = vec![(seed.to_vec(), 0_u64)];
    for right in [true, false] {
        let mut active = paths;
        for _ in 0..max_len.saturating_sub(seed.len()) {
            let mut next: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
            let mut progressed = false;
            for (fragment, score) in active {
                if fragment.len() < k - 1 || fragment.len() >= max_len {
                    next.entry(fragment)
                        .and_modify(|x| *x = (*x).max(score))
                        .or_insert(score);
                    continue;
                }
                let key = if right {
                    fragment[fragment.len() - (k - 1)..].to_vec()
                } else {
                    fragment[..k - 1].to_vec()
                };
                let counts = if right {
                    right_edges.get(&key)
                } else {
                    left_edges.get(&key)
                };
                let mut rank: Vec<_> = counts
                    .into_iter()
                    .flat_map(|x| x.iter().enumerate())
                    .filter(|(_, n)| **n >= min_count)
                    .collect();
                rank.sort_by_key(|(_, n)| Reverse(**n));
                if rank.is_empty() {
                    next.entry(fragment)
                        .and_modify(|x| *x = (*x).max(score))
                        .or_insert(score);
                    continue;
                }
                let keep =
                    if rank.len() > 1 && (*rank[0].1 as f64 / *rank[1].1 as f64) < branch_ratio {
                        rank.len().min(max_fragments)
                    } else {
                        1
                    };
                for (base, support) in rank.into_iter().take(keep) {
                    let mut child = fragment.clone();
                    if right {
                        child.push(b"ACGT"[base]);
                    } else {
                        child.insert(0, b"ACGT"[base]);
                    }
                    let child_score = score + *support as u64;
                    next.entry(child)
                        .and_modify(|x| *x = (*x).max(child_score))
                        .or_insert(child_score);
                }
                progressed = true;
            }
            active = next.into_iter().collect();
            active.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
            active.truncate(max_fragments);
            if !progressed {
                break;
            }
        }
        paths = active;
    }
    paths.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    paths.into_iter().map(|x| x.0).collect()
}

fn diagonal_similarity(query: &[u8], target: &[u8], delta: isize) -> (f64, f64) {
    let begin_q = if delta < 0 { (-delta) as usize } else { 0 };
    let begin_t = if delta > 0 { delta as usize } else { 0 };
    let overlap = query
        .len()
        .saturating_sub(begin_q)
        .min(target.len().saturating_sub(begin_t));
    if overlap == 0 {
        return (0.0, 0.0);
    }
    let same = (0..overlap)
        .filter(|&i| query[begin_q + i] == target[begin_t + i])
        .count();
    (
        same as f64 / overlap as f64,
        overlap as f64 / query.len() as f64,
    )
}
fn best_library_hit(
    fragment: &[u8],
    index: Option<&LibraryIndex>,
    k: usize,
) -> Option<(String, String, f64, f64, f64, f64)> {
    let index = index?;
    let mut offsets: AHashMap<(usize, isize), u32> = AHashMap::new();
    let mut pos = 0;
    kmers(fragment, k, |x| {
        if let Some(hits) = index.positions.get(&x) {
            for &(record, record_pos) in hits {
                *offsets
                    .entry((record, record_pos as isize - pos as isize))
                    .or_default() += 1;
            }
        }
        pos += 1;
    });
    let mut candidates: Vec<_> = offsets.into_iter().collect();
    candidates.sort_by_key(|(_, n)| Reverse(*n));
    let mut scored = Vec::new();
    for ((record_i, offset), _) in candidates.into_iter().take(16) {
        let record = &index.records[record_i];
        let (identity, coverage) = diagonal_similarity(fragment, &record.seq, offset);
        scored.push((
            record.id.clone(),
            record.class.clone(),
            identity,
            coverage,
            identity * coverage,
        ));
    }
    scored.sort_by(|left, right| right.4.partial_cmp(&left.4).unwrap());
    let best = scored.first()?.clone();
    let second = scored
        .iter()
        .skip(1)
        .find(|x| x.1 != best.1)
        .map(|x| x.4)
        .unwrap_or(0.0);
    Some((best.0, best.1, best.2, best.3, best.4, second))
}
const FAMILY_SCHEMA_VERSION: u32 = 2;
const FAMILY_K: usize = 15;
const FAMILY_WINDOW: usize = 10;
const FAMILY_MIN_JACCARD: f64 = 0.18;
const FAMILY_MIN_IDENTITY: f64 = 0.80;
const FAMILY_MIN_COVERAGE: f64 = 0.55;

fn minimizer_sketch(sequence: &[u8]) -> Vec<u64> {
    if sequence.len() < FAMILY_K {
        return Vec::new();
    }
    let mut keys = Vec::with_capacity(sequence.len() - FAMILY_K + 1);
    kmers(sequence, FAMILY_K, |key| keys.push(key));
    let width = FAMILY_WINDOW.min(keys.len());
    let mut sketch = AHashSet::new();
    for window in keys.windows(width) {
        if let Some(&key) = window.iter().min() {
            sketch.insert(key);
        }
    }
    let mut sketch: Vec<_> = sketch.into_iter().collect();
    sketch.sort_unstable();
    sketch
}

fn sketch_jaccard(left: &[u64], right: &[u64]) -> f64 {
    let (mut i, mut j, mut common) = (0, 0, 0);
    while i < left.len() && j < right.len() {
        match left[i].cmp(&right[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }
    let union = left.len() + right.len() - common;
    if union == 0 {
        0.0
    } else {
        common as f64 / union as f64
    }
}

fn local_family_alignment(left: &[u8], right: &[u8]) -> (f64, f64) {
    // Smith-Waterman is only run after minimizer candidate screening.
    let mut previous = vec![(0_i32, 0_u16, 0_u16); right.len() + 1];
    let mut best = (0_i32, 0_u16, 0_u16);
    for &a in left {
        let mut current = vec![(0_i32, 0_u16, 0_u16); right.len() + 1];
        for (j, &b) in right.iter().enumerate() {
            let diagonal = previous[j];
            let same = a == b;
            let mut choice = (
                diagonal.0 + if same { 2 } else { -2 },
                diagonal.1 + u16::from(same),
                diagonal.2 + 1,
            );
            for candidate in [
                (
                    previous[j + 1].0 - 3,
                    previous[j + 1].1,
                    previous[j + 1].2 + 1,
                ),
                (current[j].0 - 3, current[j].1, current[j].2 + 1),
                (0, 0, 0),
            ] {
                if candidate.0 > choice.0 || (candidate.0 == choice.0 && candidate.1 > choice.1) {
                    choice = candidate;
                }
            }
            current[j + 1] = choice;
            if choice.0 > best.0 || (choice.0 == best.0 && choice.1 > best.1) {
                best = choice;
            }
        }
        previous = current;
    }
    if best.2 == 0 {
        (0.0, 0.0)
    } else {
        (
            best.1 as f64 / best.2 as f64,
            best.2 as f64 / left.len().min(right.len()).max(1) as f64,
        )
    }
}

fn write_repeat_families(fragment_dir: &Path, groups: &[String], minimum_length: usize) -> R<()> {
    let mut sequences = Vec::with_capacity(groups.len());
    for group in groups {
        let path = fragment_dir.join(format!("{group}.fasta"));
        let mut candidates = if path.is_file() {
            fasta_records(&path)?
        } else {
            Vec::new()
        };
        candidates.sort_by_key(|candidate| std::cmp::Reverse(candidate.1.len()));
        sequences.push(
            candidates
                .into_iter()
                .next()
                .map(|x| x.1)
                .filter(|sequence| sequence.len() >= minimum_length)
                .unwrap_or_default(),
        );
    }
    let sketches: Vec<_> = sequences
        .iter()
        .map(|sequence| minimizer_sketch(sequence))
        .collect();
    let mut inverted: AHashMap<u64, Vec<usize>> = AHashMap::new();
    for (index, sketch) in sketches.iter().enumerate() {
        for &key in sketch {
            inverted.entry(key).or_default().push(index);
        }
    }
    let mut candidates: BTreeSet<(usize, usize)> = BTreeSet::new();
    for indexes in inverted.values() {
        if indexes.len() > 128 {
            continue;
        }
        for (offset, &left) in indexes.iter().enumerate() {
            for &right in &indexes[offset + 1..] {
                candidates.insert((left.min(right), left.max(right)));
            }
        }
    }
    let mut dsu = Dsu::new(groups.len());
    let mut edges = Vec::new();
    for (left, right) in candidates {
        let jaccard = sketch_jaccard(&sketches[left], &sketches[right]);
        if jaccard < FAMILY_MIN_JACCARD {
            continue;
        }
        let (identity, coverage) = local_family_alignment(&sequences[left], &sequences[right]);
        if identity >= FAMILY_MIN_IDENTITY && coverage >= FAMILY_MIN_COVERAGE {
            dsu.join(left, right);
            edges.push((left, right));
        }
    }
    let mut members: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for index in 0..groups.len() {
        let root = dsu.find(index);
        members.entry(root).or_default().push(index);
    }
    let mut ordered: Vec<_> = members.into_values().collect();
    ordered.sort_by_key(|members| members.iter().map(|&index| groups[index].clone()).min());
    let mut family_for = vec![String::new(); groups.len()];
    let mut size_for = vec![1_usize; groups.len()];
    for (family_index, members) in ordered.iter().enumerate() {
        let id = format!("FAM{:05}", family_index + 1);
        for &member in members {
            family_for[member] = id.clone();
            size_for[member] = members.len();
        }
    }
    let root = fragment_dir
        .parent()
        .ok_or("fragment directory has no parent")?;
    let mut output =
        BufWriter::new(File::create(root.join("repeat_families.tsv")).map_err(|e| e.to_string())?);
    writeln!(
        output,
        "family_id	equivalence_id	family_size	representative_length	accepted_family_edges"
    )
    .unwrap();
    let mut edge_count: AHashMap<usize, usize> = AHashMap::new();
    for (left, right) in edges {
        *edge_count.entry(left).or_insert(0) += 1;
        *edge_count.entry(right).or_insert(0) += 1;
    }
    for (index, group) in groups.iter().enumerate() {
        writeln!(
            output,
            "{}	{group}	{}	{}	{}",
            family_for[index],
            size_for[index],
            sequences[index].len(),
            edge_count.get(&index).copied().unwrap_or(0)
        )
        .unwrap();
    }
    Ok(())
}

fn annotation_artifact_hash(root: &Path) -> R<u64> {
    let mut state = 0xcbf29ce484222325_u64;
    for name in [
        "annotated_catalog.tsv",
        "annotation_evidence.tsv",
        "fragment_metrics.tsv",
        "repeat_families.tsv",
    ] {
        state = mix_hash(state, name.as_bytes());
        state = mix_hash(state, &file_hash(&root.join(name))?.to_le_bytes());
    }
    state = mix_hash(state, &tree_hash(&root.join("fragments"))?.to_le_bytes());
    Ok(state)
}
fn annotate(a: &A, ss: &[S]) -> R<()> {
    validate_curate(a)?;
    let curate = a.out.join("02_curate");
    let (groups, index) = seed_index_dir(a, &curate.join("library"))?;
    if groups.is_empty() {
        return Err("curated library has no equivalence groups".into());
    }
    let mut reads: Vec<Vec<Vec<u8>>> = vec![Vec::new(); groups.len()];
    let mut pairs = vec![0_u64; groups.len()];
    for sample in ss {
        let filtered = curate
            .join("candidate_recruit")
            .join(&sample.id)
            .join("filtered");
        let one = filtered.join("all_1.fq");
        if !one.exists() {
            continue;
        }
        let mut reader = PairReader::new(
            &one,
            sample
                .r2
                .as_ref()
                .map(|_| filtered.join("all_2.fq"))
                .as_deref(),
        )?;
        while let Some((left, right)) = reader.next_pair()? {
            let mut assigned = groups_for(&left.seq, a.k, &index);
            if let Some(right) = right.as_ref() {
                assigned.extend(groups_for(&right.seq, a.k, &index));
            }
            if assigned.len() != 1 {
                continue;
            }
            let group = *assigned.iter().next().unwrap();
            pairs[group] += 1;
            if reads[group].len() < 4096 {
                reads[group].push(left.seq);
                if let Some(right) = right {
                    reads[group].push(right.seq);
                }
            }
        }
    }
    let library = library_index(a.library.as_deref(), a.k)?;
    let final_d = a.out.join("03_annotate");
    let temp = a.out.join(format!(".annotate.tmp.{}", std::process::id()));
    if temp.exists() {
        fs::remove_dir_all(&temp).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(temp.join("fragments")).map_err(|e| e.to_string())?;
    let mut catalog = BufWriter::new(
        File::create(temp.join("annotated_catalog.tsv")).map_err(|e| e.to_string())?,
    );
    let mut evidence = BufWriter::new(
        File::create(temp.join("annotation_evidence.tsv")).map_err(|e| e.to_string())?,
    );
    let mut metrics =
        BufWriter::new(File::create(temp.join("fragment_metrics.tsv")).map_err(|e| e.to_string())?);
    writeln!(
        metrics,
        "equivalence_id\tfragment_id\tlength_bp\tunique_pairs\tassembly_state"
    )
    .unwrap();
    writeln!(
        catalog,
        "equivalence_id\tclass\tannotation_confidence\tdecision"
    )
    .unwrap();
    writeln!(evidence, "equivalence_id\trepresentative_fragment_id\tfragment_length\tunique_pairs\tperiod_bp\tinverted_repeat_identity\tbest_library_id\tbest_library_class\talignment_identity\talignment_coverage\tscore_delta\tfinal_class\tannotation_confidence\tdecision\treason").unwrap();
    for (i, group) in groups.iter().enumerate() {
        let seed_path = curate.join("library").join(format!("{group}.fasta"));
        let seed = fasta_records(&seed_path)?
            .into_iter()
            .next()
            .map(|x| x.1)
            .unwrap_or_default();
        let fragments = if seed.len() >= a.k {
            extend_fragments(
                &seed,
                &reads[i],
                a.k,
                a.annotation_max_fragment,
                a.assemble_min_kmer_count,
                a.assemble_branch_ratio,
                a.assemble_max_fragments,
            )
        } else {
            vec![seed.clone()]
        };
        let (fragment_id, fragment) = fragments
            .iter()
            .enumerate()
            .max_by_key(|(_, x)| x.len())
            .map(|(n, x)| (n + 1, x.clone()))
            .unwrap_or((0, Vec::new()));
        let periodic = fragments
            .iter()
            .enumerate()
            .filter_map(|(n, x)| period(&String::from_utf8_lossy(x)).map(|p| (n + 1, p)))
            .next();
        let period_bp = periodic.map(|x| x.1);
        let inverted_identity = inverted_repeat_identity(&fragment);
        let hit = fragments
            .iter()
            .enumerate()
            .filter(|(_, x)| x.len() >= a.annotation_min_fragment)
            .filter_map(|(n, x)| best_library_hit(x, library.as_ref(), a.k).map(|hit| (n + 1, hit)))
            .max_by(|left, right| left.1 .4.partial_cmp(&right.1 .4).unwrap());
        let evidence_fragment_id = periodic
            .map(|x| x.0)
            .or_else(|| hit.as_ref().map(|x| x.0))
            .unwrap_or(fragment_id);
        let evidence_fragment_len = fragments
            .get(evidence_fragment_id.saturating_sub(1))
            .map(|x| x.len())
            .unwrap_or(0);
        let (
            class,
            confidence,
            decision,
            reason,
            hit_id,
            hit_class,
            identity,
            coverage,
            score_delta,
        ): (
            String,
            String,
            String,
            String,
            String,
            String,
            f64,
            f64,
            f64,
        ) = if let Some(period_bp) = period_bp {
            let class = if period_bp <= 20 {
                "simple_repeat"
            } else {
                "satellite_candidate"
            };
            (
                class.into(),
                "high".into(),
                "retained".into(),
                "periodic_fragment".into(),
                ".".into(),
                ".".into(),
                0.0,
                0.0,
                0.0,
            )
        } else if inverted_identity >= 0.85 && pairs[i] as usize >= a.annotation_min_support {
            (
                "foldback_like_DNA".into(),
                "provisional".into(),
                "retained".into(),
                "long_inverted_repeat_no_long_orf".into(),
                ".".into(),
                ".".into(),
                inverted_identity,
                1.0,
                0.0,
            )
        } else if let Some((_, (id, class, identity, coverage, score, second_score))) = hit {
            if pairs[i] as usize >= a.annotation_min_support
                && identity >= a.annotation_min_identity
                && coverage >= a.annotation_min_coverage
                && score - second_score >= a.annotation_min_delta
            {
                (
                    class.clone(),
                    "high".into(),
                    "retained".into(),
                    "library_homology".into(),
                    id,
                    class,
                    identity,
                    coverage,
                    score - second_score,
                )
            } else {
                (
                    "unknown_interspersed_repeat".into(),
                    "low".into(),
                    "retained_flagged".into(),
                    "weak_library_homology".into(),
                    id,
                    class,
                    identity,
                    coverage,
                    score - second_score,
                )
            }
        } else if fragment.len() >= a.annotation_min_fragment
            && pairs[i] as usize >= a.annotation_min_support
        {
            (
                "unknown_interspersed_repeat".into(),
                "provisional".into(),
                "retained".into(),
                "no_confident_library_match".into(),
                ".".into(),
                ".".into(),
                0.0,
                0.0,
                0.0,
            )
        } else {
            (
                "unknown_repeat".into(),
                "low".into(),
                "retained_low_support".into(),
                "fragment_or_support_below_threshold".into(),
                ".".into(),
                ".".into(),
                0.0,
                0.0,
                0.0,
            )
        };
        if !fragment.is_empty() {
            let mut file = BufWriter::new(
                File::create(temp.join("fragments").join(format!("{group}.fasta")))
                    .map_err(|e| e.to_string())?,
            );
            for (fragment_i, sequence) in fragments.iter().enumerate() {
                writeln!(
                    file,
                    ">{group}|fragment_{}\n{}",
                    fragment_i + 1,
                    String::from_utf8_lossy(sequence)
                )
                .unwrap();
            }
        }
        for (fragment_i, sequence) in fragments.iter().enumerate() {
            let assembly_state = if sequence.len() <= seed.len() {
                "seed_only"
            } else {
                "assembled"
            };
            writeln!(
                metrics,
                "{group}\tfragment_{}\t{}\t{}\t{assembly_state}",
                fragment_i + 1,
                sequence.len(),
                pairs[i]
            )
            .unwrap();
        }
        writeln!(catalog, "{group}\t{class}\t{confidence}\t{decision}").unwrap();
        writeln!(evidence, "{group}\t{evidence_fragment_id}\t{}\t{}\t{}\t{inverted_identity:.6}\t{hit_id}\t{hit_class}\t{identity:.6}\t{coverage:.6}\t{score_delta:.6}\t{class}\t{confidence}\t{decision}\t{reason}", evidence_fragment_len, pairs[i], period_bp.map(|x| x.to_string()).unwrap_or_else(|| ".".into())).unwrap();
    }
    drop(catalog);
    drop(evidence);
    drop(metrics);
    write_repeat_families(&temp.join("fragments"), &groups, a.annotation_min_fragment)?;
    let artifact_hash = annotation_artifact_hash(&temp)?;
    let mut manifest =
        BufWriter::new(File::create(temp.join("manifest.tsv")).map_err(|e| e.to_string())?);
    writeln!(manifest, "schema_version\tstage\tcurate_manifest_hash\tte_library_hash\tmin_fragment\tmax_fragment\tmin_support\tmin_identity\tmin_coverage\tmin_delta\tassemble_min_kmer_count\tassemble_branch_ratio\tassemble_max_fragments\tfamily_schema_version\tartifact_hash").unwrap();
    writeln!(
        manifest,
        "2\tannotate\t{:016x}\t{:016x}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{}\t{artifact_hash:016x}",
        file_hash(&curate.join("manifest.tsv"))?,
        a.library
            .as_deref()
            .map(file_hash)
            .transpose()?
            .unwrap_or(0),
        a.annotation_min_fragment,
        a.annotation_max_fragment,
        a.annotation_min_support,
        a.annotation_min_identity,
        a.annotation_min_coverage,
        a.annotation_min_delta,
        a.assemble_min_kmer_count,
        a.assemble_branch_ratio,
        a.assemble_max_fragments,
        FAMILY_SCHEMA_VERSION
    )
    .unwrap();
    commit(&temp, &final_d, "annotate")
}
fn validate_annotate(a: &A) -> R<()> {
    validate_curate(a)?;
    let root = a.out.join("03_annotate");
    let fields = manifest_fields(&root.join("manifest.tsv")).map_err(|_| {
        "missing or malformed annotation output; rerun --stage annotate".to_string()
    })?;
    if manifest_value(&fields, "schema_version")? != "2"
        || manifest_value(&fields, "stage")? != "annotate"
        || manifest_hash(&fields, "curate_manifest_hash")?
            != file_hash(&a.out.join("02_curate/manifest.tsv"))?
        || manifest_hash(&fields, "te_library_hash")?
            != a.library
                .as_deref()
                .map(file_hash)
                .transpose()?
                .unwrap_or(0)
        || manifest_number::<usize>(&fields, "min_fragment")? != a.annotation_min_fragment
        || manifest_number::<usize>(&fields, "max_fragment")? != a.annotation_max_fragment
        || manifest_number::<usize>(&fields, "min_support")? != a.annotation_min_support
        || manifest_value(&fields, "min_identity")? != format!("{:.6}", a.annotation_min_identity)
        || manifest_value(&fields, "min_coverage")? != format!("{:.6}", a.annotation_min_coverage)
        || manifest_value(&fields, "min_delta")? != format!("{:.6}", a.annotation_min_delta)
        || manifest_number::<u32>(&fields, "assemble_min_kmer_count")? != a.assemble_min_kmer_count
        || manifest_value(&fields, "assemble_branch_ratio")?
            != format!("{:.6}", a.assemble_branch_ratio)
        || manifest_number::<usize>(&fields, "assemble_max_fragments")? != a.assemble_max_fragments
        || manifest_number::<u32>(&fields, "family_schema_version")? != FAMILY_SCHEMA_VERSION
        || manifest_hash(&fields, "artifact_hash")? != annotation_artifact_hash(&root)?
    {
        return Err(
            "annotation output does not match current inputs or artifacts; rerun --stage annotate"
                .into(),
        );
    }
    Ok(())
}
fn annotation_status(a: &A) -> R<AHashMap<String, (String, String, String)>> {
    let text = fs::read_to_string(a.out.join("03_annotate/annotated_catalog.tsv"))
        .map_err(|e| e.to_string())?;
    let mut out = AHashMap::new();
    for (n, line) in text.lines().enumerate() {
        if n == 0 {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() == 4 {
            out.insert(f[0].into(), (f[1].into(), f[2].into(), f[3].into()));
        }
    }
    Ok(out)
}
fn family_status(a: &A) -> R<AHashMap<String, String>> {
    let text = fs::read_to_string(a.out.join("03_annotate/repeat_families.tsv"))
        .map_err(|e| e.to_string())?;
    let mut out = AHashMap::new();
    for (line_number, line) in text.lines().enumerate() {
        if line_number == 0 {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() >= 2 {
            out.insert(fields[1].into(), fields[0].into());
        }
    }
    Ok(out)
}
fn quantify(a: &A, ss: &[S]) -> R<()> {
    validate_annotate(a)?;
    let annotate = a.out.join("03_annotate");
    let (groups, index) = seed_index_dir(a, &annotate.join("fragments"))?;
    if groups.is_empty() {
        return Err("annotation has no fragments; rerun --stage annotate".into());
    }
    let lengths = fragment_lengths(&annotate.join("fragments"))?;
    let fragment_kmers = fragment_kmer_counts(&annotate.join("fragments"), a.k, &index, &groups)?;
    let fragment_map = fragment_map_index(&annotate.join("fragments"), &groups, a.k)?;
    let status = curation_status(a)?;
    let annotation = annotation_status(a)?;
    let families = family_status(a)?;
    let denied = ledger(a.ledger.as_deref())?;
    let final_d = a.out.join("04_quantify");
    let temp = a.out.join(format!(".quantify.tmp.{}", std::process::id()));
    if temp.exists() {
        fs::remove_dir_all(&temp).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&temp).map_err(|e| e.to_string())?;
    let mut out =
        BufWriter::new(File::create(temp.join("repeat_signal.tsv")).map_err(|e| e.to_string())?);
    writeln!(out, "sample_id\ttaxon_id\tequivalence_id\tfamily_id\teffective_pairs\tspecific_pairs\tambiguous_pairs\tsignal_rpm\tsignal_rpm_ci_low\tsignal_rpm_ci_high\tquantification_mode\testimated_genome_fraction\tmapped_bases\tmean_depth\tkmer_breadth\tanchor_identity\tstate\tconfidence\tdecision\tclass\tannotation_confidence\tannotation_decision").unwrap();
    let mut coverage = BufWriter::new(
        File::create(temp.join("fragment_coverage.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(coverage, "sample_id\tequivalence_id\tfamily_id\teq_representative_length\teq_mapped_bases\teq_mean_depth\teq_kmer_breadth\tanchor_identity\tspecific_pairs\tambiguous_pairs").unwrap();
    let mut landscape =
        BufWriter::new(File::create(temp.join("repeat_landscape.tsv")).map_err(|e| e.to_string())?);
    writeln!(landscape, "sample_id\ttaxon_id\tequivalence_id\tfamily_id\tdivergence_percent_bin\tanchored_reads\tmetric").unwrap();
    let mut matrix: BTreeMap<(String, usize), Vec<(f64, String)>> = BTreeMap::new();
    let selected: Vec<_> = ss
        .par_iter()
        .map(|sample| {
            selected_pair_hashes(
                sample,
                denied.get(&sample.id).unwrap_or(&AHashSet::new()),
                a.quantify_pairs,
            )
        })
        .collect::<R<Vec<_>>>()?;
    let audits: Vec<(Audit, u64)> = ss
        .par_iter()
        .zip(selected.par_iter())
        .map(|(sample, chosen)| {
            audit_sample(
                sample,
                denied.get(&sample.id).unwrap_or(&AHashSet::new()),
                a.k,
                &index,
                &fragment_map,
                groups.len(),
                chosen.as_ref(),
            )
        })
        .collect::<R<Vec<_>>>()?;
    for (sample, (result, effective)) in ss.iter().zip(audits) {
        for (i, group) in groups.iter().enumerate() {
            let rpm = if effective == 0 {
                0.0
            } else {
                1e6 * result.specific[i] as f64 / effective as f64
            };
            let length = *lengths.get(group).unwrap_or(&0);
            let depth = if length == 0 {
                0.0
            } else {
                result.mapped_bases[i] as f64 / length as f64
            };
            let breadth = result.covered_kmers[i].len() as f64
                / *fragment_kmers.get(group).unwrap_or(&1).max(&1) as f64;
            let call = state(effective, result.specific[i], result.ambiguous[i]);
            let anchor_identity = if result.identity_n[i] > 0 {
                format!(
                    "{:.6}",
                    result.identity_sum[i] / result.identity_n[i] as f64
                )
            } else {
                "NA".into()
            };
            let (confidence, decision) = status
                .get(group)
                .cloned()
                .unwrap_or_else(|| ("unresolved".into(), "retained_low_support".into()));
            let (class, annotation_confidence, annotation_decision) =
                annotation.get(group).cloned().unwrap_or_else(|| {
                    (
                        "unknown_repeat".into(),
                        "unresolved".into(),
                        "missing_annotation".into(),
                    )
                });
            let family = families
                .get(group)
                .map(String::as_str)
                .unwrap_or("UNASSIGNED");
            let (rpm_low, rpm_high) = bootstrap_rpm_ci(
                result.specific[i],
                effective,
                a.bootstrap_replicates,
                mix_hash(pair_hash(&sample.id), group.as_bytes()),
            );
            let mode = if a.quantify_pairs == 0 {
                "all_reads"
            } else {
                "deterministic_subsample"
            };
            let fraction = if a.estimate_genome_fraction && effective > 0 {
                format!("{:.8}", result.specific[i] as f64 / effective as f64)
            } else {
                "NA".into()
            };
            writeln!(out, "{}\t{}\t{}\t{family}\t{}\t{}\t{}\t{rpm:.6}\t{rpm_low:.6}\t{rpm_high:.6}\t{mode}\t{fraction}\t{}\t{depth:.6}\t{breadth:.6}\t{anchor_identity}\t{call}\t{confidence}\t{decision}\t{class}\t{annotation_confidence}\t{annotation_decision}", sample.id, sample.taxon, group, effective, result.specific[i], result.ambiguous[i], result.mapped_bases[i]).unwrap();
            writeln!(
                coverage,
                "{}\t{}\t{family}\t{}\t{}\t{depth:.6}\t{breadth:.6}\t{anchor_identity}\t{}\t{}",
                sample.id,
                group,
                length,
                result.mapped_bases[i],
                result.specific[i],
                result.ambiguous[i]
            )
            .unwrap();
            for (bin, count) in result.identity_bins[i].iter().enumerate() {
                if *count > 0 {
                    writeln!(
                        landscape,
                        "{}\t{}\t{}\t{family}\t{bin}\t{count}\tcopy_divergence_proxy",
                        sample.id, sample.taxon, group
                    )
                    .unwrap();
                }
            }
            matrix
                .entry((sample.taxon.clone(), i))
                .or_default()
                .push((rpm, call.into()));
        }
    }
    let mut summary = BufWriter::new(
        File::create(temp.join("taxon_repeat_matrix.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        summary,
        "taxon_id\tequivalence_id\tfamily_id\tmedian_rpm\tstate"
    )
    .unwrap();
    for ((taxon, i), mut values) in matrix {
        values.retain(|x| x.1 != "NA");
        values.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let callable = values.len();
        let rpm = if callable == 0 {
            0.0
        } else {
            values[callable / 2].0
        };
        let present = values.iter().filter(|x| x.1 == "PRESENT").count();
        let call = if callable == 0 {
            "NA"
        } else if present * 2 >= callable {
            "PRESENT"
        } else {
            "ABSENT"
        };
        let family = families
            .get(&groups[i])
            .map(String::as_str)
            .unwrap_or("UNASSIGNED");
        writeln!(
            summary,
            "{taxon}\t{}\t{family}\t{rpm:.6}\t{call}",
            groups[i]
        )
        .unwrap();
    }
    commit(&temp, &final_d, "quantify")
}
struct ClusterRead {
    seq: Vec<u8>,
    mate: Option<Vec<u8>>,
}
fn interspersed(a: &A, ss: &[S]) -> R<()> {
    validate_curate(a)?;
    const LIMIT: usize = 50_000;
    const MIN_READS: usize = 20;
    let mut heap: BinaryHeap<Ranked> = BinaryHeap::new();
    for sample in ss {
        let root = a
            .out
            .join("02_curate/candidate_recruit")
            .join(&sample.id)
            .join("filtered");
        let one = root.join("all_1.fq");
        if !one.exists() {
            continue;
        }
        let mut reader = PairReader::new(
            &one,
            sample.r2.as_ref().map(|_| root.join("all_2.fq")).as_deref(),
        )?;
        while let Some((left, right)) = reader.next_pair()? {
            let rank = pair_hash(&left.id);
            let item = Ranked {
                rank,
                r1: left.seq,
                r2: right.map(|x| x.seq),
            };
            if heap.len() < LIMIT {
                heap.push(item);
            } else if rank < heap.peek().unwrap().rank {
                heap.pop();
                heap.push(item);
            }
        }
    }
    let reads: Vec<ClusterRead> = heap
        .into_iter()
        .map(|x| ClusterRead {
            seq: x.r1,
            mate: x.r2,
        })
        .collect();
    let sketches: Vec<Vec<u64>> = reads.iter().map(|r| minimizer_sketch(&r.seq)).collect();
    let mut inverted: AHashMap<u64, Vec<usize>> = AHashMap::new();
    for (i, sketch) in sketches.iter().enumerate() {
        for &key in sketch {
            let v = inverted.entry(key).or_default();
            if v.len() < 64 {
                v.push(i);
            }
        }
    }
    let mut dsu = Dsu::new(reads.len());
    for members in inverted.values() {
        for (offset, &left) in members.iter().enumerate() {
            for &right in members[offset + 1..].iter().take(16) {
                if sketch_jaccard(&sketches[left], &sketches[right]) >= 0.18 {
                    dsu.join(left, right);
                }
            }
        }
    }
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..reads.len() {
        let root = dsu.find(i);
        components.entry(root).or_default().push(i);
    }
    let mut components: Vec<_> = components
        .into_values()
        .filter(|v| v.len() >= MIN_READS)
        .collect();
    components.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a[0].cmp(&b[0])));
    components.truncate(256);
    let final_d = a.out.join("03_interspersed");
    let temp = a
        .out
        .join(format!(".interspersed.tmp.{}", std::process::id()));
    if temp.exists() {
        fs::remove_dir_all(&temp).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&temp).map_err(|e| e.to_string())?;
    let mut table =
        BufWriter::new(File::create(temp.join("clusters.tsv")).map_err(|e| e.to_string())?);
    let mut fasta =
        BufWriter::new(File::create(temp.join("consensus.fasta")).map_err(|e| e.to_string())?);
    writeln!(
        table,
        "cluster_id\tread_pairs\tconsensus_length\tstructure\tinverted_repeat_identity\tperiod_bp"
    )
    .unwrap();
    for (n, members) in components.iter().enumerate() {
        let mut seqs = Vec::new();
        for &i in members {
            seqs.push(reads[i].seq.clone());
            if let Some(mate) = &reads[i].mate {
                seqs.push(mate.clone());
            }
        }
        let seed = seqs
            .iter()
            .max_by_key(|x| x.len())
            .cloned()
            .unwrap_or_default();
        let paths = extend_fragments(
            &seed,
            &seqs,
            a.k,
            a.annotation_max_fragment,
            a.assemble_min_kmer_count,
            a.assemble_branch_ratio,
            a.assemble_max_fragments,
        );
        let consensus = paths.into_iter().max_by_key(|x| x.len()).unwrap_or(seed);
        let inv = inverted_repeat_identity(&consensus);
        let per = period(&String::from_utf8_lossy(&consensus));
        let structure = if per.is_some() {
            "tandem_repeat_candidate"
        } else if inv >= 0.80 {
            "foldback_like_DNA_candidate"
        } else {
            "interspersed_repeat_candidate"
        };
        let id = format!("IC{:05}", n + 1);
        writeln!(
            fasta,
            ">{id} reads={} structure={structure}\n{}",
            members.len(),
            String::from_utf8_lossy(&consensus)
        )
        .unwrap();
        writeln!(
            table,
            "{id}\t{}\t{}\t{structure}\t{inv:.6}\t{}",
            members.len(),
            consensus.len(),
            per.map(|x| x.to_string()).unwrap_or_else(|| ".".into())
        )
        .unwrap();
    }
    commit(&temp, &final_d, "interspersed")
}

type ComparisonRow = (
    BTreeSet<String>,
    BTreeSet<String>,
    BTreeSet<String>,
    BTreeSet<String>,
);
fn compare(a: &A) -> R<()> {
    let text = fs::read_to_string(a.out.join("04_quantify/repeat_signal.tsv"))
        .map_err(|_| "quantification output is required for compare".to_string())?;
    let mut lines = text.lines();
    let header = lines.next().ok_or("empty repeat_signal.tsv")?;
    let columns: BTreeMap<_, _> = header
        .split("\t")
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();
    let required = |name| {
        columns
            .get(name)
            .copied()
            .ok_or_else(|| format!("repeat_signal.tsv lacks {name}"))
    };
    let sample = required("sample_id")?;
    let taxon = required("taxon_id")?;
    let eq = required("equivalence_id")?;
    let family = required("family_id")?;
    let state_col = required("state")?;
    let class = required("class")?;
    let mut rows: BTreeMap<String, ComparisonRow> = BTreeMap::new();
    for line in lines {
        let f: Vec<_> = line.split("\t").collect();
        if f.len() <= state_col || f[state_col] != "PRESENT" {
            continue;
        }
        let entry = rows.entry(f[family].into()).or_default();
        entry.0.insert(f[sample].into());
        entry.1.insert(f[taxon].into());
        entry.2.insert(f[eq].into());
        entry.3.insert(f[class].into());
    }
    let final_d = a.out.join("05_compare");
    let temp = a.out.join(format!(".compare.tmp.{}", std::process::id()));
    if temp.exists() {
        fs::remove_dir_all(&temp).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&temp).map_err(|e| e.to_string())?;
    let mut out = BufWriter::new(
        File::create(temp.join("repeat_superfamilies.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(out, "superfamily_id\tsample_count\ttaxon_count\tshared_state\tequivalence_members\tclasses\tevidence") .unwrap();
    for (id, (samples, taxa, eqs, classes)) in rows {
        let state = if taxa.len() > 1 {
            "shared"
        } else if samples.len() > 1 {
            "taxon_shared"
        } else {
            "sample_specific"
        };
        writeln!(
            out,
            "{id}\t{}\t{}\t{state}\t{}\t{}\tfamily_consensus_and_read_support",
            samples.len(),
            taxa.len(),
            eqs.into_iter().collect::<Vec<_>>().join(","),
            classes.into_iter().collect::<Vec<_>>().join(",")
        )
        .unwrap();
    }
    commit(&temp, &final_d, "compare")
}
fn main() -> std::process::ExitCode {
    let a = match args() {
        Ok(x) => x,
        Err(e) => {
            eprintln!("{e}");
            return std::process::ExitCode::from(2);
        }
    };
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build_global();
    let ss = match samples(&a.samples) {
        Ok(x) => x,
        Err(e) => {
            eprintln!("{e}");
            return std::process::ExitCode::from(2);
        }
    };
    let led = match ledger(a.ledger.as_deref()) {
        Ok(x) => x,
        Err(e) => {
            eprintln!("{e}");
            return std::process::ExitCode::from(2);
        }
    };
    if a.stage == "all" || a.stage == "discover" {
        if let Err(e) = discover(&a, &ss, &led) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    if a.stage == "all" || a.stage == "curate" {
        if let Err(e) = curate(&a, &ss, &led) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    if a.stage == "all" || a.stage == "annotate" {
        if let Err(e) = annotate(&a, &ss) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    if a.stage == "all" || a.stage == "interspersed" {
        if let Err(e) = interspersed(&a, &ss) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    if a.stage == "all" || a.stage == "quantify" {
        if let Err(e) = quantify(&a, &ss) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    if a.stage == "compare" {
        if let Err(e) = compare(&a) {
            eprintln!("{e}");
            return std::process::ExitCode::from(1);
        }
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn canonical_motif_merges_reverse_complements_and_rotations() {
        assert_eq!(canonical_motif("AAG"), canonical_motif("CTT"));
        assert_eq!(canonical_motif("AAG"), canonical_motif("AGA"));
    }

    #[test]
    fn foldback_identity_detects_a_long_palindrome() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(inverted_repeat_identity(sequence) >= 0.99);
    }

    #[test]
    fn pair_hash_is_stable() {
        assert_eq!(pair_hash("read-1"), pair_hash("read-1"));
        assert_ne!(pair_hash("read-1"), pair_hash("read-2"));
    }

    #[test]
    fn fragment_index_caps_positions_per_fragment() {
        let root =
            std::env::temp_dir().join(format!("main_repeat_fragments_{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        fs::write(
            root.join("EQ00001.fasta"),
            ">a\nAAAAAAAAAAAA\n>b\nAAAAAAAAAAAA\n",
        )
        .unwrap();
        let index = fragment_map_index(&root, &["EQ00001".into()], 4).unwrap();
        fs::remove_dir_all(root).unwrap();
        let key = {
            let mut out = 0;
            kmers(b"AAAA", 4, |x| out = x);
            out
        };
        let hits = index.positions.get(&key).unwrap();
        assert!(hits.iter().any(|x| x.0 == 0));
        assert!(hits.iter().any(|x| x.0 == 1));
    }

    #[test]
    fn anchored_mapping_accepts_reverse_complement_reads() {
        let seq = b"AAGTCCTC".to_vec();
        let mut positions = AHashMap::new();
        let mut pos = 0;
        kmers(&seq, 4, |x| {
            positions.entry(x).or_insert_with(Vec::new).push((0, pos));
            pos += 1;
        });
        let map = FragmentMapIndex {
            records: vec![FragmentRecord { group: 0, seq }],
            positions,
        };
        let hit = anchored_alignment(b"AAGTCCTC", 0, &map, 4).unwrap();
        let reverse = anchored_alignment(b"GAGGACTT", 0, &map, 4).unwrap();
        assert_eq!(hit.0, 8);
        assert_eq!(reverse.0, 8);
        assert!(hit.1 > 0.99 && reverse.1 > 0.99);
    }

    #[test]
    fn bounded_local_assembly_retains_supported_bubble_paths() {
        let reads = vec![
            b"AACCGGTTA".to_vec(),
            b"AACCGGTTA".to_vec(),
            b"AACCGGTTC".to_vec(),
            b"AACCGGTTC".to_vec(),
        ];
        let paths = extend_fragments(b"AACCGGTT", &reads, 5, 12, 2, 1.5, 3);
        assert!(paths.len() >= 2);
        assert!(paths.iter().all(|x| x.len() <= 12));
    }

    #[test]
    fn family_screen_requires_shared_minimizers_and_local_identity() {
        let left = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let right = b"ACGTACGTACGTACGTACGTACGTACGTTCGT";
        let unrelated = b"TTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGG";
        assert!(
            sketch_jaccard(&minimizer_sketch(left), &minimizer_sketch(right)) >= FAMILY_MIN_JACCARD
        );
        let (identity, coverage) = local_family_alignment(left, right);
        assert!(identity >= FAMILY_MIN_IDENTITY && coverage >= FAMILY_MIN_COVERAGE);
        assert!(
            sketch_jaccard(&minimizer_sketch(left), &minimizer_sketch(unrelated))
                < FAMILY_MIN_JACCARD
        );
    }

    #[test]
    fn family_output_keeps_eqs_but_groups_similar_fragments() {
        let root =
            std::env::temp_dir().join(format!("main_repeat_families_{}", std::process::id()));
        let fragments = root.join("fragments");
        fs::create_dir_all(&fragments).unwrap();
        fs::write(
            fragments.join("EQ00001.fasta"),
            ">a\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )
        .unwrap();
        fs::write(
            fragments.join("EQ00002.fasta"),
            ">b\nACGTACGTACGTACGTACGTACGTACGTTCGT\n",
        )
        .unwrap();
        write_repeat_families(&fragments, &["EQ00001".into(), "EQ00002".into()], 16).unwrap();
        let rows = fs::read_to_string(root.join("repeat_families.tsv")).unwrap();
        fs::remove_dir_all(root).unwrap();
        let mut fields = rows
            .lines()
            .skip(1)
            .map(|row| row.split('\t').next().unwrap())
            .collect::<Vec<_>>();
        fields.sort_unstable();
        assert_eq!(fields, vec!["FAM00001", "FAM00001"]);

        let root =
            std::env::temp_dir().join(format!("main_repeat_short_families_{}", std::process::id()));
        let fragments = root.join("fragments");
        fs::create_dir_all(&fragments).unwrap();
        fs::write(
            fragments.join("EQ00001.fasta"),
            ">a\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )
        .unwrap();
        fs::write(
            fragments.join("EQ00002.fasta"),
            ">b\nACGTACGTACGTACGTACGTACGTACGTTCGT\n",
        )
        .unwrap();
        write_repeat_families(&fragments, &["EQ00001".into(), "EQ00002".into()], 80).unwrap();
        let rows = fs::read_to_string(root.join("repeat_families.tsv")).unwrap();
        fs::remove_dir_all(root).unwrap();
        let mut fields = rows
            .lines()
            .skip(1)
            .map(|row| row.split('\t').next().unwrap())
            .collect::<Vec<_>>();
        fields.sort_unstable();
        assert_eq!(fields, vec!["FAM00001", "FAM00002"]);
    }

    #[test]
    fn deterministic_subsample_uses_stable_pair_hashes() {
        let root =
            std::env::temp_dir().join(format!("main_repeat_subsample_{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let reads = root.join("reads.fq");
        fs::write(&reads, "\x40r3\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n\x40r1\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n\x40r2\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n").unwrap();
        let sample = S {
            taxon: "t".into(),
            id: "s".into(),
            r1: reads,
            r2: None,
        };
        let first = selected_pair_hashes(&sample, &AHashSet::new(), 2)
            .unwrap()
            .unwrap();
        let second = selected_pair_hashes(&sample, &AHashSet::new(), 2)
            .unwrap()
            .unwrap();
        fs::remove_dir_all(root).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.len(), 2);
    }

    #[test]
    fn bootstrap_interval_is_deterministic_and_bounded() {
        let first = bootstrap_rpm_ci(12, 100, 200, 9);
        let second = bootstrap_rpm_ci(12, 100, 200, 9);
        assert_eq!(first, second);
        assert!(first.0 <= 120_000.0 && first.1 >= 120_000.0);
    }

    #[test]
    fn library_hit_requires_a_distinguishable_class() {
        let path = std::env::temp_dir().join(format!("main_repeat_library_{}", std::process::id()));
        fs::write(
            &path,
            ">a#DNA/TcMar
ACGTACGTACGTACGT
>b#LTR/Gypsy
TTTTGGGGTTTTGGGG
",
        )
        .unwrap();
        let index = library_index(Some(&path), 4).unwrap().unwrap();
        fs::remove_file(path).unwrap();
        let hit = best_library_hit(b"ACGTACGTACGT", Some(&index), 4).unwrap();
        assert_eq!(hit.1, "DNA/TcMar");
        assert!(hit.2 > 0.99 && hit.3 > 0.99 && hit.4 - hit.5 > 0.10);
    }

    #[test]
    fn sample_manifest_has_unambiguous_three_and_four_column_forms() {
        let path = std::env::temp_dir().join(format!("main_repeat_samples_{}", std::process::id()));
        fs::write(
            &path,
            "# comment before the header
taxon_id sample_id read1 read2
t1 s1 a.fq b.fq
t2 s2 c.fq
",
        )
        .unwrap();
        let parsed = samples(&path).unwrap();
        fs::remove_file(path).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].id, "s1");
        assert_eq!(parsed[0].r2.as_deref(), Some(Path::new("b.fq")));
        assert_eq!(parsed[1].id, "s2");
        assert!(parsed[1].r2.is_none());
    }
}
