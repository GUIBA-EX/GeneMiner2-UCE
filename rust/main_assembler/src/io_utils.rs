use crate::model::{LocusTask, RefKmer};
use crate::seq::{encode_kmer, reverse_complement_kmer, valid_runs};
use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static CACHE_TEMP_ID: AtomicU64 = AtomicU64::new(0);
const CACHE_MAGIC: &[u8; 8] = b"GM2RK001";

fn clean_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .filter(|base| base.is_ascii_alphabetic())
        .map(|base| base.to_ascii_uppercase())
        .collect()
}

pub fn read_fasta(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let text = fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut title: Option<String> = None;
    let mut sequence = Vec::new();

    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = title.take() {
                records.push((previous, clean_sequence(&sequence)));
                sequence.clear();
            }
            title = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else if title.is_some() {
            sequence.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(previous) = title {
        records.push((previous, clean_sequence(&sequence)));
    }
    Ok(records)
}

pub fn read_sequences(path: &Path, fasta: bool) -> io::Result<Vec<Vec<u8>>> {
    if fasta {
        return Ok(read_fasta(path)?
            .into_iter()
            .map(|(_, sequence)| sequence)
            .collect());
    }

    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();
    if !lines.len().is_multiple_of(4) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("truncated FASTQ record in {}", path.display()),
        ));
    }
    let mut sequences = Vec::with_capacity(lines.len() / 4);
    for record in lines.chunks_exact(4) {
        if !record[0].starts_with('@') || !record[2].starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FASTQ record in {}", path.display()),
            ));
        }
        sequences.push(clean_sequence(record[1].as_bytes()));
    }
    Ok(sequences)
}

pub fn discover_references(reference: &Path) -> io::Result<Vec<LocusTask>> {
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference)? {
            let path = entry?.path();
            if path.is_file() && is_fasta(&path) {
                paths.push(path);
            }
        }
    } else if reference.is_file() && is_fasta(reference) {
        paths.push(reference.to_path_buf());
    }
    paths.sort();

    let total = paths.len();
    let mut tasks = Vec::with_capacity(total);
    for (index, path) in paths.into_iter().enumerate() {
        let key = path
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("")
            .to_string();
        let reference_count = read_fasta(&path)?.len();
        tasks.push(LocusTask {
            key,
            reference_path: path,
            reference_count,
            ordinal: index + 1,
            total,
        });
    }
    Ok(tasks)
}

pub fn is_fasta(path: &Path) -> bool {
    path.extension()
        .and_then(|value| value.to_str())
        .is_some_and(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "fa" | "fas" | "fasta"
            )
        })
}

pub fn find_filtered(output: &Path, key: &str) -> Option<(PathBuf, bool)> {
    let fasta = output.join("filtered").join(format!("{key}.fasta"));
    if fasta.is_file() {
        return Some((fasta, true));
    }
    let fastq = output.join("filtered").join(format!("{key}.fq"));
    if fastq.is_file() {
        return Some((fastq, false));
    }
    None
}

pub fn build_reference_kmers(records: &[(String, Vec<u8>)], k: usize) -> HashMap<u128, RefKmer> {
    let mut kmers = HashMap::new();
    for (_, sequence) in records {
        let total_kmers = sequence.len().saturating_sub(k) + 1;
        if total_kmers == 0 {
            continue;
        }
        for (run_start, run) in valid_runs(sequence) {
            if run.len() < k {
                continue;
            }

            // Python visits forward suffixes first.
            for start in (0..=run.len() - k).rev() {
                let encoded = encode_kmer(&run[start..start + k]).expect("valid run");
                let j = run.len() - k - start;
                let global_j = sequence.len() - run_start - run.len() + j;
                let position = (((global_j + 1) as f64 / total_kmers as f64) * 1000.0) as i32;
                insert_reference(&mut kmers, encoded, position, false);
            }

            // Then it visits reverse-complement suffixes from the original prefix.
            for start in 0..=run.len() - k {
                let encoded = encode_kmer(&run[start..start + k]).expect("valid run");
                let reverse = reverse_complement_kmer(encoded, k);
                let global_j = run_start + start;
                let position = (((global_j + 1) as f64 / total_kmers as f64) * 1000.0) as i32;
                insert_reference(&mut kmers, reverse, position, true);
            }
        }
    }
    kmers
}

fn insert_reference(
    kmers: &mut HashMap<u128, RefKmer>,
    kmer: u128,
    position: i32,
    is_reverse: bool,
) {
    if let Some(value) = kmers.get_mut(&kmer) {
        value.depth = value.depth.saturating_add(1);
    } else {
        kmers.insert(
            kmer,
            RefKmer {
                depth: 1,
                position,
                is_reverse,
            },
        );
    }
}

fn cache_path(cache_dir: &Path, reference: &Path, k: usize) -> io::Result<PathBuf> {
    let metadata = fs::metadata(reference)?;
    let modified = metadata
        .modified()
        .ok()
        .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    let mut hasher = DefaultHasher::new();
    reference
        .canonicalize()
        .unwrap_or_else(|_| reference.to_path_buf())
        .hash(&mut hasher);
    metadata.len().hash(&mut hasher);
    modified.hash(&mut hasher);
    k.hash(&mut hasher);
    let name = reference
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("reference");
    Ok(cache_dir.join(format!("{name}.k{k}.{:016x}.gm2rk", hasher.finish())))
}

pub fn load_or_build_reference_kmers(
    reference: &Path,
    records: &[(String, Vec<u8>)],
    k: usize,
    cache_dir: Option<&Path>,
) -> io::Result<HashMap<u128, RefKmer>> {
    let Some(cache_dir) = cache_dir else {
        return Ok(build_reference_kmers(records, k));
    };
    fs::create_dir_all(cache_dir)?;
    let cache = cache_path(cache_dir, reference, k)?;
    if let Ok(kmers) = read_reference_cache(&cache, k) {
        return Ok(kmers);
    }

    let kmers = build_reference_kmers(records, k);
    let temp_id = CACHE_TEMP_ID.fetch_add(1, Ordering::Relaxed);
    let temporary = cache.with_extension(format!("{}.{}.tmp", std::process::id(), temp_id));
    write_reference_cache(&temporary, k, &kmers)?;
    fs::rename(&temporary, &cache)?;
    Ok(kmers)
}

fn read_exact_array<const N: usize>(reader: &mut impl Read) -> io::Result<[u8; N]> {
    let mut bytes = [0_u8; N];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_reference_cache(path: &Path, expected_k: usize) -> io::Result<HashMap<u128, RefKmer>> {
    let mut reader = File::open(path)?;
    if &read_exact_array::<8>(&mut reader)? != CACHE_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cache magic mismatch",
        ));
    }
    let k = u32::from_le_bytes(read_exact_array::<4>(&mut reader)?) as usize;
    if k != expected_k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cache k mismatch",
        ));
    }
    let count = u64::from_le_bytes(read_exact_array::<8>(&mut reader)?) as usize;
    let mut kmers = HashMap::with_capacity(count);
    for _ in 0..count {
        let key = u128::from_le_bytes(read_exact_array::<16>(&mut reader)?);
        let depth = u32::from_le_bytes(read_exact_array::<4>(&mut reader)?);
        let position = i32::from_le_bytes(read_exact_array::<4>(&mut reader)?);
        let reverse = read_exact_array::<1>(&mut reader)?[0] != 0;
        kmers.insert(
            key,
            RefKmer {
                depth,
                position,
                is_reverse: reverse,
            },
        );
    }
    Ok(kmers)
}

fn write_reference_cache(path: &Path, k: usize, kmers: &HashMap<u128, RefKmer>) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(CACHE_MAGIC)?;
    writer.write_all(&(k as u32).to_le_bytes())?;
    writer.write_all(&(kmers.len() as u64).to_le_bytes())?;
    for (key, value) in kmers {
        writer.write_all(&key.to_le_bytes())?;
        writer.write_all(&value.depth.to_le_bytes())?;
        writer.write_all(&value.position.to_le_bytes())?;
        writer.write_all(&[u8::from(value.is_reverse)])?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_contains_both_orientations() {
        let records = vec![("x".to_string(), b"AACCGGTT".to_vec())];
        let kmers = build_reference_kmers(&records, 4);
        let forward = encode_kmer(b"AACC").unwrap();
        let reverse = reverse_complement_kmer(forward, 4);
        assert!(kmers.contains_key(&forward));
        assert!(kmers.contains_key(&reverse));
    }
}
