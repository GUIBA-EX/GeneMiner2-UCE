use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug)]
pub(crate) struct BaitLocus {
    pub(crate) name: String,
    pub(crate) records: Vec<Vec<u8>>,
}

#[derive(Clone, Debug)]
pub(crate) struct BaitCatalog {
    pub(crate) loci: Vec<BaitLocus>,
}

impl BaitCatalog {
    pub(crate) fn read(directory: &Path) -> Result<Self, String> {
        let mut paths = fs::read_dir(directory)
            .map_err(|e| format!("cannot read bait directory {}: {e}", directory.display()))?
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| path.is_file())
            .filter(|path| {
                matches!(
                    path.extension().and_then(|x| x.to_str()),
                    Some("fa" | "fas" | "fasta")
                )
            })
            .collect::<Vec<_>>();
        paths.sort();
        if paths.is_empty() {
            return Err("--panref-baits contains no FASTA files".into());
        }
        let loci = paths
            .into_iter()
            .map(|path| {
                let name = path
                    .file_stem()
                    .and_then(|value| value.to_str())
                    .ok_or_else(|| format!("invalid bait name: {}", path.display()))?
                    .to_string();
                let records = read_records(&path)?;
                if records.is_empty() {
                    return Err(format!("bait locus {name} contains no sequence"));
                }
                Ok(BaitLocus { name, records })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Self { loci })
    }
}

fn read_records(path: &Path) -> Result<Vec<Vec<u8>>, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut records = Vec::new();
    let mut sequence = Vec::new();
    let mut seen_header = false;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('>') {
            if seen_header && !sequence.is_empty() {
                records.push(std::mem::take(&mut sequence));
            }
            seen_header = true;
            continue;
        }
        if !seen_header {
            return Err(format!("FASTA sequence before header: {}", path.display()));
        }
        for base in line.trim().bytes() {
            match base.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => sequence.push(base.to_ascii_uppercase()),
                b'U' => sequence.push(b'T'),
                _ => return Err(format!("bait {} contains non-ACGT base", path.display())),
            }
        }
    }
    if !sequence.is_empty() {
        records.push(sequence);
    }
    Ok(records)
}
