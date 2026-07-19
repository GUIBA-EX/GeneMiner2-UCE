use super::bait::BaitCatalog;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::Path;

const FAST_SEED_LEN: usize = 23;
const MINIMIZER_WINDOW: usize = 12;
const AMBIGUOUS_LOCUS: u32 = u32::MAX;

#[derive(Clone, Debug)]
pub(crate) struct BaitRecord {
    pub(crate) locus: u32,
    pub(crate) sequence: Vec<u8>,
}

/// Sparse, unique bait minimizers for fast read-to-locus assignment.
///
/// This deliberately contains no FM-index, BWT, suffix array, or approximate
/// matching fallback: PanRefV2 uses a strict seed-and-core-graph policy.
#[derive(Debug)]
pub(crate) struct BaitIndex {
    fast_seeds: HashMap<u64, u32>,
    locus_count: usize,
}

impl BaitIndex {
    pub(crate) fn build_catalog(catalog: &BaitCatalog) -> Result<Self, String> {
        let records = catalog
            .loci
            .iter()
            .enumerate()
            .flat_map(|(locus, entry)| {
                entry
                    .records
                    .iter()
                    .cloned()
                    .map(move |sequence| BaitRecord {
                        locus: locus as u32,
                        sequence,
                    })
            })
            .collect::<Vec<_>>();
        Self::build(&records)
    }

    pub(crate) fn build(records: &[BaitRecord]) -> Result<Self, String> {
        if records.is_empty() {
            return Err("bait index requires at least one bait record".into());
        }
        let mut fast_seeds = HashMap::new();
        let mut locus_count = 0_usize;
        for record in records {
            locus_count = locus_count.max(record.locus as usize + 1);
            for seed in minimizer_codes(&record.sequence) {
                match fast_seeds.get_mut(&seed) {
                    None => {
                        fast_seeds.insert(seed, record.locus);
                    }
                    Some(existing) if *existing != record.locus => *existing = AMBIGUOUS_LOCUS,
                    Some(_) => {}
                }
            }
        }
        if fast_seeds.is_empty() {
            return Err(format!(
                "bait index needs at least one unambiguous {FAST_SEED_LEN}-mer"
            ));
        }
        Ok(Self {
            fast_seeds,
            locus_count,
        })
    }

    pub(crate) fn write_metadata(&self, path: &Path, loci: &[String]) -> Result<(), String> {
        let mut out = File::create(path).map_err(|e| e.to_string())?;
        writeln!(out, "locus_count\tunique_minimizers").map_err(|e| e.to_string())?;
        writeln!(out, "{}\t{}", loci.len(), self.fast_seeds.len()).map_err(|e| e.to_string())?;
        for (id, locus) in loci.iter().enumerate() {
            writeln!(out, "locus\t{}\t{}", id, locus).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub(crate) fn minimizer_loci(&self, sequence: &[u8]) -> Vec<u32> {
        minimizer_codes(sequence)
            .into_iter()
            .filter_map(|seed| self.fast_seeds.get(&seed).copied())
            .filter(|&locus| locus != AMBIGUOUS_LOCUS)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    pub(crate) fn locus_count(&self) -> usize {
        self.locus_count
    }
}

fn minimizer_codes(sequence: &[u8]) -> Vec<u64> {
    if sequence.len() < FAST_SEED_LEN {
        return Vec::new();
    }
    let mask = (1_u64 << (2 * FAST_SEED_LEN)) - 1;
    let mut codes = Vec::with_capacity(sequence.len() - FAST_SEED_LEN + 1);
    let mut code = 0_u64;
    for (position, &base) in sequence.iter().enumerate() {
        let bits = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => return Vec::new(),
        };
        code = ((code << 2) | bits) & mask;
        if position + 1 >= FAST_SEED_LEN {
            codes.push(code);
        }
    }
    let mut out = Vec::new();
    let mut previous = None;
    for window in codes.windows(MINIMIZER_WINDOW.min(codes.len())) {
        let seed = *window.iter().min().expect("nonempty minimizer window");
        if previous != Some(seed) {
            out.push(seed);
            previous = Some(seed);
        }
    }
    for seed in [
        codes[0],
        codes[codes.len() / 2],
        *codes.last().expect("nonempty codes"),
    ] {
        if !out.contains(&seed) {
            out.push(seed);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::{BaitIndex, BaitRecord};

    #[test]
    fn unique_minimizers_assign_a_single_locus() {
        let first = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT".to_vec();
        let second = b"TTGCAAGCTTAGCGATCGTACCTGAGTTCGATACCGTGAC".to_vec();
        let index = BaitIndex::build(&[
            BaitRecord {
                locus: 0,
                sequence: first.clone(),
            },
            BaitRecord {
                locus: 1,
                sequence: second.clone(),
            },
        ])
        .unwrap();
        assert_eq!(index.minimizer_loci(&first), vec![0]);
        assert_eq!(index.minimizer_loci(&second), vec![1]);
        assert_eq!(index.locus_count(), 2);
    }
}
