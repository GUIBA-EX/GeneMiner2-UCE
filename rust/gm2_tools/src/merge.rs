use crate::fasta::{read_fasta, FastaRecord};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct Alignment {
    locus: String,
    length: usize,
    sequences: HashMap<String, Vec<u8>>,
}

fn species_sort_key(name: &str) -> (i64, &str) {
    match name.split_once('_') {
        Some((number, suffix)) => number.parse().map_or((0, name), |number| (number, suffix)),
        None => (0, name),
    }
}

pub fn partition_path(output: &Path) -> PathBuf {
    let mut path = output.to_path_buf();
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or_default();
    path.set_file_name(format!("{stem}_partition.txt"));
    path
}

pub fn merge_sequences(
    input_folder: &Path,
    output_file: &Path,
    extensions: &str,
    missing: u8,
) -> io::Result<()> {
    let allowed: HashSet<String> = extensions
        .split(',')
        .map(|value| value.to_ascii_lowercase())
        .collect();
    let mut paths = Vec::new();
    for entry in fs::read_dir(input_folder)? {
        let entry = entry?;
        if !entry.file_type()?.is_file() {
            continue;
        }
        let path = entry.path();
        let extension = path
            .extension()
            .and_then(|value| value.to_str())
            .map(|value| format!(".{}", value.to_ascii_lowercase()));
        if extension.is_some_and(|value| allowed.contains(&value)) {
            paths.push(path);
        }
    }
    paths.sort_by(|left, right| left.file_name().cmp(&right.file_name()));

    let mut alignments = Vec::new();
    let mut species = HashSet::new();
    for path in paths {
        let records = read_fasta(BufReader::new(File::open(&path)?))?;
        let length = records
            .first()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty alignment"))?
            .sequence
            .len();
        let mut sequences = HashMap::new();
        for FastaRecord { name, sequence } in records {
            species.insert(name.clone());
            sequences.insert(name, sequence);
        }
        alignments.push(Alignment {
            locus: path
                .file_stem()
                .and_then(|value| value.to_str())
                .unwrap_or_default()
                .to_string(),
            length,
            sequences,
        });
    }

    let mut species: Vec<String> = species.into_iter().collect();
    species.sort_by(|left, right| {
        species_sort_key(left)
            .cmp(&species_sort_key(right))
            .then_with(|| left.cmp(right))
    });

    let mut output = BufWriter::new(File::create(output_file)?);
    for name in species {
        writeln!(output, ">{name}")?;
        for alignment in &alignments {
            if let Some(sequence) = alignment.sequences.get(&name) {
                output.write_all(sequence)?;
            } else {
                output.write_all(&vec![missing; alignment.length])?;
            }
        }
        output.write_all(b"\n")?;
    }

    let mut partitions = BufWriter::new(File::create(partition_path(output_file))?);
    partitions.write_all(b"#nexus\nbegin sets;\n")?;
    let mut end = 0_usize;
    for (index, alignment) in alignments.iter().enumerate() {
        if alignment.length == 0 {
            continue;
        }
        let start = end + 1;
        end += alignment.length;
        writeln!(
            partitions,
            "charset part{}_{} = {}-{};",
            index + 1,
            alignment.locus,
            start,
            end
        )?;
    }
    partitions.write_all(b"end;\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_suffix_replaces_fasta_extension() {
        assert_eq!(
            partition_path(Path::new("result.fasta")),
            PathBuf::from("result_partition.txt")
        );
    }

    #[test]
    fn numeric_sample_prefix_controls_order() {
        let mut names = ["10_B", "2_A", "sample"];
        names.sort_by(|left, right| species_sort_key(left).cmp(&species_sort_key(right)));
        assert_eq!(names, ["sample", "2_A", "10_B"]);
    }
}
