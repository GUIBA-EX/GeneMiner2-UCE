use crate::fasta::read_fasta;
use flate2::read::MultiGzDecoder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

type Row = BTreeMap<String, String>;

#[derive(Clone, Debug)]
pub struct Sample {
    pub name: String,
    pub reads: Vec<PathBuf>,
}

fn number(row: &Row, field: &str) -> f64 {
    row.get(field)
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0)
}

fn integer(row: &Row, field: &str) -> u64 {
    number(row, field) as u64
}

fn truth(value: Option<&String>) -> bool {
    matches!(value.map(String::as_str), Some("1" | "true" | "yes"))
}

pub fn accepted(row: &Row) -> bool {
    if row.is_empty() {
        return false;
    }
    if row
        .get("accepted")
        .is_some_and(|value| !value.trim().is_empty())
    {
        return truth(row.get("accepted"));
    }
    row.get("status").is_some_and(|value| value == "success") && !truth(row.get("low_quality"))
}

fn read_csv(path: &Path) -> io::Result<Vec<Row>> {
    if !path.is_file() {
        return Ok(Vec::new());
    }
    let mut reader = csv::Reader::from_path(path)?;
    let headers = reader.headers()?.clone();
    reader
        .records()
        .map(|record| {
            let record = record?;
            Ok(headers
                .iter()
                .zip(record.iter())
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect())
        })
        .collect()
}

fn fmt(value: f64) -> String {
    let mut text = format!("{value:.3}");
    while text.contains('.') && text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    text
}

fn mean(values: &[f64]) -> String {
    if values.is_empty() {
        String::new()
    } else {
        fmt(values.iter().sum::<f64>() / values.len() as f64)
    }
}

fn median(values: &[f64]) -> String {
    if values.is_empty() {
        return String::new();
    }
    let mut values = values.to_vec();
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        fmt(values[middle])
    } else {
        fmt((values[middle - 1] + values[middle]) / 2.0)
    }
}

fn reference_lengths(reference: &Path) -> io::Result<BTreeMap<String, u64>> {
    let mut result = BTreeMap::new();
    for entry in fs::read_dir(reference)? {
        let path = entry?.path();
        let extension = path
            .extension()
            .and_then(|value| value.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();
        if !matches!(extension.as_str(), "fa" | "fas" | "fasta") {
            continue;
        }
        let records = read_fasta(BufReader::new(File::open(&path)?))?;
        let lengths: Vec<u64> = records
            .iter()
            .map(|record| {
                record
                    .sequence
                    .iter()
                    .filter(|&&base| !matches!(base, b'-' | b'N' | b'n'))
                    .count() as u64
            })
            .collect();
        let average = if lengths.is_empty() {
            0
        } else {
            (lengths.iter().sum::<u64>() as f64 / lengths.len() as f64).round() as u64
        };
        result.insert(
            path.file_stem()
                .and_then(|value| value.to_str())
                .unwrap_or_default()
                .to_string(),
            average,
        );
    }
    Ok(result)
}

fn filtered_counts(path: &Path) -> io::Result<HashMap<String, u64>> {
    let mut result = HashMap::new();
    if !path.is_file() {
        return Ok(result);
    }
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut fields = line.split(',');
        if let (Some(name), Some(count)) = (fields.next(), fields.next()) {
            if !name.is_empty() {
                result.insert(name.to_string(), count.parse::<f64>().unwrap_or(0.0) as u64);
            }
        }
    }
    Ok(result)
}

fn fastq_reads(path: &Path) -> io::Result<u64> {
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|value| value == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::new(reader).lines().count() as u64 / 4)
}

fn write_table(path: &Path, headers: &[&str], rows: &[Row]) -> io::Result<()> {
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(headers)?;
    for row in rows {
        writer.write_record(
            headers
                .iter()
                .map(|field| row.get(*field).map(String::as_str).unwrap_or("")),
        )?;
    }
    writer.flush()
}

fn write_matrix(
    path: &Path,
    loci: &[String],
    samples: &[Sample],
    matrix: &[Vec<u64>],
    means: Option<&BTreeMap<String, u64>>,
) -> io::Result<()> {
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(std::iter::once("Species").chain(loci.iter().map(String::as_str)))?;
    if let Some(means) = means {
        writer.write_record(
            std::iter::once("MeanLength".to_string()).chain(
                loci.iter()
                    .map(|locus| means.get(locus).copied().unwrap_or(0).to_string()),
            ),
        )?;
    }
    for (sample, values) in samples.iter().zip(matrix) {
        writer.write_record(
            std::iter::once(sample.name.clone()).chain(values.iter().map(ToString::to_string)),
        )?;
    }
    writer.flush()
}

fn write_heatmap(path: &Path, matrix: &[Vec<u64>], denominators: Option<&[u64]>) -> io::Result<()> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return Ok(());
    }
    let rows = matrix.len();
    let columns = matrix[0].len();
    let scale = (1600 / columns.max(rows)).clamp(1, 8);
    let width = columns * scale;
    let height = rows * scale;
    let maximum = matrix.iter().flatten().copied().max().unwrap_or(1).max(1) as f64;
    let mut pixels = vec![255_u8; width * height * 3];
    for (row, values) in matrix.iter().enumerate() {
        for (column, &value) in values.iter().enumerate() {
            let ratio = denominators
                .map_or(value as f64 / maximum, |items| {
                    if items[column] == 0 {
                        0.0
                    } else {
                        value as f64 / items[column] as f64
                    }
                })
                .clamp(0.0, 1.0);
            let shade = (255.0 * (1.0 - ratio)) as u8;
            for y in row * scale..(row + 1) * scale {
                for x in column * scale..(column + 1) * scale {
                    let offset = (y * width + x) * 3;
                    pixels[offset..offset + 3].copy_from_slice(&[shade, shade, shade]);
                }
            }
        }
    }
    let file = File::create(path)?;
    let mut encoder = png::Encoder::new(file, width as u32, height as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    encoder
        .write_header()?
        .write_image_data(&pixels)
        .map_err(io::Error::other)
}

pub fn run(
    output: &Path,
    reference: &Path,
    samples: &[Sample],
    count_input: bool,
    heatmaps: bool,
) -> io::Result<()> {
    let refs = reference_lengths(reference)?;
    let loci: Vec<String> = refs.keys().cloned().collect();
    let global_rescue = output.join("uce_rescue_summary.csv");
    let rescue_rows = if global_rescue.is_file() {
        read_csv(&global_rescue)?
    } else {
        let mut rows = Vec::new();
        for sample in samples {
            rows.extend(read_csv(
                &output.join(&sample.name).join("uce_rescue_summary.csv"),
            )?);
        }
        rows
    };
    let mut rescue_counts: HashMap<(String, String), u64> = HashMap::new();
    for row in &rescue_rows {
        if let (Some(sample), Some(status)) = (row.get("sample"), row.get("rescue_status")) {
            *rescue_counts
                .entry((sample.clone(), status.clone()))
                .or_default() += 1;
        }
    }
    let mut lengths_matrix = Vec::new();
    let mut reads_matrix = Vec::new();
    let mut filtered_matrix = Vec::new();
    let mut sample_rows = Vec::new();
    let mut locus_values: Vec<Vec<Row>> = vec![Vec::new(); loci.len()];

    for sample in samples {
        let assembly_path = output.join(&sample.name).join("uce_assembly_summary.csv");
        let assembly: HashMap<String, Row> = read_csv(&assembly_path)?
            .into_iter()
            .filter_map(|row| row.get("locus").cloned().map(|locus| (locus, row)))
            .collect();
        let filtered =
            filtered_counts(&output.join(&sample.name).join("ref_reads_count_dict.txt"))?;
        let mut lengths = Vec::new();
        let mut read_counts = Vec::new();
        let mut statuses: HashMap<String, u64> = HashMap::new();
        let mut spans = Vec::new();
        let mut densities = Vec::new();
        for (index, locus) in loci.iter().enumerate() {
            let row = assembly.get(locus).cloned().unwrap_or_default();
            *statuses
                .entry(
                    row.get("status")
                        .cloned()
                        .unwrap_or_else(|| "missing".to_string()),
                )
                .or_default() += 1;
            let ok = accepted(&row);
            let length = if ok {
                integer(&row, "selected_contig_length")
            } else {
                0
            };
            let reads = if ok { integer(&row, "read_count") } else { 0 };
            if ok && integer(&row, "read_supported_span") > 0 {
                spans.push(integer(&row, "read_supported_span") as f64);
            }
            if length > 0 {
                densities.push(reads as f64 / length as f64);
            }
            lengths.push(length);
            read_counts.push(reads);
            locus_values[index].push(row);
        }
        let filtered_values: Vec<u64> = loci
            .iter()
            .map(|locus| filtered.get(locus).copied().unwrap_or(0))
            .collect();
        let filtered_total: u64 = filtered.values().sum();
        let input_reads = if count_input {
            sample
                .reads
                .iter()
                .filter(|path| path.is_file())
                .map(|path| fastq_reads(path))
                .collect::<io::Result<Vec<_>>>()?
                .iter()
                .sum()
        } else {
            0
        };
        let recovered: Vec<f64> = lengths
            .iter()
            .filter(|&&value| value > 0)
            .map(|&value| value as f64)
            .collect();
        let mut stats = Row::new();
        stats.insert("Name".into(), sample.name.clone());
        stats.insert(
            "InputReads".into(),
            if count_input {
                input_reads.to_string()
            } else {
                String::new()
            },
        );
        stats.insert("ReadsFiltered".into(), filtered_total.to_string());
        stats.insert(
            "PctFiltered".into(),
            if count_input && input_reads > 0 {
                fmt(filtered_total as f64 / input_reads as f64 * 100.0)
            } else {
                String::new()
            },
        );
        stats.insert(
            "LociWithFilteredReads".into(),
            filtered
                .values()
                .filter(|&&value| value > 0)
                .count()
                .to_string(),
        );
        stats.insert("LociWithContigs".into(), recovered.len().to_string());
        for (field, status) in [
            ("LociSuccess", "success"),
            ("LociLowQuality", "low quality"),
            ("LociMissing", "missing"),
            ("LociNoFilteredFile", "no filtered file"),
            ("LociNoSeed", "no seed"),
            ("LociNoContigs", "no contigs"),
            ("LociInsufficientGenomicKmers", "insufficient genomic kmers"),
        ] {
            stats.insert(
                field.into(),
                statuses.get(status).copied().unwrap_or(0).to_string(),
            );
        }
        for (field, threshold) in [
            ("LociAt25pct", 0.25),
            ("LociAt50pct", 0.5),
            ("LociAt75pct", 0.75),
            ("LociAt150pct", 1.5),
        ] {
            stats.insert(
                field.into(),
                loci.iter()
                    .enumerate()
                    .filter(|(i, locus)| {
                        refs[*locus] > 0 && lengths[*i] as f64 >= refs[*locus] as f64 * threshold
                    })
                    .count()
                    .to_string(),
            );
        }
        for (field, status) in [
            ("RescueSuccess", "success"),
            ("RescueFailedRolledBack", "failed_rolled_back"),
            ("RescueRevertedDensityDrop", "reverted_density_drop"),
            ("RescueRevertedFailed", "reverted_failed_rescue"),
        ] {
            stats.insert(
                field.into(),
                rescue_counts
                    .get(&(sample.name.clone(), status.to_string()))
                    .copied()
                    .unwrap_or(0)
                    .to_string(),
            );
        }
        stats.insert(
            "TotalBasesRecovered".into(),
            lengths.iter().sum::<u64>().to_string(),
        );
        stats.insert("MeanContigLength".into(), mean(&recovered));
        stats.insert("MedianContigLength".into(), median(&recovered));
        stats.insert("MeanReadSupportedSpan".into(), mean(&spans));
        stats.insert("MeanReadDensity".into(), mean(&densities));
        sample_rows.push(stats);
        lengths_matrix.push(lengths);
        reads_matrix.push(read_counts);
        filtered_matrix.push(filtered_values);
    }

    let sample_headers = [
        "Name",
        "InputReads",
        "ReadsFiltered",
        "PctFiltered",
        "LociWithFilteredReads",
        "LociWithContigs",
        "LociSuccess",
        "LociLowQuality",
        "LociMissing",
        "LociNoFilteredFile",
        "LociNoSeed",
        "LociNoContigs",
        "LociInsufficientGenomicKmers",
        "LociAt25pct",
        "LociAt50pct",
        "LociAt75pct",
        "LociAt150pct",
        "RescueSuccess",
        "RescueFailedRolledBack",
        "RescueRevertedDensityDrop",
        "RescueRevertedFailed",
        "TotalBasesRecovered",
        "MeanContigLength",
        "MedianContigLength",
        "MeanReadSupportedSpan",
        "MeanReadDensity",
    ];
    write_matrix(
        &output.join("uce_seq_lengths.tsv"),
        &loci,
        samples,
        &lengths_matrix,
        Some(&refs),
    )?;
    write_matrix(
        &output.join("uce_read_counts.tsv"),
        &loci,
        samples,
        &reads_matrix,
        None,
    )?;
    write_matrix(
        &output.join("uce_filtered_read_counts.tsv"),
        &loci,
        samples,
        &filtered_matrix,
        None,
    )?;
    write_table(&output.join("uce_stats.tsv"), &sample_headers, &sample_rows)?;

    let mut locus_rows = Vec::new();
    for (index, locus) in loci.iter().enumerate() {
        let rows = &locus_values[index];
        let accepted_rows: Vec<&Row> = rows.iter().filter(|row| accepted(row)).collect();
        let values = |field: &str| {
            accepted_rows
                .iter()
                .map(|row| number(row, field))
                .filter(|&value| value > 0.0)
                .collect::<Vec<_>>()
        };
        let lengths = values("selected_contig_length");
        let mut row = Row::new();
        row.insert("Locus".into(), locus.clone());
        row.insert("MeanReferenceLength".into(), refs[locus].to_string());
        row.insert("Samples".into(), samples.len().to_string());
        for (field, status) in [
            ("SuccessSamples", "success"),
            ("LowQualitySamples", "low quality"),
            ("NoFilteredFileSamples", "no filtered file"),
        ] {
            row.insert(
                field.into(),
                rows.iter()
                    .filter(|item| item.get("status").is_some_and(|value| value == status))
                    .count()
                    .to_string(),
            );
        }
        row.insert(
            "Occupancy".into(),
            if samples.is_empty() {
                String::new()
            } else {
                fmt(rows
                    .iter()
                    .filter(|item| item.get("status").is_some_and(|value| value == "success"))
                    .count() as f64
                    / samples.len() as f64)
            },
        );
        row.insert("MeanLength".into(), mean(&lengths));
        row.insert("MedianLength".into(), median(&lengths));
        row.insert(
            "MaxLength".into(),
            lengths
                .iter()
                .copied()
                .max_by(f64::total_cmp)
                .map(fmt)
                .unwrap_or_default(),
        );
        for (out, field) in [
            ("MeanReadCount", "read_count"),
            ("MeanReadSupportedSpan", "read_supported_span"),
            ("MeanFlankBalance", "flank_balance"),
            ("MeanCandidateCount", "candidate_count"),
        ] {
            row.insert(out.into(), mean(&values(field)));
        }
        locus_rows.push(row);
    }
    let locus_headers = [
        "Locus",
        "MeanReferenceLength",
        "Samples",
        "SuccessSamples",
        "LowQualitySamples",
        "NoFilteredFileSamples",
        "Occupancy",
        "MeanLength",
        "MedianLength",
        "MaxLength",
        "MeanReadCount",
        "MeanReadSupportedSpan",
        "MeanFlankBalance",
        "MeanCandidateCount",
    ];
    write_table(
        &output.join("uce_locus_stats.tsv"),
        &locus_headers,
        &locus_rows,
    )?;
    if !rescue_rows.is_empty() {
        let headers: BTreeSet<String> = rescue_rows
            .iter()
            .flat_map(|row| row.keys().cloned())
            .collect();
        let headers: Vec<String> = headers.into_iter().collect();
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(output.join("uce_rescue_stats.tsv"))?;
        writer.write_record(&headers)?;
        for row in &rescue_rows {
            writer.write_record(
                headers
                    .iter()
                    .map(|field| row.get(field).map(String::as_str).unwrap_or("")),
            )?;
        }
        writer.flush()?;
    }
    if heatmaps {
        let denominators: Vec<u64> = loci.iter().map(|locus| refs[locus]).collect();
        write_heatmap(
            &output.join("uce_recovery_heatmap.png"),
            &lengths_matrix,
            Some(&denominators),
        )?;
        write_heatmap(
            &output.join("uce_read_counts_heatmap.png"),
            &reads_matrix,
            None,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn explicit_acceptance_wins() {
        let mut row = Row::new();
        row.insert("accepted".into(), "false".into());
        row.insert("status".into(), "success".into());
        assert!(!accepted(&row));
    }
    #[test]
    fn formatting_matches_python_tables() {
        assert_eq!(fmt(1.2304), "1.23");
    }
}
