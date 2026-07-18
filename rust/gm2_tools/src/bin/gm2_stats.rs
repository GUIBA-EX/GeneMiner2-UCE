use gm2_tools::stats::{run, Sample};
use std::env;
use std::path::PathBuf;
use std::process;

fn main() {
    let mut output = None;
    let mut reference = None;
    let mut samples = Vec::new();
    let mut count_input = false;
    let mut heatmaps = true;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" => output = args.next().map(PathBuf::from),
            "--reference" => reference = args.next().map(PathBuf::from),
            "--sample" => {
                let name = args.next();
                let first = args.next();
                let second = args.next();
                match (name, first, second) {
                    (Some(name), Some(first), Some(second)) => samples.push(Sample {
                        name,
                        reads: vec![PathBuf::from(first), PathBuf::from(second)],
                    }),
                    _ => {
                        eprintln!("--sample requires NAME READ1 READ2");
                        process::exit(2);
                    }
                }
            }
            "--count-input-reads" => count_input = true,
            "--no-heatmap" => heatmaps = false,
            "-h" | "--help" => {
                println!("Usage: gm2_stats --output DIR --reference DIR [--sample NAME READ1 READ2] [--count-input-reads] [--no-heatmap]");
                return;
            }
            _ => {
                eprintln!("Unknown argument: {arg}");
                process::exit(2);
            }
        }
    }
    let required = |value: Option<PathBuf>, name: &str| {
        value.unwrap_or_else(|| {
            eprintln!("{name} is required");
            process::exit(2);
        })
    };
    let output = required(output, "--output");
    let reference = required(reference, "--reference");
    if let Err(error) = run(&output, &reference, &samples, count_input, heatmaps) {
        eprintln!("Unable to write UCE statistics: {error}");
        process::exit(1);
    }
    println!("Wrote UCE statistics to {}", output.display());
}
