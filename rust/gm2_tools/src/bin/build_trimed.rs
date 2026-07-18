use gm2_tools::trim::{run_trim, TrimMode};
use std::env;
use std::path::PathBuf;
use std::process;

fn main() {
    let mut input = None;
    let mut reference = None;
    let mut output = None;
    let mut database = None;
    let mut executable = None;
    let mut mode = 2_u8;
    let mut percentage = 50.0_f64;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if matches!(arg.as_str(), "-h" | "--help") {
            println!("Usage: build_trimed -i QUERY -r REF -o OUTPUT -b DB --executable PATH [-m 0..3] [-p PERCENT]");
            return;
        }
        let value = args.next().unwrap_or_else(|| {
            eprintln!("{arg} requires a value");
            process::exit(2);
        });
        match arg.as_str() {
            "-i" | "--input" => input = Some(PathBuf::from(value)),
            "-r" | "--ref" => reference = Some(PathBuf::from(value)),
            "-o" | "--output" => output = Some(PathBuf::from(value)),
            "-b" | "--blast-db" => database = Some(PathBuf::from(value)),
            "--executable" => executable = Some(PathBuf::from(value)),
            "-m" | "--mode" => mode = value.parse().unwrap_or(255),
            "-p" | "--pec" => percentage = value.parse().unwrap_or(f64::NAN),
            _ => {
                eprintln!("Unknown argument: {arg}");
                process::exit(2);
            }
        }
    }
    let mode = match mode {
        0 => TrimMode::All,
        1 => TrimMode::Longest,
        2 => TrimMode::Terminal,
        3 => TrimMode::Isoform,
        _ => {
            eprintln!("-m must be 0, 1, 2, or 3");
            process::exit(2);
        }
    };
    if !percentage.is_finite() {
        eprintln!("-p must be a number");
        process::exit(2);
    }
    let required = |value: Option<PathBuf>, name: &str| {
        value.unwrap_or_else(|| {
            eprintln!("{name} is required");
            process::exit(2);
        })
    };
    let input = required(input, "-i");
    let reference = required(reference, "-r");
    let output = required(output, "-o");
    let database = required(database, "-b");
    let executable = required(executable, "--executable");
    if let Err(error) = run_trim(
        &input,
        &reference,
        &output,
        &database,
        &executable,
        percentage,
        mode,
    ) {
        eprintln!("Unable to trim {}: {error}", input.display());
        process::exit(1);
    }
}
