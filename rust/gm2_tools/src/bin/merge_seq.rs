use gm2_tools::merge::merge_sequences;
use std::env;
use std::path::PathBuf;
use std::process;

fn main() {
    let mut input = None;
    let mut output = PathBuf::from("merge.fasta");
    let mut extensions = ".fasta,.fas,.fa".to_string();
    let mut missing = b'N';
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        let value = match arg.as_str() {
            "-input" | "-output" | "-exts" | "-missing" => args.next().unwrap_or_else(|| {
                eprintln!("{arg} requires a value");
                process::exit(2);
            }),
            "-h" | "--help" => {
                println!("Usage: merge_seq -input DIR [-output FILE] [-exts LIST] [-missing CHAR]");
                return;
            }
            _ => {
                eprintln!("Unknown argument: {arg}");
                process::exit(2);
            }
        };
        match arg.as_str() {
            "-input" => input = Some(PathBuf::from(value)),
            "-output" => output = PathBuf::from(value),
            "-exts" => extensions = value,
            "-missing" => {
                let bytes = value.as_bytes();
                if bytes.len() != 1 {
                    eprintln!("-missing requires one ASCII character");
                    process::exit(2);
                }
                missing = bytes[0];
            }
            _ => unreachable!(),
        }
    }
    let input = input.unwrap_or_else(|| {
        eprintln!("-input is required");
        process::exit(2);
    });
    if let Err(error) = merge_sequences(&input, &output, &extensions, missing) {
        eprintln!("Unable to merge {}: {error}", input.display());
        process::exit(1);
    }
    println!("Merging completed.");
}
