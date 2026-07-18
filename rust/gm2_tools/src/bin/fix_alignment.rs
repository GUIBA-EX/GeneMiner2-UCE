use gm2_tools::alignment::clean_file;
use std::env;
use std::path::PathBuf;
use std::process;

struct Args {
    file: PathBuf,
    minimum_sequences: usize,
    maximum_difference: f64,
}

fn parse_args() -> Result<Args, String> {
    let mut file = None;
    let mut minimum_sequences = 1;
    let mut maximum_difference = 1.0;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-f" | "--file" => file = Some(PathBuf::from(args.next().ok_or("-f requires a path")?)),
            "-n" => {
                minimum_sequences = args
                    .next()
                    .ok_or("-n requires an integer")?
                    .parse()
                    .map_err(|_| "invalid -n value")?;
            }
            "-p" => {
                maximum_difference = args
                    .next()
                    .ok_or("-p requires a number")?
                    .parse()
                    .map_err(|_| "invalid -p value")?;
            }
            "-h" | "--help" => {
                println!("Usage: fix_alignment -f FILE [-n MIN_SEQUENCES] [-p MAX_DIFFERENCE]");
                process::exit(0);
            }
            _ => return Err(format!("unknown argument: {arg}")),
        }
    }
    if !(0.0..=1.0).contains(&maximum_difference) {
        return Err("-p must be between 0 and 1".to_string());
    }
    Ok(Args {
        file: file.ok_or("-f is required")?,
        minimum_sequences,
        maximum_difference,
    })
}

fn main() {
    let args = parse_args().unwrap_or_else(|error| {
        eprintln!("Invalid argument: {error}");
        process::exit(2);
    });
    if let Err(error) = clean_file(&args.file, args.minimum_sequences, args.maximum_difference) {
        eprintln!("Unable to clean {}: {error}", args.file.display());
        process::exit(1);
    }
}
