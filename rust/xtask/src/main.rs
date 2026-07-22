use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const BINARIES: &[(&str, &str)] = &[
    ("build_consensus", "build_consensus"),
    ("MainFilterNew", "MainFilterNew"),
    ("main_refilter_new", "main_refilter_new"),
    ("uce_filter", "uce_filter"),
    ("main_assembler_original", "main_assembler-original-rust"),
    ("main_assembler", "main_assembler-rust"),
    ("main_population", "main_population"),
    ("fix_alignment", "fix_alignment"),
    ("merge_seq", "merge_seq"),
    ("build_trimed", "build_trimed"),
    ("gm2_stats", "gm2_stats"),
    ("mito_workflow", "mito_workflow"),
    ("gene_workflow", "gene_workflow"),
    ("rad_workflow", "rad_workflow"),
    ("marker_profile", "marker_profile"),
    ("main_repeat", "main_repeat"),
    ("geneminer2_cli", "geneminer2-rust"),
];

fn usage() {
    eprintln!("Usage: cargo run -p xtask -- <build|clean>");
}

fn workspace_root() -> io::Result<PathBuf> {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .map(Path::to_path_buf)
        .ok_or_else(|| io::Error::other("xtask is not inside a Cargo workspace"))
}

fn run_cargo(root: &Path, args: &[&str]) -> io::Result<()> {
    let status = Command::new("cargo")
        .args(args)
        .current_dir(root)
        .status()?;
    if status.success() {
        Ok(())
    } else {
        Err(io::Error::other("Cargo build failed"))
    }
}

#[cfg(unix)]
fn make_entrypoint(target: &Path, entrypoint: &Path) -> io::Result<()> {
    use std::os::unix::fs::symlink;
    symlink(target, entrypoint)
}

#[cfg(not(unix))]
fn make_entrypoint(target: &Path, entrypoint: &Path) -> io::Result<()> {
    std::os::windows::fs::symlink_file(target, entrypoint)
}

#[cfg(unix)]
fn mark_executable(path: &Path) -> io::Result<()> {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(path, fs::Permissions::from_mode(0o755))
}

#[cfg(not(unix))]
fn mark_executable(_: &Path) -> io::Result<()> {
    Ok(())
}

fn build(root: &Path) -> io::Result<()> {
    run_cargo(root, &["build", "--release", "--workspace", "--locked"])?;

    let cli = root.join("cli");
    let bin_dir = cli.join("bin");
    if bin_dir.exists() {
        fs::remove_dir_all(&bin_dir)?;
    }
    fs::create_dir_all(&bin_dir)?;

    let release_dir = root.join("target/release");
    for (source, destination) in BINARIES {
        let source = release_dir.join(source);
        let destination = bin_dir.join(destination);
        fs::copy(&source, &destination)?;
        mark_executable(&destination)?;
    }

    let entrypoint = cli.join("geneminer2");
    if entrypoint.exists() || entrypoint.is_symlink() {
        fs::remove_file(&entrypoint)?;
    }
    make_entrypoint(Path::new("bin/geneminer2-rust"), &entrypoint)
}

fn clean(root: &Path) -> io::Result<()> {
    run_cargo(root, &["clean"])?;
    let cli = root.join("cli");
    let bin_dir = cli.join("bin");
    if bin_dir.exists() {
        fs::remove_dir_all(bin_dir)?;
    }
    let entrypoint = cli.join("geneminer2");
    if entrypoint.exists() || entrypoint.is_symlink() {
        fs::remove_file(entrypoint)?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let root = match workspace_root() {
        Ok(root) => root,
        Err(error) => {
            eprintln!("Unable to locate workspace: {error}");
            return ExitCode::FAILURE;
        }
    };
    let result = match env::args().nth(1).as_deref() {
        Some("build") => build(&root),
        Some("clean") => clean(&root),
        _ => {
            usage();
            return ExitCode::FAILURE;
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("xtask failed: {error}");
            ExitCode::FAILURE
        }
    }
}
