use std::process::Command;

fn main() {
    println!("cargo:rustc-check-cfg=cfg(system_zlib_ng)");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_LIBDIR");

    let Ok(output) = Command::new("pkg-config")
        .args(["--libs", "zlib-ng"])
        .output()
    else {
        return;
    };
    if !output.status.success() {
        return;
    }

    for flag in String::from_utf8_lossy(&output.stdout).split_whitespace() {
        if let Some(path) = flag.strip_prefix("-L") {
            println!("cargo:rustc-link-search=native={path}");
        } else if let Some(library) = flag.strip_prefix("-l") {
            println!("cargo:rustc-link-lib={library}");
        }
    }
    println!("cargo:rustc-cfg=system_zlib_ng");
    println!("cargo:warning=MainFilter: native zlib-ng detected through pkg-config");
}
