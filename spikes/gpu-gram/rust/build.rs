fn main() {
    println!("cargo:rerun-if-env-changed=OPENBLAS_LIB_DIR");
    if let Ok(dir) = std::env::var("OPENBLAS_LIB_DIR") {
        println!("cargo:rustc-link-search=native={dir}");
    }
}
