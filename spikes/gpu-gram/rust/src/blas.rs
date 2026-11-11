//! `cblas_dsyrk` and `cblas_ssyrk` on the dosage matrix.
//!
//! macOS links Accelerate. Windows links OpenBLAS. The matrix is column-major and site-major,
//! which is the same order as the generator: element `(taxon, site)` lives at `site * taxa + taxon`.
//! One full `syrk` per shape, not one call per panel.

use std::time::Instant;

use gpu_gram::{compare_triangles, pair_count, MISSING_CODE, RELATIVE_FLOOR};

use crate::report::{stats_of, Stats};

const RUNS: usize = 3;
#[cfg(any(target_os = "macos", target_os = "windows"))]
const CBLAS_COL_MAJOR: i32 = 102;
#[cfg(any(target_os = "macos", target_os = "windows"))]
const CBLAS_NO_TRANS: i32 = 111;
#[cfg(any(target_os = "macos", target_os = "windows"))]
const CBLAS_UPPER: i32 = 121;

pub struct BlasInput<'a> {
    pub taxa: usize,
    pub sites: usize,
    pub codes: &'a [u8],
    pub reference: &'a [f64],
}

pub struct BlasOutcome {
    pub kernel: String,
    pub threads: String,
    pub pack: Stats,
    pub syrk: Stats,
    pub max_abs: f64,
    pub max_rel: f64,
    pub bit_exact: bool,
}

#[cfg_attr(target_os = "macos", link(name = "Accelerate", kind = "framework"))]
#[cfg_attr(target_os = "windows", link(name = "openblas"))]
#[cfg(any(target_os = "macos", target_os = "windows"))]
unsafe extern "C" {
    fn cblas_dsyrk(
        order: i32,
        uplo: i32,
        trans: i32,
        n: i32,
        k: i32,
        alpha: f64,
        a: *const f64,
        lda: i32,
        beta: f64,
        c: *mut f64,
        ldc: i32,
    );
    fn cblas_ssyrk(
        order: i32,
        uplo: i32,
        trans: i32,
        n: i32,
        k: i32,
        alpha: f32,
        a: *const f32,
        lda: i32,
        beta: f32,
        c: *mut f32,
        ldc: i32,
    );
}

#[cfg(target_os = "windows")]
#[link(name = "openblas")]
unsafe extern "C" {
    fn openblas_get_config() -> *const std::ffi::c_char;
    fn openblas_get_num_threads() -> i32;
    fn openblas_get_num_procs() -> i32;
    fn openblas_get_parallel() -> i32;
}

pub fn backend_label() -> &'static str {
    #[cfg(target_os = "macos")]
    {
        "Accelerate"
    }
    #[cfg(target_os = "windows")]
    {
        "OpenBLAS"
    }
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    {
        "BLAS"
    }
}

/// Environment variable the one-thread child sets before the library initializes.
pub fn thread_env_var() -> &'static str {
    #[cfg(target_os = "macos")]
    {
        "VECLIB_MAXIMUM_THREADS"
    }
    #[cfg(not(target_os = "macos"))]
    {
        "OPENBLAS_NUM_THREADS"
    }
}

/// Library and pool description. On Windows this initializes OpenBLAS, so call it after the
/// thread-count variable is set for this process.
pub fn library_line() -> String {
    #[cfg(target_os = "macos")]
    {
        "Accelerate vecLib".to_string()
    }
    #[cfg(target_os = "windows")]
    unsafe {
        let config = std::ffi::CStr::from_ptr(openblas_get_config())
            .to_string_lossy()
            .into_owned();
        let parallel = match openblas_get_parallel() {
            0 => "sequential",
            1 => "threads",
            2 => "openmp",
            other => panic!("openblas_get_parallel returned {other}"),
        };
        format!(
            "{config} procs {} threads {} parallel {parallel}",
            openblas_get_num_procs(),
            openblas_get_num_threads(),
        )
    }
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    {
        "no BLAS linked".to_string()
    }
}

pub fn time_shape(input: &BlasInput<'_>, threads: &str) -> Vec<BlasOutcome> {
    vec![time_f64(input, threads), time_f32(input, threads)]
}

fn time_f64(input: &BlasInput<'_>, threads: &str) -> BlasOutcome {
    let cells = input.taxa * input.sites;
    let mut matrix = vec![0.0f64; cells];
    let mut product = vec![0.0f64; input.taxa * input.taxa];
    let mut pack_samples = Vec::with_capacity(RUNS);
    let mut syrk_samples = Vec::with_capacity(RUNS);
    pack_f64(input.codes, &mut matrix);
    syrk_f64(input.taxa, input.sites, &matrix, &mut product);
    for _ in 0..RUNS {
        let pack_mark = Instant::now();
        pack_f64(input.codes, &mut matrix);
        pack_samples.push(elapsed(pack_mark));
        let syrk_mark = Instant::now();
        syrk_f64(input.taxa, input.sites, &matrix, &mut product);
        syrk_samples.push(elapsed(syrk_mark));
    }
    let error = upper_error_f64(&product, input.taxa, input.reference);
    BlasOutcome {
        kernel: "dsyrk".to_string(),
        threads: threads.to_string(),
        pack: stats_of(&pack_samples),
        syrk: stats_of(&syrk_samples),
        max_abs: error.0,
        max_rel: error.1,
        bit_exact: error.2,
    }
}

fn time_f32(input: &BlasInput<'_>, threads: &str) -> BlasOutcome {
    let cells = input.taxa * input.sites;
    let mut matrix = vec![0.0f32; cells];
    let mut product = vec![0.0f32; input.taxa * input.taxa];
    let mut pack_samples = Vec::with_capacity(RUNS);
    let mut syrk_samples = Vec::with_capacity(RUNS);
    pack_f32(input.codes, &mut matrix);
    syrk_f32(input.taxa, input.sites, &matrix, &mut product);
    for _ in 0..RUNS {
        let pack_mark = Instant::now();
        pack_f32(input.codes, &mut matrix);
        pack_samples.push(elapsed(pack_mark));
        let syrk_mark = Instant::now();
        syrk_f32(input.taxa, input.sites, &matrix, &mut product);
        syrk_samples.push(elapsed(syrk_mark));
    }
    let error = upper_error_f32(&product, input.taxa, input.reference);
    BlasOutcome {
        kernel: "ssyrk".to_string(),
        threads: threads.to_string(),
        pack: stats_of(&pack_samples),
        syrk: stats_of(&syrk_samples),
        max_abs: error.0,
        max_rel: error.1,
        bit_exact: error.2,
    }
}

fn pack_f64(codes: &[u8], matrix: &mut [f64]) {
    for (index, &code) in codes.iter().enumerate() {
        matrix[index] = if code == MISSING_CODE { 0.0 } else { code as f64 };
    }
}

fn pack_f32(codes: &[u8], matrix: &mut [f32]) {
    for (index, &code) in codes.iter().enumerate() {
        matrix[index] = if code == MISSING_CODE { 0.0 } else { code as f32 };
    }
}

fn syrk_f64(taxa: usize, sites: usize, matrix: &[f64], product: &mut [f64]) {
    #[cfg(any(target_os = "macos", target_os = "windows"))]
    unsafe {
        cblas_dsyrk(
            CBLAS_COL_MAJOR,
            CBLAS_UPPER,
            CBLAS_NO_TRANS,
            taxa as i32,
            sites as i32,
            1.0,
            matrix.as_ptr(),
            taxa as i32,
            0.0,
            product.as_mut_ptr(),
            taxa as i32,
        );
    }
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    {
        let _ = (taxa, sites, matrix, product);
        panic!("cblas_dsyrk is linked on macOS (Accelerate) and Windows (OpenBLAS)");
    }
}

fn syrk_f32(taxa: usize, sites: usize, matrix: &[f32], product: &mut [f32]) {
    #[cfg(any(target_os = "macos", target_os = "windows"))]
    unsafe {
        cblas_ssyrk(
            CBLAS_COL_MAJOR,
            CBLAS_UPPER,
            CBLAS_NO_TRANS,
            taxa as i32,
            sites as i32,
            1.0,
            matrix.as_ptr(),
            taxa as i32,
            0.0,
            product.as_mut_ptr(),
            taxa as i32,
        );
    }
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    {
        let _ = (taxa, sites, matrix, product);
        panic!("cblas_ssyrk is linked on macOS (Accelerate) and Windows (OpenBLAS)");
    }
}

fn upper_error_f64(product: &[f64], taxa: usize, reference: &[f64]) -> (f64, f64, bool) {
    let mut actual = vec![0.0; pair_count(taxa)];
    let mut index = 0;
    for row in 0..taxa {
        for column in row..taxa {
            actual[index] = product[column * taxa + row];
            index += 1;
        }
    }
    let error = compare_triangles(&actual, reference);
    (error.max_abs, error.max_rel, error.bit_exact)
}

fn upper_error_f32(product: &[f32], taxa: usize, reference: &[f64]) -> (f64, f64, bool) {
    let mut max_abs = 0.0;
    let mut max_rel = 0.0;
    let mut bit_exact = true;
    let mut index = 0;
    for row in 0..taxa {
        for column in row..taxa {
            let value = product[column * taxa + row] as f64;
            let expected = reference[index];
            let diff = (value - expected).abs();
            if diff > max_abs {
                max_abs = diff;
            }
            let relative = diff / expected.abs().max(RELATIVE_FLOOR);
            if relative > max_rel {
                max_rel = relative;
            }
            if value != expected {
                bit_exact = false;
            }
            index += 1;
        }
    }
    (max_abs, max_rel, bit_exact)
}

fn elapsed(started: Instant) -> f64 {
    started.elapsed().as_secs_f64()
}

pub fn write_lines(outcomes: &[(usize, usize, BlasOutcome)]) -> String {
    let mut text = String::new();
    for (taxa, sites, outcome) in outcomes {
        let pack = &outcome.pack;
        let syrk = &outcome.syrk;
        text.push_str(&format!(
            "{taxa} {sites} {} {} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9} {:.16e} {:.16e} {}\n",
            outcome.kernel,
            outcome.threads,
            pack.min_seconds,
            pack.median_seconds,
            pack.samples[0],
            pack.samples[1],
            pack.samples[2],
            syrk.min_seconds,
            syrk.median_seconds,
            syrk.samples[0],
            syrk.samples[1],
            syrk.samples[2],
            outcome.max_abs,
            outcome.max_rel,
            if outcome.bit_exact { 1 } else { 0 },
        ));
    }
    text
}

pub fn parse_lines(text: &str) -> Vec<(usize, usize, BlasOutcome)> {
    text.lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(fields.len(), 17, "blas line: {line}");
            let samples = |offset: usize| -> Stats {
                stats_of(&[
                    fields[offset + 2].parse().unwrap(),
                    fields[offset + 3].parse().unwrap(),
                    fields[offset + 4].parse().unwrap(),
                ])
            };
            (
                fields[0].parse().unwrap(),
                fields[1].parse().unwrap(),
                BlasOutcome {
                    kernel: fields[2].to_string(),
                    threads: fields[3].to_string(),
                    pack: samples(4),
                    syrk: samples(9),
                    max_abs: fields[14].parse().unwrap(),
                    max_rel: fields[15].parse().unwrap(),
                    bit_exact: fields[16] == "1",
                },
            )
        })
        .collect()
}

#[cfg(all(test, any(target_os = "macos", target_os = "windows")))]
mod tests {
    use super::*;
    use gpu_gram::{dosage_gram, triangle_index};

    #[test]
    fn hand_gram_matches_both_precisions() {
        let codes = [2u8, 0, 1, 2];
        let reference = dosage_gram(&codes, 2, 2, 2);
        assert_eq!(reference[triangle_index(0, 0, 2)], 5.0);
        let input = BlasInput { taxa: 2, sites: 2, codes: &codes, reference: &reference };
        for outcome in time_shape(&input, "test") {
            assert!(outcome.bit_exact, "{} max_abs {}", outcome.kernel, outcome.max_abs);
        }
    }
}
