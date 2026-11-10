//! Park-Miller dosage generator and the f64 Gram reference for the Rust harness.
//!
//! The numbers match `spikes/gpu-gram` (issue #112) so a wgpu run can be compared with the
//! browser spike on the same bytes. This crate is a benchmark harness. It is not a TASSEL 6
//! dependency.

use std::path::Path;

pub const PANEL_COLUMNS: usize = 1024;
pub const PARK_MILLER_SEED: i32 = 81;
pub const PARK_MILLER_MULTIPLIER: i64 = 48_271;
pub const PARK_MILLER_MODULUS: i64 = 2_147_483_647;
pub const MISSING_PERCENT: i32 = 5;
pub const MAF_MIN_PERCENT: i32 = 5;
pub const MAF_MAX_PERCENT: i32 = 50;
pub const MISSING_CODE: u8 = 3;
pub const RELATIVE_FLOOR: f64 = 1.0e-12;

pub struct Shape {
    pub taxa: usize,
    pub sites: usize,
    pub checksum: u64,
}

pub const SHAPES: [Shape; 4] = [
    Shape { taxa: 281, sites: 3_093, checksum: 4_544_405_675_908_146_118 },
    Shape { taxa: 1_000, sites: 10_000, checksum: 14_324_033_804_557_704_143 },
    Shape { taxa: 1_000, sites: 100_000, checksum: 3_888_745_768_927_589_105 },
    Shape { taxa: 5_000, sites: 20_000, checksum: 6_644_050_518_704_627_923 },
];

/// First 16 site-major codes. The first site does not depend on how many taxa follow.
pub const PREVIEW: [u8; 16] = [2, 2, 3, 1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 0, 1, 1];

/// `ParkMiller(81).next()` in the Kotlin spike.
pub const FIRST_DRAW: i32 = 3_909_951;

pub struct ParkMiller {
    state: i32,
}

impl ParkMiller {
    pub fn new(seed: i32) -> Self {
        assert!(
            (1..PARK_MILLER_MODULUS as i32).contains(&seed),
            "seed must be in 1..2147483646: {seed}"
        );
        Self { state: seed }
    }

    pub fn next(&mut self) -> i32 {
        self.state = ((self.state as i64 * PARK_MILLER_MULTIPLIER) % PARK_MILLER_MODULUS) as i32;
        self.state
    }
}

/// Site-major dosage codes. Missing is [`MISSING_CODE`]. No floating point.
pub fn generate_dosages(taxa: usize, sites: usize, seed: i32) -> Vec<u8> {
    assert!(taxa > 0 && sites > 0, "taxa and sites must be positive: {taxa} x {sites}");
    let mut random = ParkMiller::new(seed);
    let mut codes = vec![0u8; sites * taxa];
    let mut index = 0;
    let maf_span = MAF_MAX_PERCENT - MAF_MIN_PERCENT + 1;
    for _site in 0..sites {
        let maf_percent = MAF_MIN_PERCENT + (random.next() % maf_span);
        for _taxon in 0..taxa {
            let code = if random.next() % 100 < MISSING_PERCENT {
                MISSING_CODE
            } else {
                let first = if random.next() % 100 < maf_percent { 1 } else { 0 };
                let second = if random.next() % 100 < maf_percent { 1 } else { 0 };
                first + second
            };
            codes[index] = code;
            index += 1;
        }
    }
    codes
}

/// Wrapping unsigned checksum: start at 0, then `checksum = checksum * 131 + code`.
pub fn dosage_checksum(codes: &[u8]) -> u64 {
    let mut checksum = 0u64;
    for &code in codes {
        checksum = checksum.wrapping_mul(131).wrapping_add(code as u64);
    }
    checksum
}

pub fn pair_count(taxa: usize) -> usize {
    taxa * (taxa + 1) / 2
}

pub fn triangle_index(row: usize, column: usize, taxa: usize) -> usize {
    assert!(column >= row, "column {column} is above row {row}");
    row * (2 * taxa - row + 1) / 2 + (column - row)
}

/// Adds one panel's `Z Zᵀ` into a packed upper triangle. One f64 accumulator per pair.
pub fn add_panel(panel: &[f64], taxa: usize, columns: usize, stride: usize, triangle: &mut [f64]) {
    if taxa == 0 || columns == 0 {
        return;
    }
    let mut index = 0;
    for first in 0..taxa {
        let row = first * stride;
        let mut second = first;
        let limit = taxa.saturating_sub(3);
        while second < limit {
            let other0 = second * stride;
            let other1 = other0 + stride;
            let other2 = other1 + stride;
            let other3 = other2 + stride;
            let mut sum0 = 0.0;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            for column in 0..columns {
                let value = panel[row + column];
                sum0 += value * panel[other0 + column];
                sum1 += value * panel[other1 + column];
                sum2 += value * panel[other2 + column];
                sum3 += value * panel[other3 + column];
            }
            triangle[index] += sum0;
            triangle[index + 1] += sum1;
            triangle[index + 2] += sum2;
            triangle[index + 3] += sum3;
            index += 4;
            second += 4;
        }
        while second < taxa {
            let other = second * stride;
            let mut sum = 0.0;
            for column in 0..columns {
                sum += panel[row + column] * panel[other + column];
            }
            triangle[index] += sum;
            index += 1;
            second += 1;
        }
    }
}

fn fill_dosage_panel(codes: &[u8], taxa: usize, site_start: usize, columns: usize, stride: usize, panel: &mut [f64]) {
    for taxon in 0..taxa {
        let row = taxon * stride;
        for column in 0..columns {
            let code = codes[(site_start + column) * taxa + taxon];
            panel[row + column] = if code == MISSING_CODE { 0.0 } else { code as f64 };
        }
    }
}

/// Raw dosage Gram. Missing codes contribute `0`. Panel width does not change these sums.
pub fn dosage_gram(codes: &[u8], taxa: usize, sites: usize, panel_columns: usize) -> Vec<f64> {
    let mut triangle = vec![0.0; pair_count(taxa)];
    let mut panel = vec![0.0; taxa * panel_columns];
    let mut site = 0;
    while site < sites {
        let columns = panel_columns.min(sites - site);
        fill_dosage_panel(codes, taxa, site, columns, panel_columns, &mut panel);
        add_panel(&panel, taxa, columns, panel_columns, &mut triangle);
        site += columns;
    }
    triangle
}

pub struct ValueError {
    pub max_abs: f64,
    pub max_rel: f64,
    pub bit_exact: bool,
}

pub fn compare_triangles(actual: &[f64], reference: &[f64]) -> ValueError {
    assert_eq!(actual.len(), reference.len(), "triangles differ");
    let mut max_abs = 0.0;
    let mut max_rel = 0.0;
    let mut bit_exact = true;
    for (&got, &expected) in actual.iter().zip(reference) {
        let diff = (got - expected).abs();
        if diff > max_abs {
            max_abs = diff;
        }
        let relative = diff / expected.abs().max(RELATIVE_FLOOR);
        if relative > max_rel {
            max_rel = relative;
        }
        if got != expected {
            bit_exact = false;
        }
    }
    ValueError { max_abs, max_rel, bit_exact }
}

/// `../wgsl/gram_f32.wgsl`, split the way `GpuGram.kt` builds the two modules.
pub fn load_shader_sources() -> (String, String) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../wgsl/gram_f32.wgsl");
    let source = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("reading {}: {error}", path.display()));
    split_shader(&source)
}

pub fn split_shader(source: &str) -> (String, String) {
    let naive_at = source.find("// @@NAIVE").expect("gram_f32.wgsl is missing // @@NAIVE");
    let tiled_at = source.find("// @@TILED").expect("gram_f32.wgsl is missing // @@TILED");
    assert!(tiled_at > naive_at, "tiled marker precedes naive marker");
    let common = &source[..naive_at];
    (
        format!("{common}{}", &source[naive_at..tiled_at]),
        format!("{common}{}", &source[tiled_at..]),
    )
}

pub fn shader_path() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../wgsl/gram_f32.wgsl")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn park_miller_first_draw() {
        assert_eq!(ParkMiller::new(PARK_MILLER_SEED).next(), FIRST_DRAW);
    }

    #[test]
    fn triangle_index_packs_the_upper_triangle() {
        assert_eq!(triangle_index(0, 0, 3), 0);
        assert_eq!(triangle_index(0, 1, 3), 1);
        assert_eq!(triangle_index(0, 2, 3), 2);
        assert_eq!(triangle_index(1, 1, 3), 3);
        assert_eq!(triangle_index(1, 2, 3), 4);
        assert_eq!(triangle_index(2, 2, 3), 5);
    }

    #[test]
    fn hand_dosage_gram() {
        let gram = dosage_gram(&[2, 0, 1, 2], 2, 2, 2);
        assert_eq!(gram[triangle_index(0, 0, 2)], 5.0);
        assert_eq!(gram[triangle_index(0, 1, 2)], 2.0);
        assert_eq!(gram[triangle_index(1, 1, 2)], 4.0);
    }

    #[test]
    fn dosage_gram_does_not_depend_on_panel_width() {
        let codes = generate_dosages(5, 17, PARK_MILLER_SEED);
        let narrow = dosage_gram(&codes, 5, 17, 1);
        let wide = dosage_gram(&codes, 5, 17, 8);
        assert_eq!(narrow, wide);
    }

    #[test]
    fn published_checksums_and_preview() {
        for shape in SHAPES {
            let codes = generate_dosages(shape.taxa, shape.sites, PARK_MILLER_SEED);
            assert_eq!(dosage_checksum(&codes), shape.checksum, "{}x{}", shape.taxa, shape.sites);
        }
        let tutorial = generate_dosages(281, 3_093, PARK_MILLER_SEED);
        assert_eq!(&tutorial[..16], &PREVIEW);
    }

    #[test]
    fn shader_file_has_both_entry_points() {
        let (naive, tiled) = load_shader_sources();
        assert!(naive.contains("fn gram_naive"));
        assert!(!naive.contains("fn gram_tiled"));
        assert!(tiled.contains("fn gram_tiled"));
        assert!(!tiled.contains("fn gram_naive"));
        assert!(naive.contains("struct Params"));
        assert!(tiled.contains("struct Params"));
    }
}
