// Dense Float32 Gram, upper triangle only.
// triangle_start / triangle_index match triangleIndex in Spec.kt.
// @@COMMON

struct Params {
    taxa: u32,
    columns: u32,
    stride: u32,
    pair_base: u32,
}

@group(0) @binding(0) var<storage, read> panel: array<f32>;
@group(0) @binding(1) var<storage, read_write> triangle: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;

fn triangle_start(i: u32, n: u32) -> u32 {
    return (i * (2u * n - i + 1u)) / 2u;
}

fn triangle_index(i: u32, j: u32, n: u32) -> u32 {
    return triangle_start(i, n) + (j - i);
}

fn row_of(pair: u32, n: u32) -> u32 {
    var lo: u32 = 0u;
    var hi: u32 = n;
    loop {
        if (lo + 1u >= hi) {
            break;
        }
        let mid = (lo + hi) / 2u;
        if (triangle_start(mid, n) <= pair) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

fn dot_rows(i: u32, j: u32) -> f32 {
    var sum = 0.0;
    let row = i * params.stride;
    let other = j * params.stride;
    let columns = params.columns;
    for (var column = 0u; column < columns; column = column + 1u) {
        sum = sum + panel[row + column] * panel[other + column];
    }
    return sum;
}

// @@NAIVE

@compute @workgroup_size(256)
fn gram_naive(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = params.taxa;
    let pairs = n * (n + 1u) / 2u;
    let pair = gid.x + params.pair_base;
    if (pair >= pairs) {
        return;
    }
    let i = row_of(pair, n);
    let j = i + (pair - triangle_start(i, n));
    let sum = dot_rows(i, j);
    triangle[pair] = triangle[pair] + sum;
}

// @@TILED

var<workgroup> tile_a: array<f32, 256>;
var<workgroup> tile_b: array<f32, 256>;

@compute @workgroup_size(16, 16)
fn gram_tiled(
    @builtin(workgroup_id) group_id: vec3<u32>,
    @builtin(local_invocation_id) local_id: vec3<u32>,
) {
    let n = params.taxa;
    let tiles = (n + 15u) / 16u;
    let tile_pair = group_id.x + params.pair_base;
    let tile_pairs = tiles * (tiles + 1u) / 2u;
    let tile_i = row_of(tile_pair, tiles);
    let tile_j = tile_i + (tile_pair - triangle_start(tile_i, tiles));
    let i = tile_i * 16u + local_id.y;
    let j = tile_j * 16u + local_id.x;
    var sum = 0.0;
    let chunks = (params.columns + 15u) / 16u;
    for (var chunk = 0u; chunk < chunks; chunk = chunk + 1u) {
        let k = chunk * 16u + local_id.x;
        let row_a = tile_i * 16u + local_id.y;
        let row_b = tile_j * 16u + local_id.y;
        var a = 0.0;
        var b = 0.0;
        if (row_a < n && k < params.columns) {
            a = panel[row_a * params.stride + k];
        }
        if (row_b < n && k < params.columns) {
            b = panel[row_b * params.stride + k];
        }
        tile_a[local_id.y * 16u + local_id.x] = a;
        tile_b[local_id.y * 16u + local_id.x] = b;
        workgroupBarrier();
        if (tile_pair < tile_pairs && i < n && j < n && j >= i) {
            let a_row = local_id.y * 16u;
            let b_row = local_id.x * 16u;
            sum = sum + tile_a[a_row + 0u] * tile_b[b_row + 0u];
            sum = sum + tile_a[a_row + 1u] * tile_b[b_row + 1u];
            sum = sum + tile_a[a_row + 2u] * tile_b[b_row + 2u];
            sum = sum + tile_a[a_row + 3u] * tile_b[b_row + 3u];
            sum = sum + tile_a[a_row + 4u] * tile_b[b_row + 4u];
            sum = sum + tile_a[a_row + 5u] * tile_b[b_row + 5u];
            sum = sum + tile_a[a_row + 6u] * tile_b[b_row + 6u];
            sum = sum + tile_a[a_row + 7u] * tile_b[b_row + 7u];
            sum = sum + tile_a[a_row + 8u] * tile_b[b_row + 8u];
            sum = sum + tile_a[a_row + 9u] * tile_b[b_row + 9u];
            sum = sum + tile_a[a_row + 10u] * tile_b[b_row + 10u];
            sum = sum + tile_a[a_row + 11u] * tile_b[b_row + 11u];
            sum = sum + tile_a[a_row + 12u] * tile_b[b_row + 12u];
            sum = sum + tile_a[a_row + 13u] * tile_b[b_row + 13u];
            sum = sum + tile_a[a_row + 14u] * tile_b[b_row + 14u];
            sum = sum + tile_a[a_row + 15u] * tile_b[b_row + 15u];
        }
        workgroupBarrier();
    }
    if (tile_pair < tile_pairs && i < n && j < n && j >= i) {
        let index = triangle_index(i, j, n);
        triangle[index] = triangle[index] + sum;
    }
}
