# GPU Gram spike

Go/no-go for a shared GPU backend, after timing the same dense Float32 Gram on Rust `wgpu`,
Accelerate, and the CPU baselines. Epic [tassel-6#111](https://github.com/maize-genetics/tassel-6/issues/111),
stage 1b ([tassel-6#113](https://github.com/maize-genetics/tassel-6/issues/113)).

## Recommendation

**Use platform-specific WebGPU.** Keep [`spikes/gpu-gram/wgsl/gram_f32.wgsl`](../../spikes/gpu-gram/wgsl/gram_f32.wgsl)
as the shared kernel and drive it from each platform's own WebGPU. Rust stays the benchmark
harness in [`spikes/gpu-gram/rust`](../../spikes/gpu-gram/rust). It is not a TASSEL 6 dependency,
which is the same boundary as [ADR-0001](../adr/0001-kotlin-multiplatform-pure-core.md).

Three measurements decide that:

- Dosage products were bit-exact (max abs `0`) for both WGSL entry points on Apple M2 Max Metal,
  on all four shapes. PR #114 already recorded the same bit-exact dosage result in Chrome 156,
  Firefox 159, and Safari 27.0.1.
- At 5,000 x 20,000, tiled `wgpu` took 1.693 s and Accelerate `ssyrk` took 0.243 s. The GPU run
  is slower than the platform BLAS. The epic's bar for a native GPU backend is a kernel at least
  3x faster than the best CPU path at 5,000 taxa. A shared Rust backend would also need one wasm
  artifact that runs the kernel without a second module. The `wasm32-unknown-unknown` build is a
  1,447,390-byte module that imports the browser's WebGPU through `wasm-bindgen` and `web-sys`.
  It was not run, so a speed claim for Rust-in-the-browser is not measured.
- Tiled `wgpu` took 1.189 s at 1,000 x 100,000 and 1.693 s at 5,000 x 20,000. The cited Chrome
  times for the same shader and the same shapes are 0.51 s and 0.83 s. Rust is 2.3x and 2.0x
  those Chrome times, on the same GPU. It does not beat the browser's own WebGPU.

The same four dosage shapes were timed on a Windows desktop with a discrete GPU. Both WGSL
entry points were bit-exact there, and OpenBLAS `dsyrk` and `ssyrk` matched the f64 reference
with max abs `0`. Tiled `wgpu` at 5,000 x 20,000 took 1.060 s, compared with 1.693 s on the M2.
OpenBLAS `ssyrk` on that machine, 32 threads, took 0.446 s in the following run. The GPU run is
slower than that BLAS. The Windows section has the table.

**Stage 2 is a go, as the next spike.** Dosage was bit-exact here and in the three browsers.
Chrome's cited 0.51 s at 1,000 x 100,000 is 62x the wasm CPU median measured here (31.512 s),
past the epic's 5x line. The dense Float32 Gram is not the kernel to ship: on this Mac it loses
to Accelerate and to TASSEL 5's residual, and on the Windows host it loses to OpenBLAS `ssyrk`. Stage 2 is the compact kernels. File these, and do not
open them from the harness PR:

- Pure-Kotlin bit-plane popcount kernel on the JVM and on wasm
- WGSL bit-plane popcount kernel in `spikes/gpu-gram`
- WGSL port of the TASSEL 5 Endelman lookup
- Missing-data sweep at 0, 5, and 20 percent on those kernels
- PCA covariance Gram from integer count matrices

**ADR-0022 is not written.** The epic's gate requires a GPU integer kernel that is bit-identical
on every tested engine, and the 5x browser win. The 5x half is met. Dense f32 dosage sums are
not that integer kernel.

This is the first BLAS timing on this workload. The JVM and native kinship spikes
([kinship-jvm-fast-path.md](kinship-jvm-fast-path.md),
[kinship-native-fast-path.md](kinship-native-fast-path.md)) stopped before calling BLAS or the
Vector API. On this host `ssyrk` is the fastest dense Gram. On the Windows host, OpenBLAS
`ssyrk` is likewise the fastest dense Gram. Putting either behind the JVM Foreign
Function & Memory API is out of scope for this issue.

## What was timed

Host: Mac14,5 (Apple M2 Max), 12 logical CPUs (8 performance, 4 efficiency), 64 GiB,
Darwin 27.0.0 arm64. Rustc 1.87.0. `wgpu` 26.0.1 (`wgpu` 27.0.1 requires Rust 1.88, so the
harness stays on the 1.87.0 pin from the issue). The `wgpu` adapter was Apple M2 Max, Metal,
`IntegratedGpu`, not a CPU or virtual adapter. One-time device and pipeline setup was 0.444 s
and is outside the table. TASSEL 6 kinship used JDK 25.0.4.1. TASSEL 5 was v5.2.98 on JDK 21.0.7
with `-Xmx8g`.

The dosage Gram and the kinship runs are the same shapes and different draws. Dosages are the
Park-Miller generator in [`spikes/gpu-gram/README.md`](../../spikes/gpu-gram/README.md), checksums
`4544405675908146118`, `14324033804557704143`, `3888745768927589105`, and `6644050518704627923`.
`tassel dev kinship-bench` and TASSEL 5 read HapMap from `dev simulate --seed 0`. Pair-sites,
`taxa * (taxa + 1) / 2 * sites`, is how those rows share a column. Times in the table are rounded.
Pair-sites/s divides the unrounded median into that count and keeps 3 significant figures. TASSEL 6 kinship divides by
12, the logical CPU count `kinshipWorkerCount` reads. The spike f64 Gram is one thread. The
default Accelerate pool does not report its width, so that row has no per-thread rate.

TASSEL 5 seconds in the rate columns are the pipeline median minus the startup median measured
this session, 5.35 s. The HapMap read and the matrix export stay in that residual, so the rate
is slower than the kernel alone. Tutorial residuals are 0.15 s to 0.19 s, smaller than the
startup measurement, and those per-thread rates are not meaningful. At 1,000 x 10,000 the
residual is still a small fraction of a ~6 s wall time.

`wgpu` medians are end-to-end for f32 accumulation of dosage panels: staging, `writeBuffer`,
dispatch plus device poll, one readback, and the copy out. One warm-up is discarded and the
median of 3 uses the `kinship-bench` rank. Dispatch is most of the large tiled runs (0.913 s of
1.189 s at 1,000 x 100,000, and 1.579 s of 1.693 s at 5,000 x 20,000). The readback is about
1.5 ms. Accelerate medians are the `syrk` call. Packing the column-major matrix is extra: 0.037 s
for `dsyrk` and 0.009 s for `ssyrk` at 1,000 x 100,000 and at 5,000 x 20,000. Both precisions
matched the f64 reference with max abs `0`.

Browser rows are cited from PR #114, tiled linear-memory f32 dosage, setup excluded. That run
was this machine class on macOS 27.0.1. They were not re-timed here.

| Implementation | Size | Median s | Pair-sites/s | Pair-sites/s/thread | Error | Source |
|---|---:|---:|---:|---:|---|---|
| TASSEL 5, 1 thread | 281 x 3,093 | 5.54 wall | not meaningful | not meaningful | different draw | measured |
| TASSEL 5, 4 threads | 281 x 3,093 | 5.50 wall | not meaningful | not meaningful | different draw | measured |
| TASSEL 5, 8 threads | 281 x 3,093 | 5.51 wall | not meaningful | not meaningful | different draw | measured |
| TASSEL 6 kinship-bench | 281 x 3,093 | 0.024 | 5.05e9 | 4.21e8 | different draw | measured |
| Spike JVM f64 | 281 x 3,093 | 0.031 | 3.92e9 | 3.92e9 | reference | measured |
| Spike wasm f64 | 281 x 3,093 | 0.080 | 1.54e9 | 1.54e9 | reference | measured |
| Chrome 156 WebGPU | 281 x 3,093 | 0.008 | 1.53e10 | n/a | dosage bit-exact | cited, PR #114 |
| Firefox 159 WebGPU | 281 x 3,093 | 0.60 | 2.04e8 | n/a | dosage bit-exact | cited, PR #114 |
| Safari 27.0.1 WebGPU | 281 x 3,093 | 0.010 | 1.23e10 | n/a | dosage bit-exact | cited, PR #114 |
| Rust wgpu naive | 281 x 3,093 | 0.014 | 9.07e9 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 281 x 3,093 | 0.009 | 1.36e10 | n/a | max abs 0 | measured |
| Accelerate dsyrk, default | 281 x 3,093 | 0.000683 | 1.79e11 | not measured | max abs 0 | measured |
| Accelerate ssyrk, default | 281 x 3,093 | 0.000192 | 6.40e11 | not measured | max abs 0 | measured |
| Accelerate dsyrk, 1 thread | 281 x 3,093 | 0.000679 | 1.80e11 | 1.80e11 | max abs 0 | measured |
| Accelerate ssyrk, 1 thread | 281 x 3,093 | 0.000187 | 6.54e11 | 6.54e11 | max abs 0 | measured |
| TASSEL 5, 1 thread | 1,000 x 10,000 | 1.04 work (6.39 wall) | 4.81e9 | 4.81e9 | different draw | measured |
| TASSEL 5, 4 threads | 1,000 x 10,000 | 0.54 work (5.89 wall) | 9.27e9 | 2.32e9 | different draw | measured |
| TASSEL 5, 8 threads | 1,000 x 10,000 | 0.51 work (5.86 wall) | 9.81e9 | 1.23e9 | different draw | measured |
| TASSEL 6 kinship-bench | 1,000 x 10,000 | 0.373 | 1.34e10 | 1.12e9 | different draw | measured |
| Spike JVM f64 | 1,000 x 10,000 | 1.211 | 4.13e9 | 4.13e9 | reference | measured |
| Spike wasm f64 | 1,000 x 10,000 | 3.107 | 1.61e9 | 1.61e9 | reference | measured |
| Chrome 156 WebGPU | 1,000 x 10,000 | 0.054 | 9.27e10 | n/a | dosage bit-exact | cited, PR #114 |
| Firefox 159 WebGPU | 1,000 x 10,000 | 1.20 | 4.17e9 | n/a | dosage bit-exact | cited, PR #114 |
| Safari 27.0.1 WebGPU | 1,000 x 10,000 | 0.074 | 6.76e10 | n/a | dosage bit-exact | cited, PR #114 |
| Rust wgpu naive | 1,000 x 10,000 | 0.220 | 2.28e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 1,000 x 10,000 | 0.071 | 7.01e10 | n/a | max abs 0 | measured |
| Accelerate dsyrk, default | 1,000 x 10,000 | 0.017 | 2.99e11 | not measured | max abs 0 | measured |
| Accelerate ssyrk, default | 1,000 x 10,000 | 0.006 | 7.89e11 | not measured | max abs 0 | measured |
| Accelerate dsyrk, 1 thread | 1,000 x 10,000 | 0.031 | 1.64e11 | 1.64e11 | max abs 0 | measured |
| Accelerate ssyrk, 1 thread | 1,000 x 10,000 | 0.008 | 5.92e11 | 5.92e11 | max abs 0 | measured |
| TASSEL 5, 1 thread | 1,000 x 100,000 | 9.20 work (14.55 wall) | 5.44e9 | 5.44e9 | different draw | measured |
| TASSEL 5, 4 threads | 1,000 x 100,000 | 2.81 work (8.16 wall) | 1.78e10 | 4.45e9 | different draw | measured |
| TASSEL 5, 8 threads | 1,000 x 100,000 | 1.77 work (7.12 wall) | 2.83e10 | 3.54e9 | different draw | measured |
| TASSEL 6 kinship-bench | 1,000 x 100,000 | 3.265 | 1.53e10 | 1.28e9 | different draw | measured |
| Spike JVM f64 | 1,000 x 100,000 | 12.486 | 4.01e9 | 4.01e9 | reference | measured |
| Spike wasm f64 | 1,000 x 100,000 | 31.512 | 1.59e9 | 1.59e9 | reference | measured |
| Chrome 156 WebGPU | 1,000 x 100,000 | 0.51 | 9.81e10 | n/a | dosage bit-exact | cited, PR #114 |
| Firefox 159 WebGPU | 1,000 x 100,000 | 10.0 | 5.01e9 | n/a | dosage bit-exact | cited, PR #114 |
| Safari 27.0.1 WebGPU | 1,000 x 100,000 | 0.77 | 6.50e10 | n/a | dosage bit-exact | cited, PR #114 |
| Rust wgpu naive | 1,000 x 100,000 | 2.530 | 1.98e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 1,000 x 100,000 | 1.189 | 4.21e10 | n/a | max abs 0 | measured |
| Accelerate dsyrk, default | 1,000 x 100,000 | 0.173 | 2.89e11 | not measured | max abs 0 | measured |
| Accelerate ssyrk, default | 1,000 x 100,000 | 0.066 | 7.62e11 | not measured | max abs 0 | measured |
| Accelerate dsyrk, 1 thread | 1,000 x 100,000 | 0.313 | 1.60e11 | 1.60e11 | max abs 0 | measured |
| Accelerate ssyrk, 1 thread | 1,000 x 100,000 | 0.089 | 5.63e11 | 5.63e11 | max abs 0 | measured |
| TASSEL 5, 1 thread | 5,000 x 20,000 | 40.13 work (45.48 wall) | 6.23e9 | 6.23e9 | different draw | measured |
| TASSEL 5, 4 threads | 5,000 x 20,000 | 12.39 work (17.74 wall) | 2.02e10 | 5.05e9 | different draw | measured |
| TASSEL 5, 8 threads | 5,000 x 20,000 | 8.38 work (13.73 wall) | 2.98e10 | 3.73e9 | different draw | measured |
| TASSEL 6 kinship-bench | 5,000 x 20,000 | 14.281 | 1.75e10 | 1.46e9 | different draw | measured |
| Spike JVM f64 | 5,000 x 20,000 | 67.369 | 3.71e9 | 3.71e9 | reference | measured |
| Spike wasm f64 | 5,000 x 20,000 | 178.172 | 1.40e9 | 1.40e9 | reference | measured |
| Chrome 156 WebGPU | 5,000 x 20,000 | 0.83 | 3.01e11 | n/a | dosage bit-exact | cited, PR #114 |
| Firefox 159 WebGPU | 5,000 x 20,000 | 4.19 | 5.97e10 | n/a | dosage bit-exact | cited, PR #114 |
| Safari 27.0.1 WebGPU | 5,000 x 20,000 | 2.09 | 1.20e11 | n/a | dosage bit-exact | cited, PR #114 |
| Rust wgpu naive | 5,000 x 20,000 | 9.984 | 2.50e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 5,000 x 20,000 | 1.693 | 1.48e11 | n/a | max abs 0 | measured |
| Accelerate dsyrk, default | 5,000 x 20,000 | 0.795 | 3.15e11 | not measured | max abs 0 | measured |
| Accelerate ssyrk, default | 5,000 x 20,000 | 0.243 | 1.03e12 | not measured | max abs 0 | measured |
| Accelerate dsyrk, 1 thread | 5,000 x 20,000 | 1.594 | 1.57e11 | 1.57e11 | max abs 0 | measured |
| Accelerate ssyrk, 1 thread | 5,000 x 20,000 | 0.457 | 5.47e11 | 5.47e11 | max abs 0 | measured |

At 1,000 x 100,000, TASSEL 5's 1-thread residual is 5.44e9 pair-sites/s and TASSEL 6 kinship is
1.28e9 per thread. The residual still includes the read and the export, so the bitwise kernel is
at least that far ahead of the blocked dosage product. TASSEL 5's 8-thread residual, 1.77 s, is
also a shorter wall time than TASSEL 6's 3.265 s kernel.

Commands:

```bash
cd spikes/gpu-gram/rust && cargo test && cargo run --release
./gradlew :spikes:gpu-gram:cpuBench -Ptassel.spikes=true
./gradlew :spikes:gpu-gram:wasmJsNodeProductionRun -Ptassel.spikes=true
scripts/bench/run-benchmarks.sh --only kinship-kernel --runs 3
scripts/bench/run-benchmarks.sh --only startup --runs 3
scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset synthetic-100k --runs 3
```

The 1,000 x 10,000 and 5,000 x 20,000 kinship rows are the one-off commands in the spike README.
`synthetic-5k` was skipped because it needs `--heap 16g`.

## Windows discrete GPU

Host: MSI MS-7D37, Intel Core i9-14900K, 24 cores, 32 logical processors, 137,207,046,144 bytes
of memory, Windows 11 Home 10.0.26300 (26H2). Rustc 1.87.0 (`17067e9ac`, 2025-05-09). `wgpu`
26.0.1. The adapter was NVIDIA GeForce RTX 4060 Ti, Vulkan, driver 617.42, `DiscreteGpu`,
`fallback` false. Vendor `4318`, device `10243`. An Intel UHD Graphics 770 is also installed.
The harness asked for high performance and opened the 4060 Ti (16,380 MiB). One-time device and
pipeline setup was 0.709 s and is outside the table. The device kept the default limits the
kernel uses: 256 invocations per workgroup, a 128 MB storage binding, and 65,535 workgroups per
dimension. The adapter advertised higher limits.

Dosages are the Park-Miller draws in the spike README. All four checksums matched, and the first
16 codes were `2, 2, 3, 1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 0, 1, 1`. `wgpu` medians are the same
end-to-end f32 protocol as the Mac rows: staging, `writeBuffer`, dispatch plus device poll, one
readback, and the copy out. One warm-up is discarded and the median of 3 is rank 1 of the sorted
samples. Both shaders and all four shapes had max abs `0`.

Dispatch is most of the large tiled runs: 0.214 s of 0.329 s at 1,000 x 100,000, and 0.857 s of
1.060 s at 5,000 x 20,000. Readback (`mapAsync` plus the copy out) at 5,000 x 20,000 is 0.014 s.
The f64 reference in that same process is one sample on one thread: 0.024 s, 0.977 s, 9.823 s,
and 56.525 s. It is not a median, and it is the scalar loop rather than a BLAS kernel.

Tiled `wgpu` is faster here than the M2 tiled row at every shape: 0.00196 s, 0.0356 s, 0.329 s,
and 1.060 s, against 0.009 s, 0.071 s, 1.189 s, and 1.693 s. At 5,000 x 20,000, 1.060 s is 4.4
times the Mac Accelerate `ssyrk` median of 0.243 s. The cited Chrome time for that shape is
0.83 s on the M2, so the two GPU samples are different machines. Pair-sites/s below uses the
unrounded median and three significant figures, the same rule as the Mac table.

OpenBLAS 0.3.34 is the official x64 build, linked as `cblas_dsyrk` and `cblas_ssyrk`.
`openblas_get_config` reported `DYNAMIC_ARCH NO_AFFINITY Haswell MAX_THREADS=64`. The default
pool used 32 threads. The one-thread child set `OPENBLAS_NUM_THREADS=1` before the library
initialized and then reported 1 thread. Medians are the `syrk` call after one discarded warm-up.
Packing is extra: 0.060 s for `dsyrk` and 0.033 s for `ssyrk` at 1,000 x 100,000, and 0.072 s
and 0.039 s at 5,000 x 20,000. Both precisions matched the f64 reference with max abs `0` on all
four shapes. `cargo test` passed the six library tests and `hand_gram_matches_both_precisions`,
the 2 by 2 Gram at both precisions.

The `wgpu` rows are the release run described above. The OpenBLAS rows are the next release run.
That run also re-timed tiled `wgpu` at 0.00192 s, 0.0331 s, 0.329 s, and 1.085 s, still max abs
`0`. In that process, tiled `wgpu` at 5,000 x 20,000 is 1.085 s and `ssyrk` is 0.446 s. The three
`ssyrk` samples were 0.450 s, 0.446 s, and 0.392 s. Default-pool `ssyrk` is the fastest path at
every shape. Tiled `wgpu` is faster than default-pool `dsyrk` at 281 x 3,093 and at 5,000 x
20,000. Pair-sites/s/thread for the 32-thread rows divides the wall rate by 32.

| Implementation | Size | Median s | Pair-sites/s | Pair-sites/s/thread | Error | Source |
|---|---:|---:|---:|---:|---|---|
| Rust wgpu naive | 281 x 3,093 | 0.00326 | 3.76e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 281 x 3,093 | 0.00196 | 6.26e10 | n/a | max abs 0 | measured |
| OpenBLAS dsyrk, 32 threads | 281 x 3,093 | 0.00275 | 4.45e10 | 1.39e9 | max abs 0 | measured |
| OpenBLAS ssyrk, 32 threads | 281 x 3,093 | 0.00182 | 6.73e10 | 2.10e9 | max abs 0 | measured |
| OpenBLAS dsyrk, 1 thread | 281 x 3,093 | 0.00387 | 3.16e10 | 3.16e10 | max abs 0 | measured |
| OpenBLAS ssyrk, 1 thread | 281 x 3,093 | 0.00218 | 5.63e10 | 5.63e10 | max abs 0 | measured |
| Rust wgpu naive | 1,000 x 10,000 | 0.149 | 3.37e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 1,000 x 10,000 | 0.0356 | 1.40e11 | n/a | max abs 0 | measured |
| OpenBLAS dsyrk, 32 threads | 1,000 x 10,000 | 0.0198 | 2.52e11 | 7.89e9 | max abs 0 | measured |
| OpenBLAS ssyrk, 32 threads | 1,000 x 10,000 | 0.00987 | 5.07e11 | 1.58e10 | max abs 0 | measured |
| OpenBLAS dsyrk, 1 thread | 1,000 x 10,000 | 0.151 | 3.32e10 | 3.32e10 | max abs 0 | measured |
| OpenBLAS ssyrk, 1 thread | 1,000 x 10,000 | 0.0740 | 6.76e10 | 6.76e10 | max abs 0 | measured |
| Rust wgpu naive | 1,000 x 100,000 | 1.385 | 3.61e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 1,000 x 100,000 | 0.329 | 1.52e11 | n/a | max abs 0 | measured |
| OpenBLAS dsyrk, 32 threads | 1,000 x 100,000 | 0.190 | 2.63e11 | 8.21e9 | max abs 0 | measured |
| OpenBLAS ssyrk, 32 threads | 1,000 x 100,000 | 0.0994 | 5.04e11 | 1.57e10 | max abs 0 | measured |
| OpenBLAS dsyrk, 1 thread | 1,000 x 100,000 | 1.476 | 3.39e10 | 3.39e10 | max abs 0 | measured |
| OpenBLAS ssyrk, 1 thread | 1,000 x 100,000 | 0.734 | 6.82e10 | 6.82e10 | max abs 0 | measured |
| Rust wgpu naive | 5,000 x 20,000 | 3.197 | 7.82e10 | n/a | max abs 0 | measured |
| Rust wgpu tiled | 5,000 x 20,000 | 1.060 | 2.36e11 | n/a | max abs 0 | measured |
| OpenBLAS dsyrk, 32 threads | 5,000 x 20,000 | 1.502 | 1.66e11 | 5.20e9 | max abs 0 | measured |
| OpenBLAS ssyrk, 32 threads | 5,000 x 20,000 | 0.446 | 5.61e11 | 1.75e10 | max abs 0 | measured |
| OpenBLAS dsyrk, 1 thread | 5,000 x 20,000 | 6.418 | 3.90e10 | 3.90e10 | max abs 0 | measured |
| OpenBLAS ssyrk, 1 thread | 5,000 x 20,000 | 3.093 | 8.08e10 | 8.08e10 | max abs 0 | measured |

## Limitations

- Two hosts. On the Windows host, tiled `wgpu` at 5,000 x 20,000 is 1.060 s in the first run and
  1.085 s in the OpenBLAS run, against OpenBLAS `ssyrk` at 0.446 s. The discrete GPU is still
  slower than the best CPU path measured on that machine. OpenBLAS selected Haswell. A Linux
  host is still unmeasured.
- The Rust wasm module was built and not run. Its speed against the Kotlin/Wasm WebGPU path is
  not measured.
- Centered kinship was not reimplemented in Rust. The scaled Float kinship error against `8e-6`
  remains the PR #114 result (largest miss `7.749e-7`).
- TASSEL 5 residuals include the HapMap read and the `SqrMatrix` export. Kinship-bench draws are
  not the Park-Miller dosage matrices.
- The default vecLib pool width was not reported, so those rows have a wall rate and no
  per-thread rate. `VECLIB_MAXIMUM_THREADS=1` is the measured single-thread rate.
