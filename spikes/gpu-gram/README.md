# GPU Gram spike

A measurement of `K = X Xᵀ` for a dense Float32 genotype panel, driven from Kotlin/Wasm
through WebGPU. It is not part of the shipped program. Nothing here is wired into
`tassel-core`, `tassel-cli`, or `tassel-app`.

Issue: maize-genetics/tassel-6#112. The Rust harness in #113 can copy the generator below
without reading the Kotlin.

## Run

The module is compiled only when the build is asked for it:

```bash
./gradlew check
./gradlew check -Ptassel.spikes=true
./gradlew :spikes:gpu-gram:serveGpuGram -Ptassel.spikes=true
```

`check` without the property does not see this module. With the property, `check` runs the
JVM tests and the wasm Node tests (generator, f64 Gram, and the Centered IBS cross-check).
It does not launch a browser. `wasmJsBrowserTest` is disabled.

`serveGpuGram` prints a localhost URL (port 8765) and writes each posted JSON report under
`spikes/gpu-gram/build/reports/gpu-gram/`. That directory is gitignored. Open the URL in
Chrome, Firefox, and Safari. A finished report has `"complete": true` and is also written as
`chrome.json`, `firefox.json`, or `safari.json`. `?quick=1` runs only the 281 by 3,093 shape
and one sample, for a smoke test.

CPU timings, without a browser:

```bash
./gradlew :spikes:gpu-gram:cpuBench -Ptassel.spikes=true
./gradlew :spikes:gpu-gram:wasmJsNodeProductionRun -Ptassel.spikes=true
```

The JVM report is `build/reports/gpu-gram/jvm-cpu.json`. The Node task sets
`TASSEL_GPU_GRAM_OUT` to that same directory's `node-cpu.json`. A manual Node launch writes
there when its working directory is this module, or wherever `TASSEL_GPU_GRAM_OUT` points.

## Generator

Park-Miller, the same recurrence as `scripts/bench/generate-wheat-vcf.sh`. No floating point.

```text
modulus    = 2147483647
multiplier = 48271
seed       = 81
state      = seed
next():
    state = (state * 48271) % 2147483647   # 64-bit integer product
    return state
```

`seed` stays in `1 .. 2147483646`. For each site, then each taxon:

```text
mafPercent = 5 + (next() % 46)          # integer in 5..50
if next() % 100 < 5:  code = 3          # missing
else:
    code = [next() % 100 < mafPercent] + [next() % 100 < mafPercent]
                                        # 0, 1, or 2 copies
```

Storage is site-major (all taxa of site 0, then site 1, ...). The checksum is a wrapping
unsigned 64-bit integer, starting at 0:

```text
checksum = checksum * 131 + code
```

| Taxa | Sites | Checksum |
|------|------:|----------|
| 281 | 3,093 | 4544405675908146118 |
| 1,000 | 10,000 | 14324033804557704143 |
| 1,000 | 100,000 | 3888745768927589105 |
| 5,000 | 20,000 | 6644050518704627923 |

The first 16 codes of any of those shapes (the first site does not depend on later taxa) are
`2, 2, 3, 1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 0, 1, 1`.

## Gram

Missing codes contribute 0. The product is an upper triangle packed row by row, diagonal
included: `(0,0), (0,1), ... (0,n-1), (1,1), ...`. Index `(i, j)` with `j >= i` is
`i * (2n - i + 1) / 2 + (j - i)`.

The f64 reference adds panels of 1,024 columns. Each pair uses one f64 accumulator inside a
panel, then one add into the triangle. Raw `0/1/2` sums in these shapes are at most 80,000,
so they are exact in f32 and f64; panel width does not change them. Centered values are not
exact, and the panel width is part of that reference.

Centered columns follow `centeredIbs`: a biallelic site contributes its major allele,
`copies - 2p`, with `p` the Float division of the allele count by the gamete count, and a
missing genotype contributes 0. `SumPk` is `2 * sum(p * (1 - p))`. Dosages are encoded
`0 = C/C`, `1 = A/C`, `2 = A/A`, missing `= N/N`. The scaled Float kinship
(`unscaled / SumPk`, rounded to Float) matches public `centeredIbs` within `abs 8e-6` on
281 by 3,093. That bound is ADR-0020's bound on the scaled matrix, not on the unscaled Gram.

The GPU uploads the same columns as f32. `gram_f32.wgsl` has two entry points, `gram_naive`
and `gram_tiled` (16 by 16 workgroup, unrolled inner tile). Both write the upper triangle
only. Each dispatch is one panel and is finished with `onSubmittedWorkDone` before the next,
and stays under the default WebGPU limits (256 invocations, 128 MB storage bindings, 65,535
workgroups). Raised limits are requested and recorded; the kernel does not require them.

Two ways move a panel from Kotlin to the GPU. Per-element writes a JavaScript `Float32Array`
one value at a time. Linear memory copies the panel into a `kotlin.wasm.unsafe` allocation
(`wasmMemory`, 8-byte aligned) and calls `writeBuffer` once. Kotlin arrays are on the WasmGC
heap, so the linear-memory path is still a copy.

`f32` accumulation adds panel products on the GPU and reads the triangle back once.
`f64-cpu` accumulation clears that buffer every panel, reads it back, and adds it on the CPU
in f64. Firefox resolves `mapAsync` and `onSubmittedWorkDone` on a 100 ms timer (Mozilla bug
1870699). Each panel waits for `onSubmittedWorkDone` before the next, so that wait is part of
the dispatch phase, and the per-panel readback makes `f64-cpu` dominated by the timer when a
shape has many panels.

Relative error is `abs(gpu - ref) / max(abs(ref), 1e-12)`. `within8e6` is the scaled Float
kinship against `8e-6`. A miss is a finding. The bound is not loosened. A large relative
error beside a small absolute error is a near-zero reference pair.

Break-even in the JSON is the set of shapes where the fastest linear-memory f32 dosage Gram
beats the wasm CPU median, excluding one-time adapter and pipeline setup. `none` is a valid
answer.

## Report fields

A browser report has `host` (machine model and OS, from the collector's `GET /host`), `userAgent`,
`setupSeconds`, `adapter` (including `isFallbackAdapter`
and requested versus granted limits), `fallback` (`forcedFailure`, `lostDevice`, and the CPU
self-test triangle `[5, 2, 4]`), and one entry per shape. Each GPU entry names `shader`,
`transfer`, `accumulation`, and `panel` (`dosage` or `centered`), then phase medians:
staging, `writeBuffer`, dispatch plus `onSubmittedWorkDone`, `mapAsync` (including the copy
into the mappable buffer), copy back, and CPU accumulation. Dosage entries set
`dosageBitExact`. Centered entries set `centeredMaxAbs`, `centeredMaxRel`, `scaledMaxAbs`,
and `within8e6`.
