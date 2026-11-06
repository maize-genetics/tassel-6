# JVM kinship fast path spike

Go/no-go for a faster Centered IBS kernel than the pure-Kotlin blocked product, via the Vector API
or OpenBLAS `dsyrk` through the Foreign Function & Memory API (tassel-6#68).

## Recommendation

**No-go for this change.** Do not add `jdk.incubator.vector` or OpenBLAS to the build. The pure-Kotlin
kernel is the default, and on the host below it is already close to TASSEL 5's 4-thread pipeline
once startup is set aside. Neither alternative was timed: both are architectural dependencies
(AGENTS.md constraint 2), and this spike stops at the seam and the cost rather than landing either
one off by default.

A later spike can replace `kinshipTriangleUpdate` in `jvmMain` and time the two options against the
numbers below. Ship one only if it beats the pure-Kotlin median by enough to pay for the flags,
the packaging, and the fallback when the native library is missing.

## What was timed

`tassel dev kinship-bench`, median of 3 after one warm-up, file read excluded. Host: Apple M2 Max,
12 cores, 64 GiB, Darwin 27.0.0 arm64. TASSEL 6 JDK: Eclipse Temurin 25.0.4.1+1, `-Xmx8g`.
Command: `scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset synthetic-100k --runs 3`.

| Kernel | 281 x 3,093 | 1,000 x 100,000 | Max abs vs golden |
|--------|------------:|----------------:|------------------:|
| Pure Kotlin `PureKotlinSyrk` | 0.025 s | 3.207 s | 6.0e-7 on the tutorial matrix |
| Vector API | not timed | not timed | |
| OpenBLAS `dsyrk` via FFM | not timed | not timed | |

5,000 x 100,000 was not run (it needs `--heap 16g` and a long TASSEL 5 sweep). TASSEL 5's own
4-thread pipeline on the 100,000-site file was 8.11 s of wall time, of which about 5.3 s is the
startup row in [benchmarks.md](benchmarks.md), so the work is about 2.8 s. The pure-Kotlin kernel
at 3.2 s is the gap a fast path would have to close.

`SumPk` on the tutorial file was `918.729785675474` against the golden `918.7297856754799`.

## Seam

`KinshipTriangleUpdate` in `tassel-core` is the function a fast path replaces. The default is
`PureKotlinSyrk`. `kinshipTriangleUpdate` is an `internal var` in the same module, so a `jvmMain`
type can assign a different implementation without touching `commonMain`. The panel is a row-major
`DoubleArray` with an explicit stride; the triangle is the upper triangle including the diagonal,
row by row. wasmJs, macosArm64, and linuxX64 keep `PureKotlinSyrk`.

## What each option would cost

**Vector API (`jdk.incubator.vector`).** Still incubating on JDK 25. Compile and run need
`--add-modules jdk.incubator.vector`, and the JVM warns. No native library and no per-platform
packaging. Kotlin calls it as Java from `jvmMain` only. The CLI start script would gain the same
flag `scripts/bench/run-benchmarks.sh` already passes for native access
(`--enable-native-access=ALL-UNNAMED` is a different flag and does not enable the incubator module).

**OpenBLAS `dsyrk` through FFM.** FFM is final in JDK 25 (ADR-0011), so no `--enable-preview`.
`--enable-native-access` is required; the benchmark script and the desktop and CLI JVM args already
pass it. The cost is the library: find or ship OpenBLAS on macOS, Linux, and Windows, copy the
panel into an `Arena` / `MemorySegment` or keep the panel off-heap, and decide how OpenBLAS's
threads share the machine with `Dispatchers.Default`. If the library is absent the call must fall
back to `PureKotlinSyrk` and the run must still succeed.

## Why not go

The remaining gap on 1,000 x 100,000 is a fraction of a second next to TASSEL 5's 4-thread work,
and an incubator module or a shipped native library is a permanent support cost. Revisit when a
branch times both options on this host and on 5,000 taxa, where the taxa-squared term is 25 times
larger and a tuned `dsyrk` is more likely to matter.
