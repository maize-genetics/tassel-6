# Native kinship fast path spike

Go/no-go for calling a platform BLAS from the Kotlin/Native kinship kernel: Accelerate on
`macosArm64`, OpenBLAS on `linuxX64` (tassel-6#69). The R and Python bindings will run on these
targets (ADR-0007), so this is the numeric risk the roadmap names for Phase 1.

## Recommendation

**No-go.** Do not add a cinterop `.def` or an OpenBLAS link dependency. This machine compiled
`compileKotlinMacosArm64` and `compileKotlinLinuxX64`, and it could not time either native test
binary, so there is no measurement that a BLAS call is faster than `PureKotlinSyrk`. The same
`KinshipTriangleUpdate` seam the JVM spike uses is the one a native `actual` would replace. wasmJs
and the JVM keep the pure-Kotlin product.

## What was timed

Nothing on a native target.

| Target | Pure Kotlin | BLAS | Why |
|--------|-------------|------|-----|
| macosArm64 | not timed | Accelerate not called | `linkDebugTestMacosArm64` failed: `/usr/bin/xcrun xcodebuild -version` exited 72. Kotlin sources compiled. Full Xcode is not available here. |
| linuxX64 | not timed | OpenBLAS not called | `linuxX64Test` does not run on a Mac host. The Kotlin sources compiled. |

JVM pure-Kotlin, same kernel, same host (Apple M2 Max, 12 cores, 64 GiB), from
[benchmarks.md](benchmarks.md): 0.025 s on the tutorial file and 3.207 s on 1,000 x 100,000.
That is the number a native port has to beat, and it is also the Phase 1 comparison this spike
was meant to make. It is not made here. Maximum absolute difference of the JVM kernel from
`testdata/golden/kinship/centered_ibs.txt` on the tutorial file is 6.0e-7.

## What a go would cost

**Accelerate on macOS.** `cblas_dsyrk` from the Accelerate framework, linked with a cinterop
`.def`. The library is part of the OS, so there is nothing to ship. Linking the test binary
needs a Mac with full Xcode (`xcodebuild`), which is the same requirement `macosArm64Test` already
has. `macosArm64Test` is advisory in CI.

**OpenBLAS on Linux.** `cblas_dsyrk` present at link time and at run time. The CI image that runs
`linuxX64Test` would need the library and headers installed. A missing library must fall back to
`PureKotlinSyrk`.

**R and Python (ADR-0007).** A shared-library dependency shows up in the Phase 4 packages. rTASSEL
already ships prebuilt libraries, so this is not automatic grounds to refuse, but it is another
artifact to build for each platform. `mingwX64` is not a configured target.

## Why not go

There is no native timing, and the JVM kernel is already within a fraction of a second of TASSEL 5's
4-thread work on 1,000 x 100,000. Adding Accelerate or OpenBLAS before a host can show a speedup
would put a linker requirement on every native build for a gain that has not been measured.
