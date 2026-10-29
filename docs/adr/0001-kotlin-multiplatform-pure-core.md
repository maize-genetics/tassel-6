# ADR-0001: Kotlin Multiplatform with a pure-Kotlin common core

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 6 must ship a desktop app, a browser app, a CLI, and R/Python libraries
from one codebase. TASSEL 5 is JVM-only Java with htsjdk, Guava, and Swing
dependencies throughout. The earlier `tassel_wasm` prototype compiled TASSEL 5
with TeaVM and needed Wasm-specific forks and stubs for much of the dependency
graph, which does not scale.

## Decision

- Kotlin is the primary language.
- `tassel-core` is a Kotlin Multiplatform module whose `commonMain` has no
  JVM-only dependencies.
- Core targets: `jvm`, `wasmJs` (browser and Node), `macosArm64`, `linuxX64`.
- Platform-specific code (file systems, BLAS, threads) lives in platform source
  sets behind interfaces or `expect`/`actual`.

## Consequences

- Analyses are written once and tested on every target.
- TASSEL 5 code is ported, not linked: htsjdk and Guava are unavailable in
  `commonMain`, so formats such as VCF need small pure-Kotlin readers.
- JVM-only speedups (FFM + OpenBLAS, Vector API) are optional fast paths in
  `jvmMain`, never the only implementation.
- Windows Native and iOS targets can be added later without core changes.

## Alternatives considered

- TeaVM over TASSEL 5: proven possible but needs forks and stubs; no Native target.
- JVM core plus a separate JS/TS web core: two implementations to keep in parity.
- Rust core: strong Wasm and Native story, but loses the lab's Kotlin/BioKotlin
  base and the TASSEL 5 Java porting path.
