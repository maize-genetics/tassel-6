# TASSEL 6

Genotype analysis for plant genetics, as one pure-Kotlin core compiled to the JVM, WebAssembly, and
native. TASSEL 6 is the Kotlin Multiplatform successor to
[TASSEL 5](https://github.com/maize-genetics/tassel).

An analysis is written once, in `commonMain`, and the same code runs in a command-line tool, a
desktop application, and a browser tab with no server. Every ported analysis is tested against
TASSEL 5's own output, byte for byte where the format allows it.

> **Status: alpha (`6.0.0-alpha.1`).** Genotype summary, site and taxa filters, Centered IBS kinship,
> and PCA are ported, and this release reads bgzipped VCF and one region of a CSI- or TBI-indexed VCF
> on the command line, the desktop app, and in the browser. It is not a replacement for TASSEL 5 yet,
> and the [known limits](docs/poc/benchmarks.md#known-limits) are written down rather than glossed over.

## What works today

| | CLI | Desktop | Browser |
|---|---|---|---|
| Read HapMap, plain or gzipped | yes | yes | yes |
| Read VCF, the `GT` field only, plain or bgzipped | yes | yes | yes |
| Read one region of an indexed VCF (`.csi`, then `.tbi`) | yes (`--region`) | yes | yes |
| `summary`, `filter-sites`, `filter-taxa` | yes | yes | yes |
| `kinship` (Centered IBS) | yes | yes | yes |
| `pca` (scores, eigenvalues, eigenvectors) | yes | yes | yes |
| Write TASSEL 5-readable output | yes | yes | yes |

The desktop and browser apps are the same Compose Multiplatform code, and both compute locally. The
browser does not send data anywhere. A bgzipped VCF with a `.csi` or `.tbi` opens on a region: the
desktop app finds the index beside the file, and the browser takes the VCF and the index together.
With no index, the whole file is read. The CLI does the same with `--region`, and looks beside the
VCF for `.csi`, then `.tbi`.

## Quick start

TASSEL 6 needs JDK 25; Gradle provisions it if you do not have it.

```bash
git clone https://github.com/maize-genetics/tassel-6.git
cd tassel-6

# Command line
./gradlew :tassel-cli:installDist
tassel-cli/build/install/tassel/bin/tassel --help

# Desktop
./gradlew :tassel-app:run

# Browser
./gradlew :tassel-app:wasmJsBrowserDevelopmentRun
```

A worked example on the tutorial dataset in `testdata/tutorial/`:

```bash
T=tassel-cli/build/install/tassel/bin/tassel
$T filter-sites --input testdata/tutorial/mdp_genotype.hmp.txt --min-maf 0.05 --output filtered
$T kinship      --input filtered.hmp.txt --output kinship
$T pca          --input filtered.hmp.txt --output pca
```

Every subcommand and every one of its options is generated from an operation spec, so `--help` is
always current and the CLI, the UI's parameter forms, and the future language bindings cannot drift
apart.

## Parity with TASSEL 5

`testdata/golden/` holds TASSEL 5 v5.2.98's own output for every ported analysis, produced by
[`scripts/generate-goldens.sh`](scripts/generate-goldens.sh) and never edited by hand. The test
suite compares against it on every target, so the JVM, Node-hosted WebAssembly, and native builds
must all agree with TASSEL 5 and with each other.

Most comparisons are exact text. The documented exceptions are in
[`testdata/golden/MANIFEST`](testdata/golden/MANIFEST): PCA eigenvectors are compared up to sign,
and its eigenvalues without regard to row order, because TASSEL 5 gets them from an SVD that does
not sort them.

See [ADR-0009](docs/adr/0009-strict-golden-parity.md) for why the goldens are treated as an oracle
rather than as fixtures.

## Performance

Measured on an Apple M2 Max against TASSEL 5 v5.2.98 at the same heap; full tables, method, and
limits in [`docs/poc/benchmarks.md`](docs/poc/benchmarks.md).

- **Startup:** 0.18 s against TASSEL 5's 5.34 s. On tutorial-sized data this is the whole
  difference.
- **PCA on 1,000 taxa x 100,000 sites:** 57 s and 887 MiB, against 274 s and 7.7 GiB, because
  analyses stream site blocks instead of materializing matrices.
- **Centered IBS on the same set:** the first port took 19.3 s against TASSEL 5's 2.8 s of work. The
  kernel in this release is a blocked product across the cores the JVM dispatcher has, and takes 3.2 s
  on its own. The pipeline table in the benchmarks is still the earlier single-threaded measurement.

Reproduce all of it with `scripts/bench/run-benchmarks.sh`.

## Repository layout

| Path | Contents |
|------|----------|
| `tassel-core/` | Pure Kotlin model, I/O, and analyses. Targets `jvm`, `wasmJs`, `macosArm64`, and `linuxX64`. |
| `tassel-cli/` | Clikt CLI generated from the operation specs. |
| `tassel-app/` | Compose Multiplatform UI: shared code plus `desktopMain` and `wasmJsMain`. |
| `testdata/` | Tutorial inputs, TASSEL 5 golden outputs, and the htslib oracle for compressed and indexed VCF. |
| `docs/` | ADRs, the porting ledger, handoffs, and PoC notes. |
| `scripts/` | Golden generation, benchmarks, and the TASSEL 5 reference checkout. |

`commonMain` is pure Kotlin: no `java.*`, no htsjdk, no JVM-only dependency. Platform code lives in
`jvmMain`, `wasmJsMain`, or `nativeMain` behind `expect`/`actual`. This is the constraint that makes
the browser build possible, and CI enforces it.

## Building and testing

```bash
./gradlew :tassel-core:jvmTest   # fast inner loop
./gradlew :tassel-core:allTests  # JVM, WebAssembly on Node, and native
./gradlew check                  # everything CI runs
```

`macosArm64Test` needs a Mac with full Xcode; the Command Line Tools alone cannot link. On such a
machine, use `./gradlew check -x macosArm64Test`. `linuxX64Test` is the native gate in CI.

## Documentation

- [Benchmarks](docs/poc/benchmarks.md) and the [demo walkthrough](docs/poc/demo.md)
- [Changelog](CHANGELOG.md)
- [Architecture decision records](docs/adr/README.md)
- [Porting ledger](docs/porting-ledger.md): what has been ported from TASSEL 5, and from where
- [AGENTS.md](AGENTS.md): the rules this repository is developed under, by people and by agents

## Contributing

Branches follow git flow: `feature/<issue>-<slug>` off `develop`, and feature pull requests target
`develop` ([ADR-0010](docs/adr/0010-git-flow-and-versioning.md)). Commits and pull-request titles
follow [Conventional Commits](https://www.conventionalcommits.org); CI checks the title. Read
[AGENTS.md](AGENTS.md) before your first change: it is short, and it is the same brief given to the
AI agents that write much of this code.

## License

Apache License 2.0. See [LICENSE](LICENSE).
