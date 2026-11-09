# TASSEL 6 PoC benchmarks

How the TASSEL 6 proof of concept performs against TASSEL 5 on the same inputs and the same
analyses. Reproduce every JVM number in this document with:

```bash
scripts/bench/run-benchmarks.sh --runs 3
```

The script downloads the pinned TASSEL 5 standalone release, builds the TASSEL 6 CLI, generates the
synthetic inputs with `tassel dev simulate`, and prints the table below. See
[`scripts/bench/run-benchmarks.sh`](../../scripts/bench/run-benchmarks.sh) for the exact pipelines
and subcommands each row runs.

## What was measured

| | |
|---|---|
| Host | Apple M2 Max, 12 cores, 64 GiB, macOS 27.0.1 (26A434) |
| JDK | Eclipse Temurin 25.0.4.1+1, `-Xmx8g` for both tools |
| TASSEL 5 | v5.2.98 standalone, `-maxThreads 4` |
| TASSEL 6 | 6.0.0-dev at commit `68e0325` |
| Runs | Median of 3, after one warm-up run that is discarded |

Both tools get the same heap so neither wins on memory configuration alone. TASSEL 5 is pinned to
four threads for the same reason `scripts/generate-goldens.sh` pins it: it fixes Centered IBS at
1,000-site blocks, so results do not depend on the host's core count. The pipeline table below was
measured when TASSEL 6 was single threaded. The kinship kernel section later in this document is
the blocked product, which uses every core the JVM dispatcher has.

Datasets are the tutorial HapMap (281 taxa x 3,093 sites) and a synthetic 1,000 taxa x 100,000
sites HapMap from `tassel dev simulate --seed 0`. The 1,000 x 1,000,000 set the roadmap also asks
for is behind `--million`; it was not run for this document because TASSEL 5's PCA already needs
7.7 GiB at 100,000 sites and does not fit the 8 GiB budget ten times over.

## Results

Ratio above 1 means TASSEL 6 is faster. Peak is the peak resident set of the whole process.

| Dataset | Analysis | TASSEL 5 time | TASSEL 6 time | Ratio | TASSEL 5 peak | TASSEL 6 peak |
|---------|----------|--------------:|--------------:|------:|--------------:|--------------:|
| none | startup only | 5.34 s | 0.18 s | 29.67x | 105 MiB | 83 MiB |
| tutorial | load | 5.48 s | 0.25 s | 21.92x | 150 MiB | 105 MiB |
| tutorial | summary | 5.54 s | 0.30 s | 18.47x | 159 MiB | 130 MiB |
| tutorial | filter-sites | 5.52 s | 0.26 s | 21.23x | 154 MiB | 109 MiB |
| tutorial | filter-taxa | 5.53 s | 0.29 s | 19.07x | 153 MiB | 109 MiB |
| tutorial | kinship | 5.50 s | 0.33 s | 16.67x | 156 MiB | 125 MiB |
| tutorial | pca | 6.19 s | 0.52 s | 11.90x | 235 MiB | 123 MiB |
| synthetic-100k | load | 7.84 s | 1.72 s | 4.56x | 708 MiB | 820 MiB |
| synthetic-100k | summary | 8.32 s | 2.16 s | 3.85x | 864 MiB | 914 MiB |
| synthetic-100k | filter-sites | 8.00 s | 2.13 s | 3.76x | 721 MiB | 913 MiB |
| synthetic-100k | filter-taxa | 7.92 s | 2.10 s | 3.77x | 708 MiB | 906 MiB |
| synthetic-100k | kinship | 8.13 s | 19.51 s | 0.42x | 713 MiB | 811 MiB |
| synthetic-100k | pca | 274.14 s | 57.13 s | 4.80x | 7,721 MiB | 887 MiB |

### Read the tutorial rows carefully

The "startup only" row runs TASSEL 5's pipeline with no plugin and TASSEL 6 with `--version`.
Neither does any work. TASSEL 5 still takes 5.34 s, so **every tutorial row is about 97% startup**,
and the 12x-to-22x ratios there say almost nothing about the analyses. Subtracting the startup row
gives the actual work:

| Analysis | TASSEL 5 work | TASSEL 6 work |
|----------|--------------:|--------------:|
| tutorial load | 0.14 s | 0.07 s |
| tutorial summary | 0.20 s | 0.12 s |
| tutorial kinship | 0.16 s | 0.15 s |
| tutorial pca | 0.85 s | 0.34 s |

On tutorial-sized data the two tools do comparable work and TASSEL 5's cost is the JVM plus its
plugin registry. That startup cost is not a rounding error for users: it is the difference between
a tool that feels interactive and one that does not, and it is why the 100,000-site rows are the
honest comparison of the analyses themselves.

### Parity

For `load`, `filter-sites`, and `filter-taxa` on the tutorial data, the two tools' output files are
**byte-for-byte identical**. Kinship from the default kernel agrees with
`golden/kinship/centered_ibs.txt` within `abs 8e-6`, and `Centered_IBS.SumPk` within `abs 1e-4`
([ADR-0020](../adr/0020-kinship-numerics.md)). The retained TASSEL 5-order kernel is still
byte-identical, including the `SumPk` text. That split is what `CenteredIbsTest` enforces.

## Kinship kernel baseline

`centeredIbs` on its own, after one discarded warm-up, median of 3. The file is read before the
clock starts. TASSEL 5 is the whole kinship pipeline, so its columns include JVM startup (about
5.3 s on this host; see the startup row above) and the HapMap read. Reproduce with:

```bash
scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset tutorial --runs 3
scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset synthetic-100k --runs 3
scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset synthetic-5k --heap 16g --runs 3
scripts/bench/run-benchmarks.sh --only kinship-kernel --dataset synthetic-1m --heap 16g --million --runs 3
```

`synthetic-5k` is 5,000 taxa x 100,000 sites and `synthetic-1m` is 1,000 x 1,000,000. Both need
`--heap 16g`. A smaller heap skips them, unless `--dataset` names one, in which case the script
stops. The 1,000,000-site HapMap is about 2 GB of text. Neither of those two sizes was run for the
table below.

| | |
|---|---|
| Host | Apple M2 Max, 12 cores, 64 GiB, Darwin 27.0.0 arm64 |
| TASSEL 6 JDK | Eclipse Temurin 25.0.4.1+1, `-Xmx8g` (Gradle toolchain; the CLI is JDK 25) |
| TASSEL 5 JDK | Eclipse Temurin 21.0.7, `-Xmx8g` |
| TASSEL 6 | `tassel dev kinship-bench`, median of 3 after one warm-up |
| TASSEL 5 | v5.2.98 pipeline, median of 3 after one warm-up |

| Dataset | T6 kernel median | T6 peak | T5 1 thread | T5 4 threads | T5 8 threads |
|---------|-----------------:|--------:|------------:|-------------:|-------------:|
| tutorial (281 x 3,093) | 0.025 s | 146 MiB | 5.55 s | 5.52 s | 5.52 s |
| synthetic-100k (1,000 x 100,000) | 3.207 s | 845 MiB | 14.55 s | 8.11 s | 7.11 s |

On 1,000 x 100,000 the kernel is 3.2 s. The PoC kernel on the same shape was 19.3 s of work, and
TASSEL 5's 4-thread pipeline was 8.13 s of wall time, about 2.8 s of work after startup. Tutorial
rows do not separate thread counts: TASSEL 5's startup is the whole measurement.

The kernel line for the 100,000-site run was
`taxa=1000 sites=100000 panelColumns=524 panelSites=524 min=3.193361 median=3.206955`.
The tutorial run was
`taxa=281 sites=3093 panelColumns=1865 panelSites=100 min=0.021417 median=0.024996`.

On the tutorial file the written `SumPk` is `918.729785675474`, against the golden
`918.7297856754799` (absolute difference 5.9e-12). The largest absolute difference in the matrix
is 6.0e-7, inside the `abs 8e-6` contract.

## VCF reading

Whole-file VCF reads at a size both tools can load with `--heap 8g`, plus htslib on that file and on a
bgzipped wheat VCF larger than 2 GiB. Reproduce with:

```bash
scripts/bench/run-benchmarks.sh --only vcf --runs 3
```

The script writes the VCFs under `build/bench/work/` with
[`scripts/bench/generate-wheat-vcf.sh`](../../scripts/bench/generate-wheat-vcf.sh) and reuses them
when the taxa, site count, and seed match. Contig lengths are IWGSC RefSeq v2.1 (the same 21
chromosomes as the htslib fixtures, `chr1A` through `chr7D`). Every chromosome has sites, and
positions run from 1 through that chromosome's length. Genotypes are biallelic `A/C` calls from
Park-Miller seed 81. Nothing under `build/bench/` is committed.

| | |
|---|---|
| Host | Apple M2 Max, 12 cores, 64 GiB, macOS 27.0.1, Darwin 27.0.0 arm64 |
| JDK | TASSEL 6: 25.0.4.1. TASSEL 5: 21.0.7. Both `-Xmx8g` |
| TASSEL 5 | v5.2.98, `-maxThreads 4` |
| TASSEL 6 | 6.0.0-dev on `feature/90-indexed-vcf-region` (parent `1a12c3a`) |
| htslib | 1.25 (`bgzip`, `bcftools`) |
| Runs | Median of 3, after one warm-up run that is discarded |

The loadable file is 1,000 taxa by 1,000,000 sites. Both tools finish it at `--heap 8g`. Plain text
is 4,033,841,662 bytes (3.76 GiB); `bgzip -c -@ 1` is 524,156,980 bytes (0.49 GiB).

The region file is 1,000 taxa by 5,110,342 sites, streamed into `bgzip -c -@ 1` and indexed with
`bcftools index -c`. The site count comes from a 100,000-site sample of the same generator, which
compressed to 52,527,890 bytes (525.28 bytes per site), so this count would land near 2.5 GiB. The
file is 2,674,378,477 bytes (2.49 GiB), past 2 GiB. `tabix` is not used: sites on `chr3B` pass 2^29,
which a TBI cannot store.

TASSEL 6's whole-file rows are `tassel --quiet filter-sites` with no filters, which reads the VCF
and writes HapMap. The `.vcf.gz` row goes through the parallel BGZF reader. TASSEL 5 is `-vcf` and
`-exportType Hapmap` on the plain file and on the bgzipped file. Those rows include JVM startup and
the HapMap write.

TASSEL 6 decompression is `tassel dev bgzf-bench`: inflate only, one warm-up discarded, JVM startup
excluded. `--workers 1` is the single-worker row. Omitting `--workers` uses the processor count, 12
on this host. htslib decompression is `bgzip -d -@ N -o /dev/null` at 1 thread and at 4, process
startup included. Parsing is `bcftools view -Ou -o /dev/null`. htslib region rows are
`bcftools view -H -r` on the large file for `chr3B:1-1000000` (1 Mbp), `chr3B:1-50000000` (50 Mbp),
and `chr3B` (the whole contig, 851,934,019 bp), process startup included. TASSEL 6 region rows are
`tassel dev vcf-region-bench` on the same queries. The header and index are opened before the clock
starts, one warm-up is discarded, and JVM startup is excluded.

Inflate at 12 workers is 1.29 s, against 10.40 s at one worker. The whole-file `.vcf.gz` read is
20.93 s, next to 21.26 s for the plain file, so that read is no longer waiting on inflate. A 1 Mbp
region is 0.03 s once the file is open, against 0.22 s for `bcftools` including process startup.
The 50 Mbp region is 1.03 s, and all of `chr3B` is 16.87 s, against 0.81 s and 10.45 s for
`bcftools`. The whole-contig row builds an in-memory table; `bcftools` streams the lines and
discards them.

| Tool | Read | Input | Median | Peak |
|------|------|-------|-------:|-----:|
| TASSEL 6 | whole-file read | loadable .vcf | 21.26 s | 3455 MiB |
| TASSEL 6 | whole-file read | loadable .vcf.gz | 20.93 s | 3831 MiB |
| TASSEL 6 | decompress, 1 worker | loadable .vcf.gz | 10.40 s | 694 MiB |
| TASSEL 6 | decompress, 12 workers | loadable .vcf.gz | 1.29 s | 710 MiB |
| TASSEL 5 | whole-file read | loadable .vcf | 41.67 s | 5698 MiB |
| TASSEL 5 | whole-file read | loadable .vcf.gz | 42.25 s | 4050 MiB |
| htslib | decompress, 1 thread | loadable .vcf.gz | 1.88 s | 6 MiB |
| htslib | decompress, 4 threads | loadable .vcf.gz | 0.57 s | 8 MiB |
| htslib | decompress and parse | loadable .vcf.gz | 17.86 s | 8 MiB |
| TASSEL 6 | region `chr3B:1-1000000` | region .vcf.gz | 0.03 s | 747 MiB |
| htslib | region `chr3B:1-1000000` | region .vcf.gz | 0.22 s | 108 MiB |
| TASSEL 6 | region `chr3B:1-50000000` | region .vcf.gz | 1.03 s | 1187 MiB |
| htslib | region `chr3B:1-50000000` | region .vcf.gz | 0.81 s | 107 MiB |
| TASSEL 6 | region `chr3B` | region .vcf.gz | 16.87 s | 3487 MiB |
| htslib | region `chr3B` | region .vcf.gz | 10.45 s | 107 MiB |

## Known limits

### Centered IBS at scale

The PoC kernel was single-threaded and took 19.3 s of work on 1,000 x 100,000 sites, against 2.8 s
for TASSEL 5 at 4 threads. The default kernel is now the blocked dosage product in
`distance.centeredIbs` (ADR-0020), parallel on the JVM and on Native, and the baseline above is the
number to beat. `distance.centeredIbsTassel5Order` keeps the bit-exact port for tests. wasmJs stays
one thread. JVM SIMD and a native BLAS were not faster paths we shipped; see
[kinship-jvm-fast-path.md](kinship-jvm-fast-path.md) and
[kinship-native-fast-path.md](kinship-native-fast-path.md).

### The eigensolver is the PoC's slowest primitive

Measured separately in #39 on the JVM after warm-up, over five runs, for the pure-Kotlin
`SymmetricEigenDecomposition`:

| Matrix | Time | Max residual \|Av - lambda v\| |
|--------|-----:|-------------------------------:|
| 1000 x 1000 | 1,277 ms minimum, 1,518 ms median | 1.4e-13 |
| 281 x 281 | 29 ms | 9.2e-14 |

PCA only ever decomposes a taxa-sized matrix, so 281 x 281 is the size the tutorial actually hits
and 1000 x 1000 is the size the synthetic set hits. Even at 1000 x 1000 the decomposition is under
3% of the 57 s PCA run, which is why the maintainer accepted these figures for the PoC. A LAPACK
`dsyevd` fast path is out of scope. The accuracy is what matters more: this is one pure-Kotlin
solver producing the same answer on all four targets, with no native dependency.

### Memory on the synthetic set

TASSEL 6 holds 820-914 MiB where TASSEL 5 holds 708-864 MiB for the three streaming analyses, so
TASSEL 6 is slightly heavier there. The reverse is true where it counts: TASSEL 5's PCA peaks at
7.7 GiB of the 8 GiB heap against TASSEL 6's 887 MiB, because TASSEL 5 materializes intermediate
matrices that the site-block API ([ADR-0012](../adr/0012-streaming-genotype-access.md)) lets
TASSEL 6 stream.

## Cost of the responsive-UI change

ADR-0019, written in [#56](https://github.com/maize-genetics/tassel-6/pull/56) and since accepted, makes
`ProgressReporter.report` a suspending function so that long-running core work yields to the event
loop, which is what keeps the browser drawing. Yielding is not free, so it was measured on the JVM
before and after:

| Path | Before | After |
|------|-------:|------:|
| Tutorial HapMap read | 22.1 ms | 22.7 ms |
| `centeredIbs` on the tutorial data | 58.6 ms | 59.8 ms |

Both are inside the run-to-run noise. The engine path, which performs 31 real yields on this input,
also lands at 59.8 ms, so the yields themselves cost roughly 0.6 ms in total.

The table above was measured on `develop` at `68e0325`, which does **not** include the ADR-0019
change ([#56](https://github.com/maize-genetics/tassel-6/pull/56) and
[#57](https://github.com/maize-genetics/tassel-6/pull/57) were open at the time). Given deltas this
small relative to the reported times, re-running after those merge should not move any row.

## The browser

The same `commonMain` code compiled to WebAssembly, on the tutorial data. Load and PCA were measured
on 2026-11-06 against the production bundle (`./gradlew :tassel-app:wasmJsBrowserDistribution`) at
`dd6e29f`. Kinship was measured again the same day, after the heatmap was painted in slices, against
a fresh production bundle served the same way: `http://127.0.0.1:8765`, Chrome on the same Mac, tab
in the foreground. That window's animation frames were 16.7 ms apart.

Wall time is from the click to the result appearing. **Longest frame gap** is the longest interval
between two `requestAnimationFrame` callbacks during that window, which is the number that decides
whether the page feels alive: anything past about 100 ms reads as a freeze. Each step was run twice
on a freshly loaded page. The two runs did not agree to a few milliseconds, so both are reported.

| Step | Wall time | Longest frame gap | Frame gaps over 50 ms |
|------|----------:|------------------:|----------------------:|
| Load tutorial dataset, run 1 | 229 ms | 66 ms | 1 |
| Load tutorial dataset, run 2 | 169 ms | 100 ms | 1 |
| kinship, run 1 | 632 ms | 65 ms | 1 |
| kinship, run 2 | 502 ms | 34 ms | 0 |
| pca, run 1 | 638 ms | 583 ms | 1 |
| pca, run 2 | 454 ms | 417 ms | 1 |

Loading paints before the grid. A compositor frame between the click and the genotype grid showed
"Opening mdp_genotype.hmp.txt" and an indeterminate progress bar, with the Load button disabled.
The other frames in that window were 16.7 ms apart except for the one longer gap in the table.
Four later loads of the same file, on pages that had already parsed it once, had longest gaps of
50 ms, 33 ms, 17 ms, and 33 ms.

Kinship was remeasured on 2026-11-06 in a foreground Chrome window whose idle frames were 16.7 ms
apart (the tab `document.visibilityState` was `visible`). Load and PCA in the table are still the
earlier numbers; PCA was not part of this remeasure. Before the heatmap change below, the same
window spent kinship in one 317 ms gap: `centeredIbs` was already reporting about every 100 sites
and those reports were reaching the screen (the frames during the matrix stayed within 33 ms), and
the 317 ms was the first paint of the 281 by 281 heatmap, one rect per taxon pair. The heatmap is
now painted 24 rows at a time into a bitmap, yielding between slices (`yieldToEventLoop`). The two
kinship rows above are that build. Wall time is from the Run click to the last frame gap over 22 ms,
which is when the heatmap had finished and the matrix was on screen. Both gaps are under 100 ms.

`yield()` on wasm resumes as a microtask, so the wasm UI reporter suspends with `delay(1)`
(`yieldToEventLoop`), a timer and therefore a macrotask. Desktop stays on `yield()` so JVM kinship
time does not pick up a sleep. The worker question is in
[wasm-worker-spike.md](wasm-worker-spike.md).

The table this remeasure replaced, taken before the running banner yielded after each report:

| Step | Wall time | Longest frame gap | Frame gaps over 50 ms |
|------|----------:|------------------:|----------------------:|
| Load tutorial dataset | 160 ms | 33 ms | 0 |
| summary | 162 ms | 33 ms | 0 |
| kinship | 905 ms | 583 ms | 2 |
| pca | 1,216 ms | 1,083 ms | 1 |

What was measured between progress reports, which is the stretch a yield can break up if the
browser actually paints there:

| Path | Reports | Longest stretch between reports |
|------|--------:|--------------------------------:|
| `centeredIbs`, JVM, warm | 31 | 3 ms (about 60 ms for the whole matrix) |
| `centeredIbs`, wasmJs on Node, warm | 31 | 22 ms (604 ms for the whole matrix) |
| PCA Gram pass, JVM, warm | 13 | 11 ms |
| PCA Gram pass, wasmJs on Node, warm | 13 | 17 ms |
| PCA eigensolver, JVM, warm | none | about 35 ms, inside the gap before the first eigenvector report |
| PCA eigensolver, wasmJs on Node, warm | none | 139 ms, in that same gap |

`centeredIbs` reports 31 times on this file, about every 100 sites. On wasmJs under Node those steps
are about 22 ms. In the foreground Chrome window those reports already became frames; the gap that
used to cover the whole kinship run was the heatmap drawing every taxon pair in one task, and the
sliced bitmap is what brought that gap under 100 ms. The block size of the retained TASSEL 5-order
kernel, not its reporting interval, still fixes that kernel's summation order. The default kernel's
panel size was not changed, so its `Double` grouping stays the one ADR-0020 fixes from the table
shape.

PCA's passes report every 256 sites. That splits the Gram pass on Node (17 ms). It does not split
the eigensolver: `SymmetricEigenDecomposition` has no progress report, and ADR-0019 does not give it
one. On a warm Wasm run under Node that stretch is 139 ms. In Chrome the whole PCA run is still one
gap of 417 ms or 583 ms, so the per-report yields are not reaching the screen there either.

For scale, the same two analyses on the same data through the JVM CLI take 0.15 s and 0.34 s of
work. The browser walls above are about 0.6 s for kinship and about 0.5 s for PCA, so kinship is
about 4x the JVM work and PCA is under 2x. That is the expected cost of the target, not a
regression. The extra PCA reports did not slow the JVM down: four warmed tutorial runs took 181
to 193 ms (median 188 ms) with 39 reports, against 229 ms for one warmed run when each pass
reported once.

![TASSEL 6 running Centered IBS on the tutorial data in a browser](../images/browser-kinship.png)

## What is not here

- **1,000 x 1,000,000 sites.** Behind `--million`. TASSEL 5's PCA will not fit in a reasonable heap
  at that size, so the row would compare TASSEL 6 against a failure.
- **Native targets.** `linuxX64` and `macosArm64` compile the same `commonMain` code. There is no
  native CLI yet to time (Phase 1). `linuxX64Test` does not run on a Mac, and linking
  `macosArm64Test` on this machine failed in `xcrun xcodebuild -version` (exit 72), so native
  timings are in [kinship-native-fast-path.md](kinship-native-fast-path.md).
- **5,000 taxa and 1,000,000 sites.** The kinship-kernel commands are above. They were not run for
  this table.
