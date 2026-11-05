# ADR-0019: Progress reports suspend, putting a yield point inside long-running core work

- Status: Proposed
- Date: 2026-11-05

## Context

ADR-0004 made `Engine` asynchronous and put the coroutine boundary there:
`Engine.run` and `Engine.release` are `suspend`, and everything beneath them is
plain Kotlin. `ProgressReporter.report` is not a suspend function, so once
control enters a reader or an analysis kernel there is no suspension point
until that function returns.

On the JVM and on Native this is invisible, because `LocalEngine` runs
operations on `Dispatchers.Default`, off the thread drawing the UI. The browser
has one thread. Parsing the 2.8 MB tutorial HapMap file (281 taxa x 3,093
sites) holds that thread for the whole parse, so Compose never gets a frame:
the loading banner never appears and the window is frozen until the grid shows
up. Running an analysis does the same thing. These were filed as two problems
(tassel-6#54 for reads, noticed in tassel-6#22 for analyses) but they have one
cause.

The place to yield already exists and is already instrumented. Every reader and
kernel reports progress at a site-block boundary or every N sites:

| Entry point | Reports at |
|-------------|------------|
| `io.HapMap.read`, `io.Vcf.read` | every 4,096 sites while parsing |
| `io.HapMap.write` | each site block |
| `analysis.genotypeSummary` | each site block, in each of its three passes |
| `analysis.ReferenceProbability.of` | each site block |
| `analysis.principalComponents` | each site block, Gram and eigenvector passes |
| `distance.centeredIbs` | each site block |
| `filter.filterSites`, `filter.filterTaxa` | each site block |

What is missing is not a place to yield but permission to yield there.

One of those intervals is wrong for the only file the PoC demo loads. The
readers report every 4,096 sites, and the tutorial file has 3,093, so a
tutorial read reports exactly once, after the parse has finished. A yield at
the report point would therefore fire after the freeze rather than during it.
The interval is a tuning number, not a contract, but it has to be part of this
decision or the decision does not produce a visible frame.

Moving the boundary down has a cost in Phase 4. ADR-0007 builds
`bindings/capi` as a Kotlin/Native shared library behind a C API, and **a C ABI
cannot call a Kotlin suspend function.** Today a binding can call
`distance.centeredIbs` directly, or bridge coroutines once at `Engine`, which
is exactly where ADR-0004 put the seam. That trade is why this is an ADR rather
than a line in whichever pull request needed it first; `ProgressReporter` also
belongs to the operation DSL frozen by ADR-0002.

## Decision

- `ProgressReporter.report` becomes a `suspend` function. The readers, the
  writers, and the five analysis kernels listed above become `suspend fun`.
  Nothing else about them changes: the loops, the block sizes, and the numbers
  they produce stay as they are.
- Reporting progress is the only suspension point added. Kernels do not call
  `yield()`, `delay()`, or `ensureActive()` themselves; they report, and
  whoever supplies the `ProgressReporter` decides what that costs.
- The readers report every 256 sites instead of every 4,096, so a
  tutorial-sized file reports about twelve times during the parse rather than
  once after it. Block-wise kernels keep the site block as their unit.
- `LocalEngine` supplies the policy. It already wraps the caller's reporter to
  call `ensureActive()`; it now also yields before forwarding the report, so
  the single-threaded browser gets a turn at the event loop roughly as often as
  the UI gets a progress number. Front ends that want no yielding (the CLI)
  pass a reporter straight through and pay nothing.
- The suspend boundary described by ADR-0004 still exists at `Engine`. This ADR
  does not move it; it adds a second, finer one inside the work so that the
  engine's boundary can actually be used on a single-threaded target. ADR-0004
  is refined, not superseded.
- Phase 4 bindings bridge coroutines **once**, in `bindings/capi` at the C
  boundary, rather than per analysis. The C API is already handle-based and
  spec-driven (ADR-0007), so it calls `Engine`, and `Engine` is where the
  bridge belongs.

## Consequences

- The browser can paint while a file is parsed and while an analysis runs, with
  no web worker and no second thread. This is what tassel-6#54 asks for and
  what the roadmap's "long computations yield in chunks" already assumed.
- Progress stops being merely computed and starts being visible. The reports
  the kernels already emit now reach the UI during the run instead of all at
  once after it.
- Cancellation gets better for free: a cancelled run stops at the next report,
  and the suspension point makes that a real cancellation rather than a thrown
  exception from a non-suspending callback.
- Coroutines reach into the numeric layer. `tassel-core`'s analysis functions
  can no longer be called from plain Kotlin without `runBlocking` or a
  coroutine scope. Tests in `commonTest` use `runTest`, which they already do
  for anything going through `Engine`.
- `bindings/capi` cannot expose `distance.centeredIbs` directly. It must go
  through `Engine`, which ADR-0007's handle-based design already implies. The
  cost is one blocking adapter in `bindings/capi`, written once.
- Yielding is not free. Yielding once per site block on a 1M-site file is a few
  hundred dispatches, which is negligible next to the arithmetic in a block,
  but yielding per site would not be. Block size stays the unit, and the pull
  request that implements this reports measured JVM parse and kinship times
  before and after.
- `.cursor/rules/common-main.mdc` says "report progress and check cancellation
  in long loops"; that is now one call rather than two.

## Alternatives considered

- **Move the work off the UI thread instead.** Already true on desktop, where
  `Dispatchers.Default` is a real thread pool. On wasm it needs a Web Worker,
  which means shipping bytes or a built `GenotypeTable` across the worker
  boundary; `kotlinx-coroutines` does not provide this on wasm. The roadmap
  defers web workers past the PoC and tassel-6#22 lists them as out of scope.
  Deferred, not rejected forever: a worker would still want these suspension
  points.
- **Offer both a suspend and a blocking entry point for every reader and
  kernel.** Keeps the C bindings trivial. Rejected because it doubles the
  public surface of every analysis and creates two code paths that have to stay
  in numeric parity, to avoid writing one blocking adapter in `bindings/capi`.
- **Yield through the `Sequence` returned by `GenotypeTable.siteBlocks`.** A
  `Sequence` cannot suspend; `Flow` can. Changing `siteBlocks` to return a
  `Flow` would change a contract frozen by ADR-0012 and would make every
  block-wise algorithm a flow collector. Rejected.
- **Accept the freeze and label it in the UI** (a wait cursor, a disabled
  button). Cheapest, and wrong for the demo: the browser target exists to show
  the same analysis running locally in a page, and a page that stops responding
  for the length of a real dataset's parse does not show that.
