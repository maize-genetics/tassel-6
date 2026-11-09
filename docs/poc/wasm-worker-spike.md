# Wasm web worker spike

Go/no-go for running kinship in a Web Worker in the wasm build, so the page can paint on the main
thread while an analysis runs, and so the web app can later use more than one core (tassel-6#71).

## Recommendation

**Go later.** Yielding on a timer (`yieldToEventLoop` in the wasm UI) is the fix for the tutorial
freeze. A worker is a second copy of the wasm instance and a second copy of the genotypes, for
parallelism the browser does not have today. Build it when an analysis is too long for a macrotask
yield to keep the longest frame gap under 100 ms, or when shared-memory threads are actually on
the table. Do not change `Engine`. A worker engine would be a third implementation beside
`LocalEngine` and the planned `RemoteEngine` (ADR-0004).

The production bundle from `:tassel-app:wasmJsBrowserDistribution` is about 18 MB on disk: two wasm
files (8.2 MB and 3.6 MB) and `tasselApp.js` (520 KB). The tutorial HapMap next to it is 2,796,284
bytes. In an automated browser the page loaded, the tutorial grid appeared, and Centered IBS drew
its heatmap and the 281 by 281 matrix. Frame timing from that tab is not usable: `requestAnimationFrame`
sat near 300 ms, which is a background tab, not the 16.7 ms foreground cadence in
[benchmarks.md](benchmarks.md).

## How the dataset would get to the worker

Send the file bytes and parse them in the worker. A `GenotypeTable` is not transferable, and the
main thread does not need to hold genotypes if the worker owns the table and the handle. The
tutorial file `testdata/tutorial/mdp_genotype.hmp.txt` is 2.8 MB. Copying those bytes into the
worker is one `postMessage` of an `ArrayBuffer`, which the browser can transfer rather than clone.
A serialized genotype matrix of the same file is 281 x 3,093 bytes, under 1 MB, plus taxa and
positions. Parsing in the worker costs a second wasm compile of the reader and keeps the main
thread's heap free of the table. That is the better of the two options ADR-0004's handle model
already points at: the engine that runs the work is the one that stores the data.

A realistic file is larger. The 1,000 x 100,000 HapMap used in the kernel benchmark is on the
order of a gigabyte of text. Transferring it is a copy the main thread would rather not keep, which
is another reason to parse in the worker and drop the bytes.

## Does Kotlin/Wasm support a worker today

Not for free. `kotlinx-coroutines` does not create a worker or instantiate the wasm module inside
one. The production bundle from `:tassel-app:wasmJsBrowserDistribution` is a webpack output
(`tasselApp.js` plus the wasm file) loaded by Compose's `ComposeViewport` on the page's main
thread. A worker would need its own `WebAssembly.instantiate`, its own Kotlin entry that is not
`main` (the current `main` opens a viewport), and a `postMessage` protocol for `RunRequest`
parameters. `kotlinx.serialization` is already a dependency and can encode the request.
Progress has to cross back as `Progress` values and be applied on the UI side (ADR-0019).
Cancellation is a message the worker checks at the next `report`, same as `LocalEngine`'s
`ensureActive`.

The browser entry today passes `EmptyCoroutineContext` into `LocalEngine` because there is one
thread (`tassel-app` wasm `Main`). That is the line a worker engine would replace, not `Engine.run`.

## Would several workers work

One worker does not need shared memory. N workers over disjoint taxa bands would. The pure-Kotlin
product already splits the triangle into pair-balanced bands, so the split exists, but each band
reads the whole dosage panel. Shipping a panel to N workers is a copy per band unless they share
the wasm memory. Shared memory needs `SharedArrayBuffer`, which needs `Cross-Origin-Opener-Policy`
and `Cross-Origin-Embedder-Policy` on the response. The roadmap already lists those headers as a
hosting need for future wasm threads. They are not required for a single worker.

## What it would cost the page

A worker is a second wasm instance: another download of the module, another compile, and another
copy of the runtime. The first analysis pays that cold start on top of the parse. Memory is the
genotypes twice if the main thread also keeps the table for the grid; the useful design is that
the worker owns the table and the main thread keeps only what the grid asks for. The tutorial file
is small enough that the duplicate would not be the problem. A 1,000 x 100,000 table would be.

## Is it needed after the timer yield

The reports were already about every 100 sites. They did not become frames because `yield()` on
wasm is a microtask and Chrome runs `requestAnimationFrame` on a later turn. `delay(1)` in the wasm
`yieldToEventLoop` actual is a timer, which is that turn. With that yield, the tutorial kinship
frames stay inside the reports; the gap that remained was the heatmap drawing every taxon pair at
once. Painting that bitmap in slices brings the longest kinship gap to 65 ms or 34 ms
([benchmarks.md](benchmarks.md)), so a worker does not fix the tutorial freeze. It buys a place to
run longer analyses off the UI thread later. That is worth a follow-up issue, not a default in this
change.
