# ADR-0004: Hybrid web compute: in-browser Wasm first, JVM server later

- Status: Accepted
- Date: 2026-10-29

## Context

Teaching and exploration use small data (for example, the 281-taxon x
3,093-site tutorial set) and benefit from zero-install browser use. Research
datasets (millions of sites) exceed browser memory and single-thread speed.

## Decision

- The web app runs analyses in the browser through the Wasm build of
  `tassel-core` for teaching-scale data.
- Analyses are invoked through an `Engine` interface. The PoC ships only
  `LocalEngine`. A later `tassel-server` (Ktor, JVM) provides `RemoteEngine`
  for large data, using the same operation specs over the wire.

## Consequences

- No server is required for the PoC or the teaching release.
- The `Engine` contract must be asynchronous, report progress, support
  cancellation, and pass data by handle so remote execution fits later.
- Browser limits (memory, single main thread) must be handled with block-wise
  algorithms and, where needed, web workers.

## Alternatives considered

- Browser-only: caps TASSEL 6 web at teaching scale permanently.
- Server-only web app: needs hosting from day one and loses offline teaching use.
