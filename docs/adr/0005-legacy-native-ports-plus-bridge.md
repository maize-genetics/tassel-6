# ADR-0005: Legacy analyses: native ports plus a TASSEL 5 bridge

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 5 has well over a hundred plugins. Porting all of them before release
is not feasible, but users rely on many of them, including through rTASSEL.

## Decision

- High-use analyses are ported natively into `tassel-core`, with golden parity
  tests (ADR-0009). The porting ledger tracks the choice for each class.
- The long tail is reachable through `legacy/tassel5-bridge`, a JVM-only module
  that runs TASSEL 5 plugins in an isolated classloader and converts data at
  the boundary.
- The bridge is never a dependency of `tassel-core`.

## Consequences

- Isolation prevents clashes between TASSEL 5's dependencies (htsjdk 2.x, old
  Guava) and TASSEL 6's.
- Bridged analyses are available on the JVM desktop, CLI, and server only, not
  in the browser or Native builds. The UI and CLI must say so clearly.
- Each bridged analysis is a candidate for a later native port; usage data
  should drive the order.

## Alternatives considered

- Port everything natively first: delays any release by a long time.
- Depend on TASSEL 5 directly on the JVM: dependency conflicts and a JVM-only core.
- Drop the long tail: breaks existing user workflows.
