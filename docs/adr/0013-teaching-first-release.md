# ADR-0013: Teaching-first public release

- Status: Accepted
- Date: 2026-10-29

## Context

A complete TASSEL 5 replacement is a multi-year effort. A useful public
release is needed much sooner, both for users and to keep the project visible.
TASSEL is heavily used in courses and workshops, where installation is the most
common obstacle.

## Decision

The first public release targets teaching: a zero-install web app (plus the
desktop app) covering data loading, summaries, filters, kinship, PCA, and
association basics on teaching-scale data, with bundled tutorial datasets.

## Consequences

- Phase ordering favors the web MVP, UI polish, and tutorial content over the
  long tail of research analyses.
- Browser-only compute is enough for the first release (ADR-0004).
- Research-scale features (server tier, out-of-core storage) come later.

## Alternatives considered

- Research-first release: higher bar, longer time to first release.
- CLI-first release: reaches existing power users but misses the teaching audience.
