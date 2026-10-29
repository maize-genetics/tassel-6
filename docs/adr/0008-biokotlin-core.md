# ADR-0008: Reuse BioKotlin through a multiplatform `biokotlin-core`

- Status: Accepted
- Date: 2026-10-29

## Context

BioKotlin already models genomic positions, ranges, sample gametes, and IUPAC
codes. It is JVM-only and depends on htsjdk and Guava, so it cannot be used in
`tassel-core`'s `commonMain`.

## Decision

- Carve a Kotlin Multiplatform `biokotlin-core` out of BioKotlin upstream,
  containing pure types such as `Position`, `PositionRange`, `SampleGamete`, and
  IUPAC data.
- `tassel-core` depends on `biokotlin-core` in `commonMain`.
- The full BioKotlin library is used only from JVM source sets and JVM-only
  modules (server, bridge).

## Consequences

- Shared vocabulary between TASSEL 6, PHG, and BioKotlin.
- Requires coordinated changes in the BioKotlin repo. Until `biokotlin-core` is
  published, `tassel-core` may hold minimal local equivalents, which must be
  replaced once it is available. Track this in the porting ledger.

## Alternatives considered

- Copy the types into TASSEL 6: fast, but forks the vocabulary.
- Depend on full BioKotlin: impossible in `commonMain`.
