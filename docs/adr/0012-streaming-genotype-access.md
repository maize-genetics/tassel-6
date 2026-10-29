# ADR-0012: Streaming site-block genotype access; out-of-core storage later

- Status: Accepted (interfaces frozen in the `core-model-io` task)
- Date: 2026-10-29

## Context

TASSEL 5 loads whole genotype tables into memory (`SuperByteMatrix`), which
limits dataset size and is impossible in the browser at research scale. The
long-term goal is out-of-core analysis of very large datasets.

## Decision

- `GenotypeTable` exposes genotypes as an iterable sequence of site blocks
  (contiguous sites x all taxa, or a taxa subset), in addition to random access.
- Analyses are written against site blocks so they accumulate results without
  needing the full matrix in memory.
- The PoC backs `GenotypeTable` with in-memory byte arrays. Chunked, memory-mapped,
  or file-backed storage is added later behind the same interface.
- Genotype byte encoding follows TASSEL 5's `NucleotideAlignmentConstants`.

## Consequences

- Kinship and PCA (Gram-matrix accumulation) and summaries fit naturally.
- Analyses that genuinely need random access can still use it but should say so.
- Blocks are also the natural unit for progress reporting, cancellation, and
  parallel work.

## Alternatives considered

- Whole-matrix API like TASSEL 5: simplest to port, blocks out-of-core later.
- Fully lazy per-genotype access: too slow without batching.
