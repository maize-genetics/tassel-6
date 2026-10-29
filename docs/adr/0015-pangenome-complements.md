# ADR-0015: Pangenome complements: consume grits/ropebwt3-phg outputs and add a viewer

- Status: Accepted
- Date: 2026-10-29

## Context

The lab's pangenome tools (grits, ropebwt3-phg, PHG v2) produce PS4G files,
haplotype VCFs (h.vcf), and imputed VCFs. TASSEL is where users run downstream
genetic analyses.

## Decision

- TASSEL 6 reads PS4G, h.vcf, and imputed VCF outputs and runs its analyses
  (summary, filters, kinship, PCA, association) on them.
- TASSEL 6 adds a PS4G/haplotype viewer in the Compose app.
- TASSEL 6 does not re-implement alignment, indexing, or imputation; those stay
  in the pangenome tools.

## Consequences

- Readers for these formats belong in `tassel-core` (pure Kotlin), so they work
  in the browser too.
- Format changes in grits/PHG need coordination; pin format versions in tests.

## Alternatives considered

- Embed the pangenome pipelines: duplicates Rust/C tools and grows scope.
- No pangenome support: misses the lab's main new data source.
