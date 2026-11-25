# ADR-0027: Association results: columnar records, sinks, legacy reports, and provenance

- Status: Accepted
- Date: 2026-11-25

## Context

A GWAS produces one result per (trait, marker, test). At 1,000,000 markers and several traits that
is millions of rows. TASSEL 6's only table type, `table.ResultTable`, stores
`List<List<Any?>>` with boxed cells, which is fine for summaries and PCA tables and roughly
100 bytes per cell for a GWAS: a 1M-marker, 20-column, one-trait run would need about 2 GB before
anything is written.

The frozen contracts shape the answer. `Operation<I, O>` returns one value, and `Engine.run`
returns a `DataHandle` to an output the engine keeps in memory (ADR-0002, ADR-0004). There is no
streaming output type in the DSL. A `file(kind = FileKind.Output)` parameter exists, but
`commonMain` has no portable file system: browsers have none, and the CLI writes outputs itself.

TASSEL 5 reports GLM and MLM results as `TableReport`s with fixed column names, which users,
rTASSEL, and downstream scripts parse. Those layouts are a compatibility surface:

- GLM statistics: `Trait, Marker, Chr, Pos, marker_F, p, marker_Rsq, add_F, add_p, dom_F, dom_p,
  marker_df, marker_MS, error_df, error_MS, model_df, model_MS, minorObs` (`perm_p` after `p` with
  permutations; `addEffect, domEffect` appended with `appendAddDom`). GLM genotype effects:
  `Trait, Marker, Chr, Pos, Obs, Allele, Estimate`.
- MLM statistics: `Trait, Marker, Chr, Pos, df, F, p, add_effect, add_F, add_p, dom_effect, dom_F,
  dom_p, errordf, MarkerR2, Genetic Var, Residual Var, -2LnLikelihood`, with a `None` null-model
  row first per trait. MLM effects: `Trait, Marker, Locus, Site, Allele, Effect, Obs`. Compression:
  `Trait, # groups, Compression, -2LnLk, Var_genetic, Var_error`.
- Exported through `DoubleFormat` (five decimals, E-notation below 0.001 and for every negative
  number); written straight to files with `Double.toString`.

TASSEL 5 also keeps tables small with `maxP` (GLM drops NaN p when `maxP < 1`; MLM keeps NaN rows)
and streams to files with GLM's `saveToFile`, `siteFile`, and `alleleFile` and MLM's
`-mlmOutputFile`, which is the precedent for bounded-memory output.

## Decision

### Records

- New package `net.maizegenetics.tassel.association` defines `AssociationRecord` fields as a
  fixed, documented schema rather than free-form cells:
  - identity: trait, marker index, marker name, chromosome, position;
  - test: model (`glm`, `mlm`), test kind (`genotype` omnibus, `additive`, `dominance`, ...),
    marker encoding, effect allele, comparison allele(s);
  - counts: observations, unique taxa, allele count and frequency for the effect allele, minor
    allele count;
  - estimates: effect, standard error, confidence interval bounds;
  - statistics: statistic value, numerator df, denominator df, p, log10 p, marker R²;
  - variance components (MLM): genetic variance, residual variance, heritability, method, and
    whether the null-model estimate was reused (P3D) or re-estimated;
  - status: `ok`, `skipped` with a typed reason (monomorphic, below minimum class size, below MAF
    or MAC, insufficient df, not estimable, non-finite), or `failed` with a message.
- Effect-level detail (one row per genotype class or allele) is a second record type,
  `AssociationEffect`, keyed by trait and marker.

### Storage and flow

- Kernels write to an `AssociationSink` (`suspend fun accept(chunk)`), one chunk per site block,
  in deterministic (trait, marker) order regardless of worker count.
- The default sink is `AssociationTable`: primitive-array columns (`DoubleArray`, `IntArray`,
  interned strings), growing by chunk, about 8 bytes per numeric cell. It is the operation output,
  so `Engine`, `RunRequest`, and the DSL do not change.
- Output filters (`maxP`, `ok`-only, top-k per trait) are applied by the sink, so a filtered run
  stores only what it keeps. Summary counts (tested, skipped by reason, failed) are always kept in
  full.
- **Bounded-memory streaming to files** (plan decision Q4, chosen 2026-11-25). The `glm` and `mlm`
  specs declare TASSEL 5's `saveToFile`, `siteFile`, and `alleleFile` (and MLM's output base) as
  `file(kind = FileKind.Output)` parameters with those legacy names. When they are set, the kernel
  writes through a file-backed `AssociationSink` that uses kotlinx-io `SystemFileSystem`
  (ADR-0017) on the JVM, native, and Node targets, flushing per site block, and the returned
  `AssociationTable` keeps only summaries and `maxP` survivors. In the browser, where there is no
  file system, the parameters are refused with a clear message. The spec stays the single source
  for every front end (ADR-0002), and the CLI does not need a separate streaming path.

### Reports

- **TASSEL 5-compatible reports** are writers over `AssociationTable`: the GLM statistics and
  genotype-effect tables and the MLM statistics, effect, and compression tables, with TASSEL 5's
  column names, order, and number formatting, so existing scripts keep working. These writers are
  the parity surface for goldens (ADR-0029). Their columns never change; new information goes
  elsewhere. A column whose TASSEL 5 value is wrong (for example MLM `errordf`, which is not the
  df used for `p`) is listed in the departure table and either reproduced for compatibility or
  corrected. The maintainer chose to correct them (plan decision Q6, 2026-11-25): the compatible
  report prints the correct value, and the departure table names each affected column and why.
- **Modern reports** are TSV (and CSV) with the full record schema, explicit `log10_p`, effect
  alleles, and status. Columns are versioned by a `#schema` header line. Columnar formats (Arrow,
  Parquet) are a later, JVM-first writer behind the same table and need their own dependency ADR.
- **Front ends never recompute statistics.** The CLI, app, MCP, and bindings read the table.

### Provenance

Every `AssociationTable` carries an `AnalysisProvenance`: TASSEL 6 version, operation id, bound
parameter values, input names and dimensions (taxa, markers, observations), model specification
in canonical text, kinship source and scaling, random seed (permutations), and the counts of
excluded observations and markers by reason. Modern reports write it as `##key=value` lines, as
`DistanceMatrix` does; legacy reports stay byte-compatible and omit it.

## Consequences

- A 1M-marker, one-trait GLM table costs on the order of 150 MB without a p filter and almost
  nothing with one, instead of gigabytes.
- `ResultTable` is unchanged; `AssociationTable` can produce one for small displays.
- Plots (Manhattan, QQ) and multiple-testing corrections read the same columns.
- Streaming output (#209) adds a file-backed sink behind the same `AssociationSink` interface; the
  in-memory and streamed reports must be byte-identical.

## Alternatives considered

- **Reuse `ResultTable`.** No new type, and memory use rules it out at GWAS scale. Rejected.
- **Change `Operation` to emit a `Flow` of rows.** Natural for streaming, and it changes a frozen
  contract and every front end. Rejected for this epic; revisit with the out-of-core store.
- **Only TASSEL 5 layouts.** Zero design work, and no place for log p, alleles, status, or
  provenance. Rejected.
- **Only modern layouts.** Cleaner, and it breaks rTASSEL and user scripts. Rejected.
