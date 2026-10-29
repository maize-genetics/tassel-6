# ADR-0002: Operation specs as the single source of truth

- Status: Accepted (DSL details are frozen in the `operation-specs` task)
- Date: 2026-10-29

## Context

In TASSEL 5, `PluginParameter` definitions drive the GUI dialogs, the
`run_pipeline.pl` flags, and generated R code (`GenerateRCode`). That reuse was
one of TASSEL 5's strengths. TASSEL 6 has more front ends: CLI, Compose UI, MCP
tools, a C API for R/Python, docs, and a legacy-flag translator. Hand-writing
each would drift.

## Decision

Every user-facing analysis is declared once as a typed operation spec in
`tassel-core`: id, input/output types, typed parameters with defaults, ranges,
help text, and the TASSEL 5 legacy plugin and parameter names. Front ends are
generated or driven from the spec registry.

```kotlin
val FilterSites = operation<GenotypeTable, GenotypeTable>(
    id = "filter-sites",
    legacyPlugin = "FilterSiteBuilderPlugin",
) {
    val minMaf by double(default = 0.0, range = 0.0..0.5, legacy = "siteMinAlleleFreq",
        help = "Minimum minor allele frequency")
    val minCount by int(default = 0, legacy = "siteMinCount",
        help = "Minimum number of non-missing taxa")
    run { input, progress -> filterSites(input, minMaf(), minCount(), progress) }
}
```

## Consequences

- Adding an analysis means writing the kernel plus one spec; CLI, UI form, MCP
  tool, and docs follow.
- The legacy translator and rTASSEL 1.0 compatibility come from the `legacy`
  names rather than separate mapping tables.
- The spec DSL is a frozen contract once the `operation-specs` task lands;
  changes need an ADR.
- Specs must stay serializable so the remote engine and MCP can describe them.

## Alternatives considered

- Annotation processing (KSP) over plain functions: less explicit, harder on
  Wasm/Native, and weaker for ranges and legacy names.
- Separate hand-written CLI and UI: simple at first, guaranteed drift later.
