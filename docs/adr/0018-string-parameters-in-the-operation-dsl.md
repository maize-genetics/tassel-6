# ADR-0018: String and string-list parameters in the operation DSL

- Status: Accepted (amends ADR-0002; the DSL stays frozen, this adds two types)
- Date: 2026-11-04

## Context

ADR-0002 makes the operation spec the single source for the CLI, UI forms, MCP
tools, and the legacy translator, and the `operation-specs` task (#11) froze the
DSL with five parameter factories: `int`, `double`, `boolean`, `enum`, and
`file`. AGENTS.md constraint 4 says changing a frozen contract needs an ADR.

Three TASSEL 5 parameters cannot be declared with those five (#16, #42):
`startChr` and `endChr` on `FilterSiteBuilderPlugin` are chromosome names, which
are not integers (`1`, `chr1`, `Mt`), and `taxaList` on
`FilterTaxaBuilderPlugin` is an inline comma-separated list of taxon names.
`file(...)` is not a workaround for either: neither value is a path, and a front
end that saw `file` would offer a file chooser.

The consequences were visible in four places. `filterSitesOperation` threw for
the `POSITIONS` range type, so three of the ten golden parameter sets in #16
were reachable only by calling the `SiteFilter` and `TaxaFilter` kernels
directly. The generated CLI (#20) had no `--start-chr`, `--end-chr`, or
`--taxa-list`. The generated UI forms (#22) could not offer those fields. And
`tassel translate` (ADR-0006) resolves TASSEL 5 flags through the `legacy` names
on spec parameters, so it could not translate `-taxaList` either.

## Decision

Add two parameter types to the DSL, additively. Nothing about the existing five
changes.

- `OperationScope` gains `string(default, legacy, help)` / `string(legacy, help)`
  and `stringList(default, legacy, help)` / `stringList(legacy, help)`, in the
  shape every other factory already has: one overload with a default and one
  without, whose value is `null` unless supplied, plus `required()`.
- `ParameterDescriptor` gains `StringParameterDescriptor` and
  `StringListParameterDescriptor`, serialized with the `type` discriminators
  `string` and `string-list`. A string list's `default` is a JSON array.
- A string list is written as **one comma-separated value**, which is what
  TASSEL 5's `AbstractPlugin.convert` does for a `TaxaList` on the command line
  (`-taxaList W22,B73,CML247,33-16`). Entries are not trimmed, as in TASSEL 5;
  a blank entry is rejected by `Operation.bind` rather than passed on, because
  it can only come from a stray comma and would match nothing downstream.
- The generated CLI registers them as plain text options
  (`--start-chr=<text>`, `--taxa-list=<text,text>`), like every other generated
  option, so `Operation.bind` stays the single validator.

`filter-sites` therefore declares `startChr` and `endChr`, and `filter-taxa`
declares `taxaList`, each with its TASSEL 5 `legacy` name.

## Consequences

- The three golden parameter sets that only the kernels could reach
  (`filter-sites/position-range-5-6`, `filter-taxa/include-list`,
  `filter-taxa/exclude-list`) now run through the specs, and
  `filterSitesOperation` no longer throws for `POSITIONS`.
- `--start-chr`, `--end-chr`, and `--taxa-list` appear in the CLI with no
  hand-written code, and `tassel translate` can reach `-startChr`, `-endChr`,
  and `-taxaList` through `Operation.parameterByLegacyName`.
- `string` is the DSL's escape hatch, so a parameter that deserves a narrower
  type can be declared with it instead. Reviewers should push back on a `string`
  whose values are a fixed set (that is `enum`) or a path (that is `file`).
- A string list has no element-level validation: a spec cannot say "these must
  be taxon names in the table". The kernels still do that, as they did before.
- TASSEL 5 also accepts a `.txt` or `.json` file wherever it takes a `TaxaList`.
  That form is not ported; it needs file reading, which a spec parameter does
  not do, and no golden uses it.

## Alternatives considered

- A `chromosome` type and a taxa-list type that parse to `Chromosome` and
  `TaxaList`: more precise, but it puts a model type in the serialized
  descriptor that every front end would have to understand, and `Chromosome` is
  a local stand-in to be replaced by `biokotlin-core` (ADR-0008).
- Reusing `file(...)` for both: no DSL change, but it lies to every front end
  about what the value is, and it was rejected in #45 for that reason.
- A repeatable CLI option (`--taxa-list A --taxa-list B`) instead of one
  comma-separated value: idiomatic for a CLI, but the engine takes command-line
  text for every parameter, so the CLI would own a spelling the spec does not,
  which is the drift ADR-0002 exists to prevent.
