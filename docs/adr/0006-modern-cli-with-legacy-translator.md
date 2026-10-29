# ADR-0006: Modern CLI with a TASSEL 5 pipeline translator

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 5's `run_pipeline.pl` uses a positional, stack-based flag syntax
(`-h file -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin -export`)
that is hard to read and validate. Many published workflows and scripts depend
on it.

## Decision

- `tassel-cli` is a new Clikt-based CLI with kebab-case subcommands generated
  from operation specs, following the phg_v2 CLI conventions.
- `tassel translate` converts `run_pipeline.pl` command lines into equivalent
  TASSEL 6 commands, using the `legacy` names on spec parameters. Plugins that
  have no native port are routed to the bridge (ADR-0005) or reported.

## Consequences

- New users get a conventional CLI with `--help` from the specs.
- Existing scripts have a migration path without TASSEL 6 carrying the old
  syntax as its primary interface.
- The translator's coverage is a measurable compatibility metric.

## Alternatives considered

- Re-implement `run_pipeline.pl` syntax as the primary CLI: preserves a poor UX.
- New CLI only: abandons existing workflows.
