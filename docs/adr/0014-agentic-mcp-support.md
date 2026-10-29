# ADR-0014: Agent support: local MCP server, in-app assistant, CLI skill

- Status: Accepted
- Date: 2026-10-29

## Context

Users increasingly run analyses through AI agents. The lab already has a
`tassel_cli_skill` for TASSEL 5 and a Koog-based agent pattern in `sra_koog`.

## Decision

- `tassel mcp` runs a local stdio MCP server, built with the Kotlin MCP SDK,
  that exposes operation specs as tools (ADR-0002) and data handles as resources.
- The app includes an optional in-app assistant built with Koog that uses the
  same tools through the `Engine`.
- The `tassel-cli` agent skill is regenerated from the specs on each release.

## Consequences

- MCP tools, CLI help, and the skill stay in sync automatically.
- Large data never passes through the model: tools take and return handles,
  files, and summaries.
- The in-app assistant needs a user-supplied model API key or local model and
  must work without network access when disabled.

## Alternatives considered

- Remote hosted MCP server: needs hosting and data upload; may come with the server tier.
- Skill only: no structured tool calls.
