# ADR-0016: Agent-driven development workflow

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 6 is built by one maintainer with AI agents. Agent sessions have limited
context and no memory between chats, so project knowledge must live in the
repository, and work must be split into tasks that one session can finish and
verify.

## Decision

- Repository memory: `AGENTS.md`, scoped `.cursor/rules/`, ADRs, the porting
  ledger, and handoffs written into PR bodies.
- One task = one GitHub issue = one fresh chat = one `feature/*` branch = one PR
  into `develop`. Issues use the agent-task template.
- Work proceeds in waves: foundations (skeleton, CI, goldens), then frozen core
  contracts, then parallel ports and front ends against those contracts, then
  integration.
- At most two to four agent PRs are open at once. Headless ports may run as
  cloud agents or in worktrees; UI and integration work runs locally.
- Golden test data and contracts are protected; agents stop and ask rather than
  changing them.
- Automated review (Bugbot) and CI run on every PR. Only the maintainer merges.

## Consequences

- Every session starts from the same kickoff: read `AGENTS.md`, the issue, the
  cited roadmap section and ADRs, and the latest handoff.
- The cost is writing good issues; the payoff is sessions that do not need
  re-explaining.

## Alternatives considered

- Long-running single sessions: lose context and drift.
- Many parallel agents without frozen contracts: merge conflicts and
  inconsistent APIs.
