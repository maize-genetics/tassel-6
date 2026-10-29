# Handoff template

Paste this into the PR body at the end of every session (the PR template
includes it). For work that spans several PRs or sessions, also save a copy as
`docs/handoff/<issue>-<slug>.md` and update it each session.

Write for a fresh agent with no memory of this session: no private shorthand,
and name files, classes, and commands explicitly.

```markdown
## Handoff

**Status:** done | partial | blocked

**What changed**
- <file or module>: <what and why, one line each>

**Verification**
- Command: `./gradlew ...`
- Result: <pass/fail counts per target; paste the relevant output below>

**Decisions made**
- <decision, and the ADR it follows or the ADR it needs>

**Deviations from TASSEL 5**
- <none, or: behavior, reason, and whether golden tests cover it>

**Open questions for the maintainer**
- <question needing a human decision>

**Next step**
- <the single most useful next task, specific enough to paste into an issue>

**Out of scope, noticed**
- <problems seen but intentionally not fixed here>
```
