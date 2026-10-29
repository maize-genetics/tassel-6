<!--
PR title must be a Conventional Commit: feat(core): port Centered IBS kinship
Target branch: develop (release/* and hotfix/* follow git flow; see ADR-0010).
-->

Closes #

## Summary

<!-- What this PR does and why, in two or three sentences. -->

## Checklist

- [ ] Scope matches the linked issue; nothing unrelated is included
- [ ] No `java.*` or other JVM-only dependencies in any `commonMain`
- [ ] No new dependencies, or they are listed here with a reason:
- [ ] No changes under `testdata/golden/` (or this PR has the `golden-update` label and explains why)
- [ ] No changes to frozen contracts (or the ADR is linked)
- [ ] `docs/porting-ledger.md` updated, if anything was ported, bridged, or retired

## Verification

<!-- Real command output, not a description of it. -->

```text
$ ./gradlew :tassel-core:allTests

```

## Handoff

**Status:** done | partial | blocked

**What changed**
-

**Decisions made**
-

**Deviations from TASSEL 5**
- None

**Open questions for the maintainer**
-

**Next step**
-

**Out of scope, noticed**
-
