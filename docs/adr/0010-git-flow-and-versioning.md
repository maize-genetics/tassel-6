# ADR-0010: Git flow branching, Conventional Commits, SemVer

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 5 syncs `main` to `develop` and publishes nightlies. TASSEL 6 will be
built by one maintainer with several AI agents working in parallel, so branch
discipline and readable history matter more than usual.

## Decision

- Full git flow: `main` holds released versions, `develop` is integration,
  `feature/<issue>-<slug>` branches target `develop`, `release/<version>`
  branches stabilize releases, and `hotfix/<slug>` branches fix `main`.
- `main` and `develop` are protected: PRs only, CI required, no force pushes.
- Commits and PR titles follow Conventional Commits; a CI check enforces PR titles.
- Versions follow SemVer: `6.0.0-poc.1`, `6.0.0-alpha.N`, `6.0.0-beta.N`,
  `6.0.0`; nightly builds are `6.0.0-dev.YYYYMMDD`.
- Changelogs are generated with git-cliff from commit history.

## Consequences

- Release notes come for free from commit messages.
- Releases merge `release/*` into `main`, tag, and merge back into `develop`.
- Squash-merge feature PRs so each PR is one conventional commit on `develop`.

## Alternatives considered

- Trunk-based development: simpler, but the user wants release branches and a
  stable `main` matching TASSEL 5 practice.
