# ADR-0003: Compose Multiplatform for desktop and web

- Status: Accepted
- Date: 2026-10-29

## Context

TASSEL 5's Swing GUI is dated and JVM-only. TASSEL 6 needs a modern desktop app
and a browser app, built and maintained mostly by one developer with AI agents.

## Decision

`tassel-app` uses Compose Multiplatform with a shared `commonMain` UI, a
`desktopMain` (JVM) entry point, and a `wasmJsMain` (browser) entry point. The
UI talks to analyses only through an `Engine` interface (see ADR-0004).
Visualizations (genotype grid, kinship heatmap, PCA scatter, later the
haplotype viewer) are drawn with Compose `Canvas`.

## Consequences

- One UI codebase in Kotlin, sharing types with the core.
- Analysis forms are generated from operation specs (ADR-0002).
- Compose for Web (Wasm) is newer than desktop Compose; expect rough edges in
  text input, accessibility, and initial bundle size. Track them in issues.
- Large data views must be virtualized and canvas-drawn to perform in the browser.

## Alternatives considered

- Kotlin core plus a React/TypeScript web UI: richer web ecosystem, but two UI
  codebases and a JS interop boundary.
- JavaFX or Swing refresh: desktop only.
- Tauri/Electron shell over a web UI: adds a second language and runtime.
