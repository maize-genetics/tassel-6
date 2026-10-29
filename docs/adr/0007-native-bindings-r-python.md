# ADR-0007: R and Python bindings through a Kotlin/Native C API

- Status: Accepted
- Date: 2026-10-29

## Context

rTASSEL 1.0 calls TASSEL 5 through rJava, which makes installation fragile
(JDK discovery, memory flags) and slow at the boundary. Python users have no
supported binding.

## Decision

- `bindings/capi` builds `tassel-core` as a Kotlin/Native shared library with a
  small, stable C API (opaque handles, operation-spec-driven calls, error codes).
- `bindings/r` is rTASSEL 2.0, which keeps the rTASSEL 1.0 user-facing API where
  feasible, and `bindings/python` is the Python package. Both live in this
  monorepo and call the C API.
- Analyses not yet ported natively are delegated to a local `tassel serve` JVM
  process (which can use the bridge), behind the same R/Python functions.

## Consequences

- No JVM is required for natively ported analyses in R or Python.
- The C API must be versioned and kept small; spec-driven calls avoid one C
  function per analysis.
- Native builds per OS and architecture are needed for CRAN/PyPI packaging.
- The delegation path requires a JVM only when a legacy analysis is used.

## Alternatives considered

- Keep rJava for R and use JPype for Python: keeps the installation problems.
- Bindings only through a server: needs a running process for everything.
