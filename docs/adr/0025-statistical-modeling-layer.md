# ADR-0025: A decomposition-based statistical modeling layer in `commonMain`

- Status: Proposed
- Date: 2026-11-25

## Context

Association mapping (GLM, MLM, and later stepwise, logistic, and genomic prediction) needs dense
linear algebra and linear-model machinery that TASSEL 6 does not have. Today `tassel-core` has one
decomposition, `math.SymmetricEigenDecomposition` (a JAMA `tred2`/`tql2` port used by PCA), plus
the Cephes special functions and the F, t, and chi-square distributions from epic #155. There is no
matrix type, no QR, no Cholesky, no linear solver, no design-matrix builder, and no model
specification.

TASSEL 5 builds its linear models on `matrixalgebra.DoubleMatrix`, backed by Colt, EJML, or BLAS
through a factory, and on `stats.linearmodels.SweepFastLinearModel`, which fits by sweeping `X'X`
and computes a generalized inverse. Several TASSEL 5 association paths invert matrices explicitly.
None of that can be linked from `commonMain` (ADR-0001), and copying the sweep design would also
copy its numerical weaknesses on collinear covariates.

The planning brief for association mapping (epic to be opened) asks for decomposition-based
numerics, no explicit inversion, rank-deficient model handling, reusable fixed-effect
decompositions (Frisch–Waugh–Lovell), and stable log p-values, without overengineering the first
release.

## Decision

### Packages

- `net.maizegenetics.tassel.math` (existing) gains the dense kernels. Everything is `Double`,
  column-major `DoubleArray` storage, `kotlin.math` only, no new dependencies:
  - `DenseMatrix` (rows, columns, column-major data) and `DenseVector` helpers. Views and
    workspaces are explicit; nothing allocates per marker in hot loops.
  - `HouseholderQr` with column pivoting and a rank-revealing tolerance
    `tol = max(n, p) * eps * |r11|` (the LAPACK `dgeqp3` convention), exposing rank, the pivot,
    `Q'b` and `Qb` applications, and triangular solves. `Q` is never formed unless asked for.
  - `Cholesky` (lower triangle; reports the failing pivot as a value instead of throwing deep in a
    loop), and `LDL'` only if a later issue needs it.
  - `SymmetricEigenDecomposition` stays the eigensolver. A values-and-vectors-into-workspace overload
    is added when MLM needs repeated decompositions.
  - Triangular solves and `BLAS`-style level-1/2/3 kernels (`dot`, `axpy`, `gemv`, `syrk`, `gemm`)
    written so a later platform backend can replace them behind an interface without changing
    callers.
- `net.maizegenetics.tassel.stats.model` (new) holds the modeling layer:
  - `ModelSpec`: response, fixed terms, and (for mixed models) random terms. Terms are a sealed
    hierarchy: `Intercept`, `Covariate(name)`, `Factor(name, coding)`, `MarkerTerm(encoding)`, and
    `RandomEffect(kinship)`. A formula string (`Height ~ Marker + PC1 + PC2 + Env`) is an optional
    front end that parses into a `ModelSpec`; it is a separate issue.
  - `Contrast` coding for factors: `treatment` (TASSEL 5's drop-last-level convention, so effect
    reports line up) and `sum`. Level order is the `FactorAttribute` order (sorted, ADR-0023).
    TASSEL 5 numbers nuisance-factor levels by first appearance instead; that changes only which
    level is the reference, not marker tests, marker effects, or BLUEs. Marker classes keep
    TASSEL 5's order (ADR-0026).
  - `DesignMatrix`: built from a `PhenotypeTable` and an observation mask, recording which columns
    belong to which term, so sums of squares and Wald tests can be computed per term.
  - `LinearModelFit` (OLS): coefficients, rank, residuals, fitted values, residual sum of squares,
    and residual degrees of freedom, all through the pivoted QR. Effective model df is the
    detected rank, never the column count.
  - Tests: `FTest` (nested models, any numerator df), `WaldTest` (contrast `L`), and `TTest`
    (single coefficient), each returning statistic, df, p, and log10 p.
  - `GeneralizedLeastSquares`: the same fit after a whitening transform supplied by a covariance
    model; MLM (ADR-0028) uses it.
- Effect standard errors come from triangular solves against `R`, never from `(X'X)^-1`.

### Numerics rules

- **No explicit inverse** in production paths. Where TASSEL 5 inverts, the port solves, and the
  difference is covered by the parity tolerance in ADR-0029.
- **Rank deficiency is reported, not thrown.** A collinear covariate reduces the model rank and the
  model df; the fit records which columns were aliased. A marker term that adds no rank yields a
  "not estimable" result for that marker rather than an exception.
- **Frisch–Waugh–Lovell reuse.** A marker test on top of a fixed null design `X0` may project the
  marker and response onto the orthogonal complement of `X0` with a cached QR of `X0`, as long as
  the marker's observation set equals the cached one. A cached projection is never applied to a
  different observation set. When a marker drops rows (missing genotypes), the kernel either
  appends one indicator column per dropped row to `X0` (an exact identity: estimates, residual sum
  of squares, rank, and df equal those of the subset fit, for OLS and for GLS with known `V`) or
  uses a QR cached for that missingness pattern, choosing by cost. Both give the subset fit's
  answer to rounding.
- **p-values are computed in the upper tail and in log space.** Every test reports `p` and
  `log10(p)`. `log10(p)` comes from a log-space evaluation of the regularized incomplete beta or
  gamma function, so it stays finite when `p` underflows to `0.0` (`FDistribution.logSf`,
  `StudentTDistribution.logSf`, `ChiSquaredDistribution.logSf`; a prerequisite issue).
- **Determinism.** Results do not depend on thread count or platform beyond the documented
  tolerances: reductions inside a model fit are sequential, and parallelism is across markers or
  traits, with results written in a fixed order (ADR-0020's rule, extended).
- **A TASSEL 5-order reference path may be kept** for a kernel whose default is reorganized for
  speed, as ADR-0020 did for kinship, so differential tests can separate "wrong rule" from
  "different summation order".

### Scope guard

The first release implements what GLM and P3D MLM need: OLS and GLS fits, nested F tests, Wald and
t tests, treatment and sum contrasts, and the kernels above. Interactions, nested and random
factors beyond one kinship term, iterative solvers, sparse matrices, and logistic regression are
designed for (the term hierarchy is sealed but open to new subclasses) and left to later issues.

## Consequences

- One linear-model implementation serves GLM, MLM, stepwise, and future methods, on every target.
- Decomposition-based fits are slower per call than a sweep on tiny models, and that is acceptable:
  the GLM hot loop is FWL residualization, which is `O(n * p)` per marker after one `O(n * p^2)` QR.
- TASSEL 5 outputs that depend on its sweep's handling of near-collinear designs can differ in the
  last digits; ADR-0029 sets the tolerances and the procedure for documenting such departures.
- A BLAS, Panama FFM, or GPU backend (Phase 3, ADR-0004) can replace the level-3 kernels later
  without touching model code.

## Alternatives considered

- **Port `SweepFastLinearModel` and the `DoubleMatrix` factory.** Fastest path to parity, and it
  carries explicit inversion, the sweep's collinearity handling, and a JVM-shaped factory into the
  shared core. Rejected; the sweep's results remain the parity oracle.
- **Normal equations with Cholesky.** Simple and fast, and squares the condition number, which
  matters for PCs plus factors plus markers. Rejected as the default; kept for small, well-scaled
  systems where a later issue measures a benefit.
- **A third-party multiplatform linear algebra library.** None is pure Kotlin, maintained, and
  available on jvm, wasmJs, and native with the decompositions needed. Rejected under ADR-0001 and
  AGENTS.md constraint 2.
- **Formula strings as the only model interface.** Pleasant for users, and it makes the operation
  spec and MCP schema depend on a parser. Rejected; typed specs first, formulas as a front end.
