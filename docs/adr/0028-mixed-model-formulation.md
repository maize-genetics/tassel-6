# ADR-0028: Mixed linear model formulation, REML, P3D, and failure semantics

- Status: Accepted
- Date: 2026-11-25

## Context

TASSEL 5's MLM (`analysis/association/MLMPlugin.java`; the pipeline's `-mlm` runs
`WeightedMLMPlugin`, which extends it) fits `y = X b + Z u + e` with `Var(u) = sigma_g^2 K` and
`Var(e) = sigma_e^2 I`, estimating variance components with its EMMA port
(`stats/EMMA/EMMAforDoubleMatrix.java`) and testing markers in
`CompressedMLMusingDoubleMatrix.java`. The inventory in `docs/llm_notes/ASSOCIATION_PLAN.md`
section 2 records how, and what goes wrong:

- **EMMA** eigendecomposes `Z K Z'`, adds a "bend" when its smallest eigenvalue is below 0.01,
  eigendecomposes the projected `S (Z K Z' + bend I) S`, and maximizes a REML log-likelihood in
  `delta = sigma_e^2 / sigma_g^2` over 100 log-spaced points on `[1e-5, 1e5]`, then Newton–Raphson
  inside each bracket (`|d1| < 1e-10`, 50 iterations, nested rescans). Boundary solutions are
  accepted silently. The reported `-2LnLikelihood` includes a constant computed with integer
  division and omits `log|X'X|`.
- **P3D** (the default) fixes `sigma_g^2` and `sigma_e^2` from the null model, inverts `V` by a full
  sweep, removes each marker's missing rows by reverse sweeps (exact in theory, error-accumulating
  in practice), and tests the marker with a Wald F, `(M b)' (M (X'V^-1 X)^-1 M')^-1 (M b) / df_m`,
  against `F(df_m, n - ncol(X))`. Explicit inverses ignore their failure flag.
- **Per-marker estimation** (`EachMarker`) reruns EMMA with the marker in `X`, at `O(n^3)` per
  marker, and silently drops markers whose `X` is rank-deficient.
- **Compression** clusters taxa by UPGMA on `c - K` and replaces `K` by group means; the pipeline
  default is no compression, while the GUI and documentation describe `Optimum`.
- **Reports** have known defects: `errordf` does not equal the df used for `p`, there is no
  heritability column, and a failed fit prints nothing useful.

The brief asks for uncompressed REML MLM with P3D first, per-marker re-estimation as a distinct
option, explicit handling of singular or invalid kinship, convergence diagnostics, documented
heritability assumptions, correct replicates, and never a silent fall back to GLM.

## Decision

### Model

- `y = X b + Z u + e`, `u ~ N(0, sigma_g^2 K)`, `e ~ N(0, sigma_e^2 I)`, so
  `V = sigma_g^2 (Z K Z' + delta I)`. `Z` maps observations to the kinship taxa (ADR-0026);
  unreplicated data use `Z = I` without forming it.
- One random effect in this epic. The term hierarchy (ADR-0025) leaves room for more.

### Estimation

- **Spectral form.** Eigendecompose `Z K Z' = U diag(lambda) U'` once per trait observation set,
  with progress reports. Rotate `y` and `X` by `U'`. The REML log-likelihood at any `delta` is
  evaluated in the rotated basis in `O(n p^2)` (the FaST-LMM/GEMMA form). It is the same function
  EMMA maximizes, up to the constants TASSEL 5 drops, so the optimum is the same.
- **Search.** TASSEL 5's grid (100 log10-spaced points on `[1e-5, 1e5]`) and Newton refinement
  inside each bracket, with its tolerance and iteration cap, are the default for parity. The best
  of the grid maximum and the refined optima wins, as in TASSEL 5. A later issue may add a wider
  grid or Brent's method as an option; it may not change the default without the ADR-0029
  departure procedure.
- **Missing rows by indicator columns.** An observation dropped for one test (a missing genotype)
  is modeled by adding an indicator column for it to `X` rather than by refactoring `V`. For known
  `V` this gives exactly the GLS estimates, residual quadratic form, and df of the subset analysis,
  and as a function of `delta` the REML log-likelihood equals the subset's (the `log|V|` and
  `log|X'V^-1 X|` terms combine to the subset's, and the residual df match). So one eigendecomposition serves every marker, P3D or re-estimated, with no
  imputation. When a marker has many missing rows, a per-pattern decomposition may be cheaper; the
  choice is a performance heuristic and must not change results beyond ADR-0029 tolerances.
- **Kinship checks** (ADR-0026) run first. TASSEL 5's "bend" for a slightly indefinite `K` is
  reproduced only as an explicit, recorded option needed for parity on matrices TASSEL 5 bends; a
  materially non-PSD kinship is an error.

### Tests

- **P3D (default).** `delta`, `sigma_g^2`, and `sigma_e^2` come from the null model. The marker
  statistic is TASSEL 5's Wald F with the variance components fixed, against `F(df_m, n - p)`
  where `n - p` is the residual df after the marker and any indicator columns. In the rotated
  basis this is the GLM FWL statistic with `sigma_g^2` in place of the residual mean square.
- **Re-estimated (`EachMarker`).** REML is rerun with the marker in `X` (reusing `U`), and the test
  uses that marker's variance components. Rank-deficient markers produce a "not estimable" record,
  not a silent drop.
- Effects, additive and dominance contrasts, and TASSEL 5's GLS `MarkerR2` (Buse) are reported as
  TASSEL 5 reports them in the compatible report.

### Reported quantities and their meaning

- Every record carries `sigma_g^2`, `sigma_e^2`, `h^2 = sigma_g^2 / (sigma_g^2 + sigma_e^2)`, the
  REML log-likelihood with all constants (and TASSEL 5's `-2LnLikelihood` value in the compatible
  report), the variance method (`p3d` or `re-estimated`), and status flags: `boundary-low`,
  `boundary-high` (`delta` at a grid end), `not-converged`, `not-estimable`.
- **Kinship scaling.** `h^2` is relative to `K` as given. Provenance records the kinship source
  and its mean diagonal and `trace(K) / n`, because two kinships of different scale give different
  `sigma_g^2` and `h^2` and the same marker tests. No automatic rescaling.
- **Failure.** A null model that fails (non-finite likelihood everywhere, invalid kinship, no
  residual df) yields failed records for that trait with the reason. MLM never substitutes GLM.

### Compression

Compressed MLM (UPGMA clustering, group kinship, the `Optimum` and custom levels) is a later issue
in this epic (#199) that reproduces TASSEL 5 for parity. The TASSEL 6 default is no compression,
which is TASSEL 5's pipeline default (its GUI and documentation suggest `Optimum`). Compressed MLM
follows the first MLM release rather than blocking it (plan decision Q5, 2026-11-25).

## Consequences

- P3D MLM costs one `O(n^3)` decomposition per trait plus an `O(n^2)` rotation per marker (blocked
  as a matrix product), and re-estimation costs `O(n p^2)` per likelihood evaluation instead of
  `O(n^3)` per marker.
- Missing genotypes are handled exactly and without imputation at the price of a few extra columns
  per marker.
- Reported `-2LnLikelihood` in the compatible report reproduces TASSEL 5's constant; the modern
  report gives the full REML log-likelihood. The departure table explains both.

## Alternatives considered

- **Port EMMA's projected decomposition per marker.** Parity by construction, and `O(n^3)` per
  marker for re-estimation. Kept only as a test-only reference if parity requires it.
- **Mean-impute missing genotypes** (GEMMA, GCTA). Fast and simple, and it changes results without
  the user asking (ADR-0026). Rejected as a default.
- **Refactor `V` per missingness pattern.** Exact, and `O(n^3)` per pattern. Kept as a heuristic
  for markers with many missing rows.
- **AI-REML or EM** instead of a grid. Faster to converge with several random effects; with one
  random effect the 1-D search is cheap and matches TASSEL 5. Deferred.
