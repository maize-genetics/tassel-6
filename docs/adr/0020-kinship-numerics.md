# ADR-0020: Kinship numerics, tolerance parity, and parallelism

- Status: Accepted
- Date: 2026-11-06

## Context

Centered IBS in `distance.centeredIbs` was a bit-exact transcription of TASSEL 5
`EndelmanDistanceMatrix`: `Float` sums, site blocks of
`max(numberOfSites / 4, 1000)`, and a right-to-left reduction. That is what
makes `CenteredIbsTest.writesTheGoldenByteForByte` pass, and it is also why the
kernel cannot use a blocked matrix product or a host-sized thread count
without moving the last bits.

ADR-0009 sets the kinship golden tolerance at `abs 8e-6` and says an
intentional departure from TASSEL 5 needs an ADR. This is that ADR for
kinship. It refines ADR-0009; it does not supersede it. The maintainer's
decision, recorded in tassel-6#62, is that a modernized kinship algorithm
must pass the existing goldens within `abs 8e-6`, including the
`Centered_IBS.SumPk` annotation, and may become the default. Byte-for-byte
output is no longer required for the default kinship path.

No file under `testdata/golden/` changes. The golden text stays the TASSEL 5
`v5.2.98` output.

## Decision

- **The default kernel may change summation order.** `kinship` calls the
  blocked dosage product. Agreement with `testdata/golden/kinship/centered_ibs.txt`
  is absolute tolerance `8e-6` on every matrix value, the tolerance ADR-0009
  and the MANIFEST already name.
- **The TASSEL 5 summation order stays available to tests.** The bit-exact
  kernel remains `internal` in `commonMain` (`centeredIbsTassel5Order`) so
  every target can still run it. It is not a `method` or `algorithm` parameter
  on the operation spec. Published bit-exact reproduction is TASSEL 5;
  TASSEL 6's default guarantees `abs 8e-6`.
- **`Centered_IBS.SumPk` is compared as a number.** Parse the annotation and
  require absolute difference at most `1e-4` from the golden
  `918.7297856754799`. That bound covers a `Double` frequency against TASSEL 5's
  `Float` frequency over a few thousand sites, and a dropped site moves the sum
  by about `0.1`. The writer stays `doubleToString`. Exact annotation text is
  required only of `centeredIbsTassel5Order`.
- **`DistanceMatrix` stays `Float`-backed.** Pair sums accumulate in `Double`
  in the default kernel, are divided by `SumPk`, and are rounded to `Float`
  when the matrix is built. The text format is the precision ceiling.
- **Parallelism is not a user-facing parameter.** Workers come from the
  platform dispatcher (`Dispatchers.Default`), one worker on wasmJs. Results
  must not depend on the host core count: the default kernel's panel size is a
  function of the taxa and site counts only, and the retained kernel keeps
  `GOLDEN_MAX_THREADS = 4` so its block size stays the golden one. This
  diverges from TASSEL 5, where `-maxThreads` also changes the block size and
  therefore the last bits. A `threads` spec parameter (legacy name
  `maxThreads`) waits until a BLAS path needs a user cap.
- **Exact-text tests stay on the retained kernel.**
  `CenteredIbsTest.writesTheGoldenByteForByte` calls `centeredIbsTassel5Order`.
  The CLI golden check (`KinshipEndToEndTest`) follows the default kernel and
  uses `abs 8e-6` plus the `SumPk` bound above. `matchesHandComputedValuesOnTwoTaxa`
  stays an exact string: one site, `p = 0.5`, both kernels.

## Consequences

- tassel-6#67 can replace the default kernel. tassel-6#65 and tassel-6#66 keep
  the retained kernel byte-identical, including its `Float` reduction.
- `CenteredIbsTest.matchesGoldenSumPk` and `KinshipEndToEndTest` no longer
  require the default writer to reproduce TASSEL 5's `SumPk` characters.
  Loosening a tolerance to make a failing port pass is still forbidden;
  `1e-4` and `8e-6` are this ADR's contract, not a test-only escape hatch.
- Two kernels remain until a later issue deletes the retained one. The cost is
  maintenance. The benefit is a differential test: missing data, monomorphic
  sites, extra alleles, and `maxAlleles` of 2 and 6 are compared against the
  bit-exact kernel rather than against a second golden.
- Host core count does not change kinship bytes. TASSEL 5's `-maxThreads`
  sweep is a timing tool only.

## Alternatives considered

- **Drop the TASSEL 5 order entirely.** Smaller surface, and no way to tell a
  wrong allele rule from ordinary `Double` noise. Rejected.
- **Expose the old order as a spec parameter.** Keeps a reproducibility flag
  on the CLI, and commits the project to two supported kernels. Rejected;
  tests are enough.
- **Require exact `SumPk` text from the default kernel.** Forces the new
  accumulator to replay TASSEL 5's `Float` frequency and block reduction, which
  is the constraint this ADR exists to lift. Rejected.
- **Store the matrix as `Double`.** The writer and every TASSEL 5 consumer
  would change. Rejected.
- **A `threads` parameter now.** Nothing in the operation needs it for results
  to be reproducible, and TASSEL 5's legacy name also changes numerics.
  Rejected until a fast path has a pool to cap.
