# ADR-0026: Association inputs, taxon alignment, missing data, and marker encoding

- Status: Proposed
- Date: 2026-11-25

## Context

An association analysis combines four inputs that arrive separately and must agree on taxa: a
genotype table, one or more phenotype tables (traits, covariates, factors, possibly replicated),
optional principal components, and, for mixed models, a kinship matrix. Epic #155 delivered the
pieces: `PhenotypeTable` (ADR-0023), `DataSet` and joins (ADR-0024), and `GenotypePhenotype`, which
maps observations to genotype rows without copying genotypes. PCA already returns its scores as a
covariate `PhenotypeTable`, so it can be joined like any other covariate file.

What is missing is a rule set for how a model sees these inputs:

- which observations a test uses when the trait, a covariate, or the marker is missing;
- how replicated observations relate to unique taxa, especially for the kinship random effect;
- how kinship is read, validated, and subset to the analysis taxa;
- how a genotype becomes model columns. `GenotypeTable` (frozen, ADR-0012) stores diploid
  nucleotide bytes. It has allele identity (`AlleleStats`, `Position.knownVariants`) but no
  dosages, probabilities, or haplotype alleles, and future marker sources (dosages, genotype
  probabilities, SVs, haplotypes, CNVs, k-mers, PHG-derived calls) will not all be nucleotide
  bytes.

TASSEL 5 (inventory in `docs/llm_notes/ASSOCIATION_PLAN.md` section 2):

- treats each genotype class as a factor level (the omnibus "marker" test) in both GLM and MLM,
  with the last class in its sort order (usually the heterozygote) as the zero reference, and adds
  additive and dominance F tests for biallelic sites. GLM sorts classes as IUPAC strings; MLM sorts
  homozygotes by byte, then heterozygotes;
- can test `ReferenceProbability` or `AlleleProbability` as numeric covariates in GLM (the latter
  path is broken);
- drops each marker's observations with missing calls, in GLM by refitting and in MLM by
  reverse-sweeping rows out of `V^-1`; it never imputes in GLM or MLM;
- sets genotype classes smaller than `minClassSize` to missing (GLM only), and has no MAF or
  minimum-count filter in GLM or MLM;
- **refuses replicated taxa in GLM** (`AbstractFixedEffectLM.testTaxaReplication` throws). The
  taxa-within-marker error term its manual describes is unreachable code. MLM handles replicates
  through `Z`, and `PhenotypeLM` through a taxon factor;
- counts phased heterozygotes (`A/C`, `C/A`) as different classes in its GLM filters but as one
  class in the model.

## Decision

### Inputs and alignment

- `glm` and `mlm` take a `DataSet` (ADR-0024). Accepted contents: one `GenotypePhenotype`, or one
  `GenotypeTable` plus one or more `PhenotypeTable`s (joined with `join-genotype-phenotype`'s
  intersect rules first), and for `mlm` exactly one `DistanceMatrix` as kinship. A `glm` run for
  BLUEs (`phenoOnly`) takes phenotypes alone. Anything else is an error naming what was found. The
  CLI already accepts `DataSet` inputs; the app needs kinship to become a selectable dataset (D2 in
  the plan).
- `DataSet` gains `distanceMatrices()` and `genotypePhenotypes()` accessors. This is additive;
  `Operation`, `RunRequest`, and `Engine` do not change.
- Alignment is by exact, case-sensitive taxon name throughout (ADR-0024). Kinship is subset to the
  analysis taxa by name, in observation order of first appearance, without copying the full matrix
  more than once. A phenotype taxon absent from the kinship matrix is excluded and counted, never
  silently given a zero row.
- Duplicated taxon names in a genotype table or kinship matrix are an error. Replicated
  observations in a phenotype are normal (ADR-0023).

### Observations, missing data, and counts

- **Per trait**, the base observation set is the rows where the trait, every covariate, and every
  factor in the model are present. Excluded rows are counted by reason.
- **Per marker**, observations with a missing genotype (or with a genotype class set to missing by
  `minClassSize`) are dropped for that marker's test only. TASSEL 5 does this, and so does every
  GLM; it is why ADR-0025's FWL path handles dropped rows with indicator columns or a
  pattern-keyed cache.
- **No silent imputation.** GLM and MLM drop a marker's missing rows exactly, as TASSEL 5 does.
  Mean or other imputation happens only when the analysis configuration names it
  (`missingGenotypes = impute-mean`), and the provenance records it. ADR-0025 and ADR-0028 describe
  how dropped rows are handled without refactoring per marker.
- Every result row reports both the **observation count** and the **unique taxon count**.

### Replicated observations

- **GLM with replicated taxa is new TASSEL 6 behavior** (TASSEL 5 throws). The model adds taxa
  nested within marker classes, and the marker is tested against the taxa-within-marker mean
  square with that term's df, the design TASSEL 5's manual describes. The report names the error
  term. It is validated against R (`aov` with an error stratum), not a TASSEL 5 golden. Averaging
  replicates first (BLUEs, then GLM) and MLM remain the alternatives users are pointed to. Decision
  Q1 in the plan confirms this before B4 implements it.
- MLM models replicates through an incidence matrix `Z` mapping observations to taxa, so `K` stays
  taxa by taxa (`V = sigma_g^2 Z K Z' + sigma_e^2 I`). Unreplicated data use `Z = I` without
  forming it.

### Kinship

- A core reader for TASSEL 5's square-matrix text (`##` annotations, taxa count, labeled rows), the
  format `DistanceMatrix.writeDelimited` already writes, returning a `DistanceMatrix`. `Float`
  storage stays (ADR-0020), so a kinship file read from TASSEL 5 gives TASSEL 6 the same values
  TASSEL 5 used.
- Validation before fitting: square, symmetric within a tolerance, finite, unique labels, and
  positive-semidefinite up to a stated eigenvalue tolerance. A failure names the problem and
  stops; it is never repaired silently (bending or adding a ridge is an explicit option, if added).

### Marker encodings

- An internal `MarkerSource` turns site blocks into model columns for the analysis observations. It
  is the seam for future marker kinds; `GenotypeTable` is its first and only implementation in this
  epic. Encodings:
  - `genotype`: each observed diploid genotype class is a factor level (TASSEL 5's omnibus test,
    `k - 1` df). This is the TASSEL 5-compatible default for `glm` and `mlm` (decision Q3). Phased
    heterozygotes are one class everywhere, filters included (a departure only for phased input).
    The class order and reference class follow TASSEL 5 per analysis (IUPAC string order for GLM,
    homozygotes-then-heterozygotes for MLM), so effect reports line up.
  - `additive`: count of a named effect allele (0, 1, 2), one df. The effect allele is the minor
    allele among analysis observations unless the configuration names the major allele or a
    specific allele (decision Q2; TASSEL 6 has no reference allele). The record states the effect
    and comparison alleles explicitly.
  - `additive-dominance`: additive plus a heterozygote indicator, reported as the TASSEL 5 `add_`
    and `dom_` tests for biallelic sites. TASSEL 5's `add_F` codes every heterozygote as 1 for
    every allele, which is wrong at multiallelic sites; TASSEL 6 counts allele copies (a recorded
    departure where it changes a golden).
  - `dosage`: a numeric value per observation from a numeric marker source (reference
    probability today; dosages and probabilities later).
- **Multiallelic sites.** `genotype` handles any number of classes. `additive` applies to biallelic
  sites; for multiallelic sites it is either skipped with a typed reason or expanded to one column
  per non-reference allele, chosen by configuration, never silently collapsed.
- **The omnibus test and the 1-df additive test are different tests** and are labeled as such in
  every report. Neither is computed and reported under the other's name.
- Allele labels in effect reports use TASSEL 5's IUPAC genotype strings in the compatible report
  and explicit allele strings (from `Position.knownVariants` when present) in the modern report.

### Filters

`minClassSize` and `biallelicOnly` keep TASSEL 5's GLM semantics, evaluated on the trait's
analysis observations. Minimum MAF and minimum minor allele count are new (TASSEL 5 has none;
users filter sites first), default to off, and are evaluated on the same observations, so a marker
can be tested for one trait and skipped for another. A skipped marker is reported with its reason,
not dropped silently.

## Consequences

- `GenotypeTable` stays frozen; new marker kinds plug in behind `MarkerSource` without a contract
  change.
- Kinship becomes a first-class input on every front end.
- Counts and exclusion reasons make the analysis auditable, at the cost of a few more columns.

## Alternatives considered

- **Add dosage and allele accessors to `GenotypeTable`.** Changes a frozen contract for something an
  adapter can do. Rejected.
- **Impute missing genotypes by default, as many GWAS tools do.** Faster for MLM, and it changes
  results without the user asking. Rejected as a silent default.
- **Align by taxon position rather than name.** Fragile with unsorted or partly overlapping inputs.
  Rejected.
