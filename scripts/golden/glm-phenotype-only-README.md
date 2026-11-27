# glm-phenotype-only

TASSEL 5 `v5.2.98` `FixedEffectLMPlugin -phenoOnly true` (`PhenotypeLM`): for every trait, the model
mean + factors + covariates + taxa, the best linear unbiased estimate (BLUE) of each taxon, and the
test of the taxon effect. Each run writes the BLUEs as a phenotype (`NAME.blues.txt`, floats, compared
as a phenotype with rows in any order, because TASSEL 5 sorts the taxa) and the ANOVA report
(`NAME.anova.txt` rounded as `-export` writes it, `NAME.anova.full.txt` with `Double.toString`).

| Run | Input | Covers |
|---|---|---|
| `earht_with_rep` | `phenotype/earht_with_rep.txt` (30 taxa, 58 rows, factor `rep`) | Replicated taxa; TASSEL 5's `PhenotypeLMTest` checks the same BLUEs against R |
| `mdp_traits` | `tutorial/mdp_traits.txt` (three traits with different missingness) | Unreplicated taxa: error df 0, TASSEL 5 prints rounding noise (`F = -0`, `errorMS = Infinity`) |
| `synthetic` | `association/glm_phenotype_only.txt` (12 taxa, 1 to 3 rows each, factor `env`, covariate `cov`) | Factor and covariate together (BLUEs at covariate 0), a taxon with no value of one trait, missing values in a trait and a covariate |

TASSEL 5 lists a phenotype's attributes of one type in the iteration order of a Guava `HashMultimap` of
column indices, so two `data` columns at indices 4 and 3 come out as 4 then 3 (`y2` before `y1`). TASSEL 6
keeps the file's order. The synthetic input puts its two traits at indices 1 and 2, where TASSEL 5 also keeps
file order, so every golden compares row for row.

The `mdp_traits` golden is the one place TASSEL 6 departs on purpose (departure D5 in the
`[glm-departures]` section of `MANIFEST`): with no residual degrees of freedom TASSEL 6 writes `NaN`
for `F`, `p`, and `errorMS`, where TASSEL 5 writes the sign of its sweep's rounding noise.
