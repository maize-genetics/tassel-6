# glm-synthetic

TASSEL 5 `v5.2.98` `FixedEffectLMPlugin` run on small hand-built inputs in `testdata/association/`.
The inputs are static files (a fixed-seed Python script built them once); the generator is not
part of the repository, so edit the files only through a new issue. Every report is kept rounded
(`*.txt`, compared with `format`) and at full precision (`*.full.txt`); the commands and
tolerances are in `testdata/golden/MANIFEST`.

## Genotypes: `glm_synthetic.hmp.txt` (40 taxa `S01`..`S40`, 10 sites)

| Site | What it exercises |
|------|-------------------|
| `common_biallelic` | A/C, three classes with a heterozygote: the ordinary `add_F`/`dom_F` path |
| `no_hets` | two homozygous classes only: `dom_F` is `NaN`, one marker df |
| `multiallelic` | A/C/G with AA, CC, GG, AC, AG: five classes, marker df 4; `-biallelicOnly` drops it |
| `indel` | `+`/`-` calls, a non-ACGT allele coding |
| `monomorphic` | one class: `marker_F` and `p` are `NaN`, `marker_df` is 0; `-biallelicOnly` and `-minClassSize 3` drop the site |
| `rare_het` | a heterozygous class of 2 taxa: `-minClassSize 3` sets it missing, leaving two classes (marker df 1) |
| `rare_homozygote` | a homozygous class of 2 taxa: `-minClassSize 3` sets it missing for `y1`, leaving two classes |
| `missing_calls` | 8 `NN` calls: observations dropped per site, so error df differs from other sites |
| `mostly_missing` | 31 `NN` calls: tiny error df (3 for `y1`, 1 for `y2`) and a near-saturated model |
| `two_class_het` | AA and AG only: a heterozygote without the other homozygote |

## Phenotypes

| File | Roles | What it exercises |
|------|-------|-------------------|
| `glm_synthetic_traits.txt` | `y1`, `y2` data; `grp` factor (east/west/north); `c1`, `c2` covariates | the full model; `c1` is missing for `S07` and `y2` for `S04`, `S13`, `S31`, so each trait has its own observation set |
| `glm_synthetic_factor.txt` | `y1` data; `grp` factor | a factor and no covariates |
| `glm_synthetic_collinear.txt` | `y1` data; `c1`, `c2`, `c3 = c1 + c2` covariates | an exactly collinear covariate: TASSEL 5 sweeps it out, and the error df shows it (34 at an ordinary three-class site: 40 taxa less the mean, `c1`, `c2`, and two marker df) |

## Runs

`model`, `min-class-size-3`, `biallelic-only`, and `append-add-dom` use `glm_synthetic_traits.txt`;
`factor-only` and `collinear` use the other two files. Each run has `*.stats` and `*.effects`
reports, rounded and `.full`.
