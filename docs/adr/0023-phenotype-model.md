# ADR-0023: Phenotype data model and TASSEL 5 file compatibility

- Status: Proposed
- Date: 2026-11-22

## Context

GLM, MLM, LD, and the teaching release (ADR-0013) all need phenotype data. Until this ADR,
TASSEL 6 had only a write-only `table.PhenotypeTable` that PCA used for its output, holding each
column as a `List<Any?>`. Epic #155 builds the phenotype foundations: model, file formats, joins,
filters, transforms, summaries, and F, t, and chi-square distributions. This ADR fixes the model
and the reading rules. ADR-0024 (#161) will fix joins and multi-input operations.

TASSEL 5 (`net.maizegenetics.phenotype`, pinned at `v5.2.98` by ADR-0009) stores a phenotype as a
list of attributes:
- `NumericAttribute`: a `float[]` plus a missing `BitSet`.
- `CategoricalAttribute`: level codes into a `TreeSet` of labels.
- `TaxaAttribute`: one taxon per observation.

Its implementation has behavior that users would not choose and analyses cannot rely on:
- The `BitSet` and NaN can disagree. A `<Header>` file leaves NaN cells unflagged, and transforms
  create NaN and -Infinity without flagging them.
- `CategoricalAttribute` never marks a value missing, so `?` is an ordinary level.
- `TaxaAttribute` ignores the name it is given and is always `Taxa`.
- `areTaxaReplicated()` is inverted.
- Attribute lookup is by identity, which breaks on filtered phenotypes.
- The same file convention means different things in different formats. In `<Phenotype>`, `NA`,
  `NaN`, and `.` are missing and any other non-number is an error. In `<Trait>`, `-99` and `-999`
  are missing and any unparseable token is *silently* missing.

The maintainer chose, for epic #155, a clean redesign. TASSEL 5 is the reference for the file
formats only: TASSEL 5 files must load with the same values, and TASSEL 6 output must load in
TASSEL 5. Everything else may be defined fresh.

## Decision

### Model

- **Package and types.** `net.maizegenetics.tassel.phenotype` holds `PhenotypeTable` and a sealed
  `PhenotypeAttribute` with three kinds: `TaxaAttribute`, `NumericAttribute`, `FactorAttribute`.
  `AttributeType` keeps TASSEL 5's names (`taxa`, `factor`, `covariate`, `data`), because they are
  the words written in the file. The class moves here from `table`; its file bytes do not change.
- **Rows and replicates.** One row is one observation. A taxon may have many observations;
  replicates and environments are told apart by factor attributes, as in TASSEL 5.
- **Table invariants.** A table has exactly one taxa attribute and unique attribute names, and
  every attribute has the same number of observations. The constructor rejects anything else.
  TASSEL 5 checks none of these.
- **Immutability.** Tables and attributes are immutable. `selectRows` and `selectAttributes`
  return new tables and copy the values they keep; phenotypes are small next to genotypes, so
  copying is cheaper to reason about than views. `selectAttributes` keeps the taxa attribute even
  when it is not named.
- **Numeric values are `Float`, read as `Double`.** Numeric attributes store a `FloatArray`, as
  TASSEL 5 does, so values match TASSEL 5 bit for bit, and so do the inputs later analyses see.
  Analyses read them through `value(observation): Double`. AGENTS.md's "`Double` throughout"
  applies to the arithmetic, not to this storage.
- **Missing numbers are NaN, and only NaN.** There is no separate missing mask. Constructors store
  any infinite value as NaN, so "no usable value" has one representation.
- **Missing factor values are code -1**, written as `?`, TASSEL 5's spelling, so TASSEL 5 can
  read the file back; it sees `?` as a level. `?` cannot be a level.
- **Factor levels.** Levels are unique and sorted by `String.compareTo`, as TASSEL 5 sorts them
  (`"10"` before `"2"`). Selecting rows keeps only the levels still used.
- **The taxa attribute keeps its name** from the file or the caller. `Taxa` is used when nothing
  names it, as for PCA output.

### `<Phenotype>` format

- **Reading.** `PhenotypeFiles.read(source, options, progress)` returns a `PhenotypeRead` holding
  the table and any warnings, as `VcfRead` does (ADR-0017 `Source`, ADR-0019 progress).
- **Header lines.** The tag matches `<phenotype` in any case. Types match by TASSEL 5's
  case-insensitive prefixes `tax`, `fac`, `cov`, and `dat`.
- **Cells** are split on runs of Java `\s` whitespace, and numbers are parsed exactly as Java's
  `Float.parseFloat` parses them (`JavaFloat`, bit-exact on every target).
- **More forgiving than TASSEL 5:** lines are trimmed, a byte-order mark is dropped, and blank
  lines and lines starting with `#` are skipped.
- **Stricter than TASSEL 5:**
  - a type other than the four above is an error;
  - a repeated attribute name is an error;
  - zero or several taxa columns is an error;
  - every error names the line.
- **Missing tokens.** The defaults are `NA` and `NaN` in any case, `.`, and `?`, in numeric and
  factor columns alike. `-999` is a number in a `<Phenotype>` file, as in TASSEL 5.
  `PhenotypeReadOptions.missingTokens` replaces the defaults with an exact set.
- **Writing.** `writeDelimited` writes TASSEL 5's `PhenotypeUtils.write` layout. Numbers are
  written with Java's `Float.toString`, a missing number as `NaN`, and a missing level as `?`. The
  writer refuses names, taxa, or levels that would not read back as themselves: whitespace, or a
  level spelled like a default missing token.

### Legacy formats

`<Trait>`, `<Covariate>`, `<Factor>`, and `<Header name=x>` files are read by the same
`PhenotypeFiles.read` (issue #157). The format is detected from the first tag after the
byte-order mark, blank lines, and `#` lines: `<Phenotype` is the format above, `<Trait`,
`<Covariate`, `<Factor`, or `<Header` is legacy, and anything else is an error naming the accepted
tags. TASSEL 5's prefixes (`trai`, `cov`, `fac`, `head`) and case-insensitivity are kept.

- **`<Trait>`** names the columns after the taxon. A `<Covariate>` line makes them covariates,
  otherwise they are data; a `<Factor>` line makes them factors. A `<Covariate>`, `<Factor>`, or
  `<Header>` without a `<Trait>` line is an error. The taxa column is named `Taxa`.
- **`<Header name=x>`** lines reshape the file to long format. Each header line becomes a factor
  column `x`. Each data line gives one observation per distinct combination of header levels, in
  first-appearance order, and traits are in first-appearance order too. A (trait, combination)
  pair the file does not give is missing. TASSEL 5 sorts rows and columns, so the goldens are
  compared with `rows-unordered`.
- **Missing tokens.** The defaults are `-999`, `-99` (written exactly, so `-999.0` is a number),
  `NA` and `NaN` in any case, `.`, and `?`, in numeric and factor columns alike. Any other
  unparseable token is an error naming the line and column. `PhenotypeReadOptions.missingTokens`
  replaces the defaults with an exact set.

### Parity scope (refines ADR-0009 for phenotypes)

- **Format goldens.** TASSEL 5 exports of its own test inputs are goldens under
  `testdata/golden/phenotype-*`. They are compared at a new `phenotype` tolerance: the tag and types
  lines identical, the names line identical except the taxa column's name, and every data cell
  identical as text. Identical text means identical floats, because both sides write
  `Float.toString`. The tolerance allows `rows-unordered` for formats whose row order TASSEL 6
  defines differently.
- **What else gets parity.** Transform numerics on valid inputs and the distributions get goldens
  in their own tasks.
- **Joins and new behavior** (filters, summaries, missing-value handling) get hand-written expected
  outputs, not TASSEL 5 goldens.

### Deviations from TASSEL 5 (so far)

| TASSEL 5 | TASSEL 6 |
|---|---|
| Taxa column always named `Taxa` | Name kept from the file |
| Missing `BitSet` alongside NaN, which can disagree | NaN is the only missing marker |
| `Infinity` stored and not missing | Read as missing, with a warning |
| Factors have no missing value; `?` is a level | Code -1, written `?`; `NA`, `NaN`, `.`, and `?` read as missing |
| Leading whitespace, BOM, blank lines, or `#` lines break a `<Phenotype>` file | Trimmed or skipped |
| Duplicate names and unknown types accepted or misread | Errors |
| Hexadecimal floats (`0x1p3`) accepted | Not accepted; reported as not a number |
| `<Trait>`: an unparseable token is silently missing | Error naming the line and column |
| `<Header>`: `-999` and `-99` are numbers, and a NaN cell is unflagged | Missing, as in `<Trait>` files |
| `<Trait>`/`<Header>` factor cells are always levels | `NA`, `NaN`, `.`, `?`, `-999`, and `-99` are missing |
| `<Header>` rows and columns sorted | First-appearance order |
| `<Header>`: short headers fail with an index error, long ones are truncated, a repeated (trait, level) overwrites | Errors |
| Tags may follow data; `<Factor>` with `<Covariate>` or `<Header>` ignores one tag | Tags must precede data; those combinations are errors |
| `TransformDataPlugin` keeps an unflagged NaN and `-Infinity` (log of 0 or of a negative, a power with no real result) | Any NaN or infinite result is missing, and the result's warnings count them per trait |
| `TransformDataPlugin` standardizes with a `Float` mean and `Float` arithmetic, so a standardized value can be off by up to ~8e-7 | `Double` throughout, then narrowed to `Float`; goldens are compared within 1e-6 relative or absolute |
| A trait name that is not a data attribute is silently ignored; `-factors` without `-standardize` is ignored | Errors naming the attribute and its type |
| A row with a missing factor value is left out of every group and keeps its unstandardized value | The row's result is missing, and the warning says why |
| `-standardize` of one value or of a constant trait gives NaN, kept unflagged | Missing, with a warning |
| `Math.pow(NaN, 0)` turns a missing value into 1 | A missing value stays missing |

## Consequences

- One representation of missing values makes filters, summaries, and analyses simpler and
  removes a whole class of TASSEL 5 bugs (NaN without a flag, flags without NaN).
- Exported text still matches TASSEL 5 for every file TASSEL 5 reads correctly, because a value
  TASSEL 5 leaves unflagged is written as `NaN` by both.
- `Float` storage keeps later GLM and MLM inputs identical to TASSEL 5's. Analyses widen to
  `Double` once, on read.
- A factor cell spelled `NA` is missing in TASSEL 6 but a level in TASSEL 5, so an analysis on such
  a file can differ from TASSEL 5. This is intended.
- `formatCell(Any?)` and the untyped `PhenotypeTable.Attribute` are gone. PCA, the CLI, and the app
  use the typed attributes.
- Follow-ups: filters (#158), summary (#159), transforms (#160, done), and joins
  under ADR-0024 (#161, #162, #163).

## Alternatives considered

- **Keep a missing mask beside NaN, as TASSEL 5 does.** Two sources of truth disagree; that is
  most of TASSEL 5's missing-value bugs.
- **`Double` storage.** Simpler arithmetic, but exports and GLM inputs would differ from TASSEL 5
  in the last digits, and goldens would need tolerances where text identity is free today.
- **Bug-for-bug parity.** Rejected by the maintainer for this epic: users would inherit
  first-key-only joins and unflagged missing values.
- **Levels in first-appearance order.** Matches some tools, but TASSEL 5 and R sort, and level
  order decides the reference level in later linear models.
