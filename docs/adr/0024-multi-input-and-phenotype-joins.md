# ADR-0024: Multi-input operations (`DataSet`) and phenotype join semantics

- Status: Proposed
- Date: 2026-11-22

## Context

Epic #155 needs three joins: phenotype with phenotype (intersect and union, #162), and genotype with
phenotype (#163). Each takes more than one table. Two things block that today.

- **One input per operation.** `Operation<I, O>` has a single `I`, and `RunRequest.input` is one
  `DataHandle?`. `Operation`, `Engine`, and the DSL are frozen (ADR-0002, AGENTS.md constraint 4).
- **No agreed join semantics.** TASSEL 5 (`PhenotypeBuilder.joinPhenotypes` and `mergeTwoPhenotypes`,
  pinned at `v5.2.98` by ADR-0009) joins by comparing only the alphabetically first shared factor,
  treats shared covariates as keys and then fails with "both phenotypes have values", treats factor
  `?` as an ordinary level, and sorts its output. ADR-0023 chose a clean redesign of everything
  except file formats, and says joins get fresh semantics and hand-written expected outputs.

TASSEL 5 solved multi-input with its GUI tree: the user selects several data nodes, and a plugin
receives them as one `net.maizegenetics.plugindef.DataSet` of `Datum` items. The CLI mirrors it with
`-union` and `-intersect` in `TasselPipeline`.

## Decision

### Multi-input: a `DataSet` value

- **Type.** `net.maizegenetics.tassel.dataset.DataSet(items: List<Any>)` is a plain core value in
  `commonMain`. It keeps its items in order and offers typed accessors, `phenotypes():
  List<PhenotypeTable>` and `genotypes(): List<GenotypeTable>`, each in `DataSet` order. Items are
  immutable core values; a `DataSet` inside a `DataSet` is flattened on construction. Accessors do
  not validate counts; each operation checks its own arity and throws `IllegalArgumentException`
  naming what it needs and what it got.
- **Name.** `DataSet` is TASSEL 5's word for "the selected data", so TASSEL 5 users and the legacy
  translator already speak it. Rejected names are under Alternatives.
- **Operations.** A join is `operation<DataSet, PhenotypeTable>(...)` or
  `operation<DataSet, GenotypePhenotype>(...)`. No change to `operation`, `Operation`, or the
  parameter DSL.
- **Engine.** `LocalEngine.put(DataSet(...))` gives one `DataHandle`, and `RunRequest.input` carries
  it. `LocalEngine.put` already takes `Any`, and `run` checks the input's type against the
  operation's, so a `DataSet` is just another type.
- **Descriptors.** `OperationDescriptor.inputType` is the string `DataSet`. That is how a front end
  learns that an operation wants several selected items. No descriptor field is added.
- **Front ends.** The CLI makes `--input` repeatable for operations whose input type is `DataSet`
  and builds the `DataSet` from the files in order (#166). The app builds it from a multi-selection,
  in selection order (#167). Both are front ends, not frozen contracts.
- **Why this is additive and not a contract change.** `I` was already any `Any`; `RunRequest` and
  `Engine` are untouched; no existing operation, descriptor, or serialized form changes meaning.
  AGENTS.md constraint 4 protects those signatures, and none is edited.

### Join keys

- The key is the **taxon** plus **every factor attribute present in both inputs**.
- Taxa match on exact, case-sensitive text. The taxa column's name is ignored for matching.
- Factors match by **name** (exact) and by level **label**, never by level code, which differs
  between tables.
- Covariate and data attributes are never keys, even when both inputs have one with that name; they
  fall under Collisions below.
- A **missing key value** (a factor cell `?`) matches nothing, not even another missing value. A row
  with one is dropped by an intersect, and kept unmatched by a union.
- A name that is a factor on one side and another type on the other is not a key; it is a
  collision.

### Intersect and union

- **Intersect** is an inner join on the full key: it keeps only matched row pairs.
- **Union** is a full outer join: it keeps every matched pair, every unmatched left row, and every
  unmatched right row. The cells the other side would have supplied are missing. A key factor on an
  unmatched row takes that row's own value; the taxon likewise.
- Joining one input with nothing else is an error (a `DataSet` needs at least two phenotypes).

### Row and column order

- **Rows.** Left rows in their order, each followed by its matches in right order. Union then
  appends unmatched right rows in right order.
- **Columns.** The left columns in their order, then the right side's columns that are new. There
  is one taxa column, named after the left one, and each key factor appears once, in its left
  position. Output is never sorted. This replaces TASSEL 5's sorted output, so expected tables are
  written by hand and ordered.

### Many-to-many keys

Replicates make a key repeat. One-to-many is allowed: a left row meets each of its several right
rows (or the reverse). A key that repeats on **both** sides is ambiguous, and the `multiple`
parameter decides:

- `error` (default): fail, naming the keys, as taxon plus factor levels. Long lists are cut after
  ten keys with a count of the rest.
- `cartesian`: emit every left-right pair, left-major, in the row order above.

A key with a missing value is never ambiguous, since it matches nothing.

### Collisions

A **shared non-key name** is a name that is not a key and appears in both inputs. The `collision`
parameter decides:

- `error` (default): fail, listing the names.
- `suffix`: rename each side to `name.1` and `name.2`, numbered by input position in the
  `DataSet` (1-based, among its phenotypes). The result keeps each in its own column position.
  If a renamed column would equal an existing name, that is an error.
- `coalesce`: one column, at the left position, taking the non-missing value. If both values are
  non-missing and differ, that is an error naming the key and the column. Factors compare by label,
  numbers by `Float` equality.

A shared name whose attribute types differ (factor against covariate, say) is an error unless
`collision = suffix`; `coalesce` cannot merge across types. Parameter spelling follows the DSL's
enum rule (constant names matched exactly); #162 fixes the spelling and the operation ids, which
this ADR proposes as `intersect-phenotypes` and `union-phenotypes`.

### N-way joins

A join over more than two phenotypes is a **left fold over `DataSet` order**: the result of joining
the first two is joined with the third, and so on. The key at each step is the taxon plus the
factors present in both the accumulated result and the next input. `multiple` and `collision` apply
at each step. For `suffix`, collisions are resolved first over the original inputs, so a name held
by inputs 1 and 3 becomes `name.1` and `name.3` whatever the fold does, and the result does not
depend on how the fold is grouped.

### Genotype and phenotype

- **Value.** `GenotypePhenotype(genotypes: GenotypeTable, phenotype: PhenotypeTable, genotypeRow:
  IntArray)`. Phenotype rows drive the result. `genotypeRow[i]` is the genotype-table taxon index
  for phenotype row `i`, or `-1` when its taxon has none.
- **Input.** A `DataSet` with exactly one genotype table and exactly one phenotype table, in any
  order. Taxa match exactly, case-sensitively. A taxon repeated in the genotype table is an error.
- **Intersect** keeps the phenotype rows whose taxon is in the genotype table, so no `genotypeRow`
  entry is `-1`.
- **Union** keeps every phenotype row, with `-1` for those with no genotype, and then appends one
  all-missing phenotype row for each genotype-only taxon, in genotype order. The appended rows have
  missing factor and numeric values, and their taxon names come from the genotype table.
- **No copying.** The genotypes are held by reference and read through the site-block API
  (ADR-0012). `GenotypePhenotype` is a new type, not a `GenotypeTable`, so the frozen
  `GenotypeTable` contract is untouched. The phenotype table is new (selected or extended) and
  follows ADR-0023's immutability.
- **Out of scope.** Joining two genotype tables by taxa or sites (TASSEL 5
  `UnionAlignmentPlugin`, `IntersectionAlignmentPlugin`) is a different operation and needs its own
  ADR or task.

### Deviations from TASSEL 5

| TASSEL 5 | TASSEL 6 | Why |
|---|---|---|
| Output rows and columns sorted | Defined order: left, then right | Predictable results; sorting hid which side a row came from |
| Only the alphabetically first shared factor is compared | Taxon plus every shared factor | Matching on one factor of several merges unrelated observations |
| A shared covariate is a key, then "both phenotypes have values" fails | Covariates are never keys; handled by `collision` | The error was the only outcome for a legitimate case |
| `?` is a level, so missing keys match each other | A missing key matches nothing | A missing value is not evidence that two rows belong together |
| Many-to-many keys join silently (or by accident of order) | `multiple = error` by default | An ambiguous join should be a choice, not a surprise |
| Only two inputs, joined by pipeline order | N-way left fold over `DataSet` order | One rule for any number of inputs |

## Consequences

- The frozen `Operation`, `Engine`, `RunRequest`, and DSL stay as they are, and every join is
  declared once like any other operation, so the CLI, app, MCP tools, and bindings get it from the
  spec (ADR-0002).
- Front ends must learn that `inputType = DataSet` means "several items" (#166, #167). Until then,
  such operations stay unregistered, as the epic plans.
- A `RemoteEngine` (ADR-0004) cannot build a `DataSet` from several handles without a new engine
  operation. That would be an additive `Engine` method, so it needs its own ADR when the remote
  engine arrives. `LocalEngine` can offer a convenience that bundles existing handles; that is not
  on the `Engine` interface.
- Because argument order is item order, left and right matter for column and row order. Front ends
  show the order they send.
- Join results get hand-written expected outputs (ADR-0023), not TASSEL 5 goldens. A TASSEL 5
  output may serve as a row-unordered cross-check where both semantics agree.
- The legacy translator maps `-union` and `-intersect` to these operations later; that mapping is
  not designed here.
- Follow-ups: `DataSet` and phenotype joins (#162), genotype-phenotype (#163), CLI (#166), app
  (#167).

## Alternatives considered

- **`Operation<List<Any>, O>` or a vararg input.** Changes the DSL's single-`I` contract and loses
  the typed accessors; an ADR superseding ADR-0002 for no gain.
- **`RunRequest.inputs: List<DataHandle>`.** A frozen-contract change across the engine, every
  transport, and every front end, to solve what one wrapper type solves.
- **A parameter naming the second input.** Parameters are text (ADR-0002), and a table cannot be
  sent as text; it would also need handle names inside parameter strings.
- **Other names (`DataBundle`, `Selection`, `InputSet`).** `DataSet` is TASSEL 5's term; the others
  are new vocabulary.
- **Match on the first shared factor, as TASSEL 5 does.** Rejected for the reasons in the table.
- **Join key chosen by the caller.** Flexible, but needs a text-valued list parameter and gives the
  wrong answer silently when the caller forgets a factor. Using every shared factor is the safe
  default; a `keys` override can be added later without changing this decision.
- **Last-wins or first-wins collisions.** Silent data loss; `coalesce` keeps the information and
  fails when the inputs disagree.
