# Porting ledger

Tracks every TASSEL 5 class or plugin that TASSEL 6 ports, bridges, redesigns,
or retires. Update the row in the same PR that changes its status.

Paths for TASSEL 5 are relative to `src/main/java/net/maizegenetics/` in the
TASSEL 5 repo (`reference/tassel5/` locally). Pinned reference version: see
ADR-0009.

**Strategy**
- `port`: reimplemented natively in `tassel-core` with golden parity tests.
- `redesign`: replaced by a TASSEL 6 concept, with behavior kept where it matters.
- `bridge`: run through `legacy/tassel5-bridge` (JVM only).
- `retire`: intentionally not carried forward (give the reason).

**Status**: `planned`, `in progress`, `parity` (golden tests pass on all core
targets), `bridged`, `retired`.

## PoC (`6.0.0-poc.1`)

| TASSEL 5 source | TASSEL 6 target | Strategy | Golden test | Status | Issue |
|-----------------|-----------------|----------|-------------|--------|-------|
| `dna/snp/NucleotideAlignmentConstants.java` | `genotype.Nucleotide`, `Allele`, `Genotype` | port | encoding tables vs TASSEL 5 `v5.2.98` (`NucleotideTest`) | parity | #10 |
| `util/SuperByteMatrix.java`, `dna/snp/genotypecall/ByteGenotypeCallTable.java` | `genotype.GenotypeTable`, `SiteBlock`, `InMemoryGenotypeTable` (ADR-0012) | redesign | via all I/O goldens | in progress | #10 |
| `dna/snp/genotypecall/AlleleFreqCache.java`, `dna/snp/genotypecall/Stats.java` | `genotype.AlleleStats` | port | cases vs TASSEL 5 (`AlleleStatsTest`); all 3,093 sites of `summary/sites.full.txt` (`GenotypeSummaryTest`) | parity | #10 |
| `taxa/TaxaList*.java`, `taxa/Taxon.java`, `dna/map/PositionList*.java`, `dna/map/Position.java`, `dna/map/Chromosome.java` | `taxa.TaxaList`, `taxa.Taxon`, `position.PositionList`, `position.Position`, `position.Chromosome` | redesign | names and ordering vs TASSEL 5 (`PositionTest`); the HapMap golden (`HapMapTest`) | parity | #10 |
| `dna/snp/io/BuilderFromHapMap.java` | `io.HapMap.read` | port | loads `mdp_genotype.hmp.txt` as 281 x 3,093 and writes back to the golden (`HapMapTest`) | parity | #13 |
| `dna/snp/ExportUtils.java` (HapMap writer) | `io.HapMap.write` | port | `golden/hapmap/mdp_genotype.hmp.txt`, exact (`HapMapTest`) | parity | #13 |
| `dna/snp/ExportUtils.java` (`writeToVCF`, plain text) | `io.Vcf.write`, internal `vcfSiteRecord` | port | `golden/vcf-export/mdp_genotype.vcf`, `include-list.vcf`, and `vcf_export_edges.vcf`, exact (`VcfWriteTest`) | parity | #127 |
| `dna/snp/ExportUtils.java` (`writeToVCF` of a `.vcf.gz`, via `Utils.getBufferedWriter` / `GZIPOutputStream`). TASSEL 5 has no BGZF writer | `io.BgzfWriter` at one fixed DEFLATE level, parallel by default. `Vcf.write` is still plain text | redesign | #118 and #129 payloads round-trip through `BgzfReader` with the level-0 ISIZE list, and the compressed tutorial is 983,502 bytes with a recorded SHA-256 (`BgzfDeflateTest`). One worker and N workers write that same digest (`BgzfParallelTest.parallelEncodeIsIndependentOfWorkerCount`). `Inflater` accepts each tutorial payload on the JVM (`BgzfWriterJvmTest`). Stored blocks remain the level-0 check (`BgzfWriterTest`) | parity | #128, #131, #132 |
| `dna/snp/io/BuilderFromVCF.java` (`HeaderPositions`, `ProcessVCFBlock`, the `GT` field only) | `io.Vcf.read`, byte-level `VcfHeader` and `VcfRecordParser`; allele and `GT` rules in `io.variantsOf`, `io.splitVariants`, `io.allelesOf`, and `io.call` | port | loads `mdp_genotype_w_depth.vcf` as 281 x 3,093 and writes back to `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down, exact (`VcfTest`, `VcfEndToEndTest`) | parity | #14, #83, #121 |
| `util/Utils.java` (`getBufferedReader`, the `.gz` wrap) | gzip `RawSource` (BGZF members one block at a time, plain gzip through a 32 KiB window) | port | `htslib/tutorial/mdp_genotype_w_depth.vcf.gz` and a `GZIPOutputStream` of the tutorial VCF versus `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down, exact (`GzipSourceTest`, `GzipStreamJvmTest`, `VcfEndToEndTest`); `htslib/variants/plain.vcf.gz` and `concatenated.vcf.gz` versus `source.vcf` (`GzipSourceTest`) | parity | #84, #85 |
| `analysis/data/FileLoadPlugin.java` ("Make Best Guess" for `.hmp.txt`, `.hmp.txt.gz`, `.vcf`, `.vcf.gz` only). TASSEL 5's GUI loads whole files; it has no region open | `app.data.GenotypeFile.readGenotypes` for HapMap and plain VCF. A `.vcf.gz` is opened with its `.csi` or `.tbi` and one region is read through `IndexedVcf` (`IndexedOpen`, `RegionPane`). A `.bcf` is not guessed as HapMap | port | suffix dispatch in `GenotypeFilesTest`; desktop `.vcf` and bgzipped `.vcf.gz` of `mdp_genotype_w_depth` are the same 281 x 3,093 table (`DesktopGenotypeFilesTest`). Pairing and region validation in `IndexedOpenTest` and `RegionLoadTest`. Desktop `.csi`-then-`.tbi` lookup, and `chr3B:599000000-601000000` of `htslib/wheat/wheat.vcf.gz` matching `IndexedVcf.read` (`DesktopGenotypeFilesTest`). Other TASSEL 5 file types are not guessed. HapMap text named `.bcf` fails with the BCF magic message (`GenotypeFilesTest`) | parity | #92, #93, #126 |
| none; TASSEL 5 has no region queries and no CSI or TBI | `io.IndexedVcf` | redesign | htslib 1.25 `.record.txt` regions equal `Vcf.read` of the header plus those lines, on CSI and TBI (`IndexedVcfTest`) | parity | #90 |
| none; TASSEL 5 has no CSI reader | `io.index.readGenomicIndex` with caller-supplied names for a CSI whose `l_aux` is 0 | redesign | the #118 BCF CSI files against `bcftools index -s` and the region notes (`BcfCsiIndexTest`) | parity | #120 |
| none; TASSEL 5 has no BCF reader | `io.readBcfHeader` | redesign | the #118 BCF headers against `bcftools view` dictionaries, including an `IDX` gap and PASS at 0 (`BcfHeaderTest`) | parity | #121 |
| none; TASSEL 5 has no BCF reader | `io.readBcfRecord` | redesign | the #118 edge, boundary, and no-sample records against `VcfRecordParser` on htslib's text (`BcfRecordTest`); direct GT against the text path, including truncation before a call error (`BcfGenotypeTest`) | parity | #122, #145 |
| none; TASSEL 5 has no BCF reader | `io.Bcf.read` | redesign | the tutorial BCF writes `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down, and the #118 BCF files equal `Vcf.read` of htslib's text, including a missing EOF warning and a truncated record (`BcfReadTest`) | parity | #123 |
| none; TASSEL 5 has no BCF writer | `io.Bcf.write` | redesign | inflated bytes and the uncompressed BGZF block list of `htslib/bcf/golden/mdp_genotype.bcf`, `include-list.bcf`, and `vcf_export_edges.bcf` (`BcfWriteTest`). On-disk members use the #131 DEFLATE level | parity | #130, #131 |
| none; TASSEL 5 has no BCF region query | `io.IndexedGenotypes`, `io.IndexedBcf` | redesign | #118 regions equal `IndexedVcf.read` of the sibling `.vcf.gz`; whole-contig reads equal `Bcf.read`; ADR-0022 no-match cases, including wheat `chr6D` as an empty table (`IndexedBcfTest`). `chr3B:1000` reads the 620-byte header, the 14,392-byte index, and only the BGZF blocks of its one chunk | parity | #124 |
| none; TASSEL 5 has no BCF input and no `--region` | `cli.readGenotypes` for a whole BCF and one CSI region | redesign | `filter-sites` of `htslib/bcf/tutorial/mdp_genotype_w_depth.bcf` writes `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down; `chr3B:599000000-601000000` of `wheat.bcf` matches the same region of `wheat.vcf.gz` (`BcfCliTest`) | parity (JVM, the CLI's only target) | #125 |
| `dna/snp/ExportUtils.java` (`-export` suffix) and `util/Utils.java` (`addSuffixIfNeeded`, case-insensitive `.hmp.txt`, `.vcf`, `.vcf.gz`). TASSEL 5 has no BCF export | `cli.outputParts` / `cli.outputPath`: a genotype `--output` of `.hmp.txt`, `.vcf`, `.vcf.gz`, or `.bcf` selects that writer. `.vcf.gz` and `.bcf` are BGZF at the fixed DEFLATE level, not `GZIPOutputStream`. There is no `convert` operation (ADR-0022) | redesign | `filter-sites` of `tutorial/mdp_genotype.hmp.txt` writes `golden/vcf-export/mdp_genotype.vcf`, exact; `.vcf.gz` inflates to that text; `.bcf` reads back equal through `Bcf.read` (`GenotypeExportCliTest`). HapMap output is unchanged (`VcfEndToEndTest`, `KinshipEndToEndTest`) | parity (JVM, the CLI's only target) | #133 |
| none; TASSEL 5's GUI has no BCF open | `app.data.GenotypeFile.readGenotypes` for a whole BCF and `IndexedOpen` for one CSI region | redesign | `.bcf` pairs with `.bcf.csi` only (`IndexedOpenTest`). Desktop lookup uses that `.csi` and ignores a `.tbi`. The tutorial BCF is the same 281 x 3,093 table as the plain VCF, `small.raw.bcf` matches `small.view.vcf`, and `chr3B:599000000-601000000` of `wheat.bcf` matches `IndexedVcf.read` of `wheat.vcf.gz` while the pane keeps header `chr6D` (`DesktopGenotypeFilesTest`) | parity (desktop JVM and wasmJs) | #126 |
| none; TASSEL 5's GUI has no BCF export | `app.export.writeGenotype`: HapMap, plain VCF, `Vcf.writeBgzf`, and `Bcf.write`. The browser download is a byte-array Blob | redesign | file names in `GenotypeExportTest`; each format's bytes against the core writer, with `.vcf.gz` inflated back to `Vcf.write`, in `ExportTest`; the browser copy keeps bytes 0 through 255 and a buffer past one 32,768-byte chunk (`JsBytesTest`) | parity (desktop JVM and wasmJs) | #134 |
| `dna/snp/io/BuilderFromVCF.java` depth, quality, and likelihood fields (`keepDepth`, `AD`, `DP`, `GQ`, `PL`) and `dna/snp/score/AlleleDepthBuilder.java` | none; `io.Vcf.read` reads past every `FORMAT` field after `GT` | retire | n/a; the HapMap golden carries no depth | retired | #14 |
| `taxa/TaxaListIOUtils.parseVCFHeadersIntoMap`, `taxa/Taxon` annotations (the VCF's `##SAMPLE=` headers) | none; read and dropped by `io.Vcf.read` | retire | TASSEL 6's VCF-to-HapMap output is missing the golden's 281 `##SAMPLE=` lines, so the golden is compared from the `rs#` row down (`VcfTest`, maintainer decision 2026-11-04) | retired | #14 |
| `dna/map/GeneralPosition` `WHICH_ALLELE.Reference` and `Alternate`, and its `INFO` annotations | none; `Position` has no reference, alternate, or annotation slot, so only `knownVariants` survives | retire | a VCF-backed site summary prints `N` for `Ref` and `Alt`, as #43 hardcodes (`GenotypeSummaryTest`) | retired | #14, #43 |
| `dna/snp/ExportUtils.writeToVCF` INFO, QUAL, FILTER, depth (`AD`, `DP`, `GQ`, `PL`), the reference allele, and `##SAMPLE` on VCF-to-VCF export | none; `io.Vcf.write` emits QUAL `.`, FILTER `PASS`, INFO `.`, and FORMAT `GT` because the model has no annotations or depth (#14). The `##Tassel` line is written because there is no reference allele | retire | the three `golden/vcf-export/*.vcf` files have those fixed columns and no `##SAMPLE` (`VcfWriteTest`) | retired | #14, #127 |
| `analysis/data/GenotypeSummaryPlugin.java` | `analysis.genotypeSummary`, `ops.summaryOperation` | port | all eight `summary/*` goldens, rounded and full precision (`GenotypeSummaryTest`) | parity | #15 |
| `analysis/filter/FilterSiteBuilderPlugin.java`, `dna/snp/FilterSite.java`, `analysis/filter/FilterBySites.java` | `filter.filterSites`, `filter.SiteFilter`, `filter-sites` operation | port | all six `golden/filter-sites/*.hmp.txt`, exact, through the kernel and (for `position-range-5-6`) through the spec (`FilterSitesTest`) | parity | #16, #45 |
| `analysis/filter/FilterTaxaBuilderPlugin.java`, `dna/snp/FilterTaxa.java`, `analysis/filter/FilterByTaxa.java` | `filter.filterTaxa`, `filter.TaxaFilter`, `filter-taxa` operation | port | all four `golden/filter-taxa/*.hmp.txt`, exact, through the kernel and (for `include-list` and `exclude-list`) through the spec (`FilterTaxaTest`) | parity | #16, #45 |
| `matrixalgebra/decomposition/SymmetricEigenvalueDecomposition.java` | `math.SymmetricEigenDecomposition` | port | R `prcomp` values from `PrinCompTest`, through the Gram matrix (`PrinCompOracleTest`); the PCA goldens (`PrincipalComponentsTest`) | parity | #17 |
| `analysis/distance/KinshipPlugin.java`, `analysis/distance/EndelmanDistanceMatrix.java` | `kinship` operation (Centered IBS). Default `distance.centeredIbs` is a blocked dosage product (ADR-0020); `distance.centeredIbsTassel5Order` keeps the bit-exact port | port | `golden/kinship/centered_ibs.txt` within `abs 8e-6`, `SumPk` within `abs 1e-4` (`CenteredIbsTest`); byte-exact text from the retained kernel | parity | #18, #62 |
| `taxa/distance/DistanceMatrix.java`, `DistanceMatrixBuilder.java`, `WriteDistanceMatrix.java` (`saveDelimitedDistanceMatrix`) | `distance.DistanceMatrix` and its `writeDelimited` | port | `golden/kinship/centered_ibs.txt` written from its values (`CenteredIbsTest`) | parity | #18 |
| `analysis/data/PrincipalComponentsPlugin.java`, `stats/PCA/PrinComp.java` | `pca` operation (Gram-matrix eigen), `analysis.principalComponents` | port | `golden/pca/eigenvectors*` and `golden/pca/pcs*` whole; `golden/pca/eigenvalues*` at their `rows-unordered` MANIFEST tolerance, because TASSEL 5's SVD does not sort (`PrincipalComponentsTest`, #47) | parity | #19 |
| `analysis/numericaltransform/NumericalGenotypePlugin.java` (`setAlternateMinorAllelesToMinor`), `dna/snp/GenotypeTableUtils.convertGenotypeToFloatProbability`, `analysis/numericaltransform/ImputationPlugin.java` and `ImputationByMean.java` (the `by_mean` path) | `analysis.ReferenceProbability` | port | the PCA goldens, which the taxon-mean imputation moves by 7.6% if it is wrong (`PrincipalComponentsTest`) | parity | #19 |
| `dna/snp/score/SiteScoreUtil.java`, `dna/snp/score/ReferenceProbabilityBuilder.java` | `analysis.SiteScore` | port | the 255-step round trip (`SiteScoreTest`); the PCA goldens, which need it to reach `rel 1e-9` | parity | #19 |
| `phenotype/Phenotype.java`, `phenotype/PhenotypeUtils.java` (`write`) | `phenotype.PhenotypeTable` (moved from `table` in #156) and its `writeDelimited` | redesign | `golden/pca/pcs.txt` layout and values (`PhenotypeTableTest`, `PrincipalComponentsTest`); the `phenotype-format` goldens written back byte for byte (`PhenotypeFilesTest`) | parity | #19, #156 |
| `util/TableReport.java`, `util/SimpleTableReport.java` | `table.ResultTable`, `table.ResultTableBuilder` | redesign | via the summary and PCA table goldens (#15, #19) | in progress | #38 |
| `util/TableReportUtils.java` (`saveDelimitedTableReport`) | `table.ResultTable.writeDelimited` | port | `summary/alleles.txt` written from its values (`ResultTableTest`) | parity | #38 |
| `util/DoubleFormat.java` | `table.DoubleFormat` | port | every numeric cell of `summary/*.full.txt` and `pca/*.full.txt` against its rounded twin (`DoubleFormatGoldenTest`) | parity | #38 |
| `java.lang.Float.toString`, `Double.toString` (how TASSEL writes float matrix cells and quantitative annotations) | `table.floatToString`, `table.doubleToString` | port | the kinship golden's cells and `SumPk` (`CenteredIbsTest`); 400,000 values against the JVM (`ShortestDecimalJvmTest`) | parity | #18 |
| `plugindef/PluginParameter.java`, `plugindef/AbstractPlugin.java` | operation spec DSL (ADR-0002, ADR-0018) | redesign | n/a | in progress | #11, #45 |
| `pipeline/TasselPipeline.java` (`-h`, `-export`, `-exportType`, and the per-plugin flags it builds from `PluginParameter`) | `tassel-cli`: one subcommand per registered operation, generated from the specs (ADR-0002, ADR-0006) | redesign | `golden/kinship/centered_ibs.txt` written by `tassel kinship`, exact (`KinshipEndToEndTest`); `golden/vcf/mdp_genotype_w_depth.hmp.txt` written by `tassel filter-sites` from a `.vcf` input (`VcfEndToEndTest`) | parity (JVM, the CLI's only target) | #20, #14, #133 |
| `pipeline/TasselPipeline.java` | `tassel translate` (ADR-0006) | redesign | translated commands give golden outputs | planned (stretch) | |

### Filter parameters not ported yet (#16)

The two filter rows above are at parity for every parameter the goldens exercise. These TASSEL 5
parameters are not ported. Nothing in the PoC needs them and no golden covers them.

| TASSEL 5 parameter | Plugin | Why not |
|--------------------|--------|---------|
| `siteNames` | `FilterSiteBuilderPlugin` | Site-name lists; out of scope in #16 |
| `positionList`, `bedFile`, `chrPosFile` | `FilterSiteBuilderPlugin` | Position-list and BED filters; out of scope in #16 |
| `includeSites` | `FilterSiteBuilderPlugin` | Only has a meaning with the four parameters above |
| `removeSitesWithIndels` | `FilterSiteBuilderPlugin` | Not in #16's parameter table; one predicate on `AlleleStats.hasIndel` when wanted |
| `filterName` | both | Names the TASSEL 5 `Datum` and the gzipped JSON filter definition, neither of which TASSEL 6 has |

`startChr` and `endChr` (`FilterSiteBuilderPlugin`) and `taxaList`
(`FilterTaxaBuilderPlugin`) were in this list until #45. They are now declared on the
`filter-sites` and `filter-taxa` specs with the `string` and `stringList` parameter types ADR-0018
added, so all ten golden parameter sets are reachable through the specs, the CLI, and the legacy
translator. `taxaList` takes only TASSEL 5's inline comma-separated form; its `.txt`, `.json`, and
`.json.gz` file forms are not ported, because a spec parameter does not read files and no golden
uses them.

## Phase 2: phenotypes and statistics (epic #155)

ADR-0023 records the model and the format rules. TASSEL 5 is the reference for the phenotype file
formats; join, filter, and missing-value semantics are a redesign, so those rows are tested against
hand-written expected tables rather than TASSEL 5 goldens.

| TASSEL 5 source | TASSEL 6 target | Strategy | Golden test | Status | Issue |
|-----------------|-----------------|----------|-------------|--------|-------|
| `phenotype/CorePhenotype.java`, `NumericAttribute.java`, `CategoricalAttribute.java`, `TaxaAttribute.java`, `FilterPhenotype.java` (row subsets) | `phenotype.PhenotypeTable` (`selectRows`, `selectAttributes`), `phenotype.TaxaAttribute`, `NumericAttribute`, `FactorAttribute`. `Float` storage, NaN the only missing number, factor missing code -1, sorted levels, unique names, one taxa column (ADR-0023) | redesign | model invariants (`PhenotypeTableTest`); the `phenotype-format` goldens (`PhenotypeFilesTest`) | parity | #156 |
| `phenotype/PhenotypeBuilder.java` (`importFile` detection, `importPhenotypeFile`) | `phenotype.PhenotypeFiles.read`, `PhenotypeReadOptions`, `PhenotypeRead` | port | `golden/phenotype-format/basic_phenotype_test.txt` and `earht_with_rep.txt` at the `phenotype` tolerance, and each golden written back byte for byte (`PhenotypeFilesTest`) | parity | #156 |
| `java.lang.Float.parseFloat` (how TASSEL 5 parses every phenotype cell) | internal `phenotype.JavaFloat.parse`; hexadecimal literals are not accepted | port | bit-exact against `Float.parseFloat` on 55,000 inputs, including decimals at and either side of float halfway points (`JavaFloatJvmTest`); a JVM-generated corpus on every target (`JavaFloatTest`) | parity | #156 |
| `phenotype/PhenotypeBuilder.java` (`importTraits`, `processTraits`, `processFactors`, `processTraitsAndFactors`) | `phenotype.PhenotypeFiles.read` (`LegacyPhenotypeReader`) for `<Trait>`, `<Covariate>`, `<Factor>`, `<Header name=x>`; unparseable cells are errors, `<Header>` order is first-appearance (ADR-0023) | port | `golden/phenotype-legacy/` for the tutorial traits and population structure and the `PhenotypeTests` and `CandidateTests` files at `phenotype rows-unordered` (`LegacyPhenotypeFilesTest`) | parity | #157 |
| `phenotype/PhenotypeBuilder.java` (`keepTaxa`, `removeTaxa`, `keepAttributes`, `changeAttributeType`, `removeMissingObservations`), `analysis/filter/FilterTraitsPlugin.java` | `phenotype.filterPhenotype`, `PhenotypeFilter`, and the `filter-phenotype` spec (declared, registered by #166). `excludeLast` keeps TASSEL 5's `excludeLastTrait` name but drops the last non-taxa column; factor-level and numeric-range filters are new | redesign | hand-written expected tables (`FilterPhenotypeTest`) | ported | #158 |
| none; TASSEL 5 has no phenotype summary (conventions from `TransformDataPlugin.meanStdDev`, `GenotypeSummaryPlugin`) | `analysis.summarizePhenotype`, `PhenotypeSummary`, and the `phenotype-summary` spec (declared, registered by #166) | redesign | hand-computed values (`PhenotypeSummaryTest`) | ported | #159 |
| `analysis/numericaltransform/TransformDataPlugin.java` (`powerTransform`, `logTransform`, `standardize`, `subsets`) | `analysis.transformPhenotype`, `PhenotypeTransform`, `LogBase`, and the `transform-phenotype` spec (declared, registered by #166). Power wins over log; math in `Double` then narrowed to `Float`; NaN and infinite results are missing with warnings (ADR-0023) | port | `golden/phenotype-transform/` (log natural and base 10, power 0.5, standardize, standardize by `rep`) at `phenotype rel 1e-6 abs 1e-6`; hand-computed cases (`TransformPhenotypeTest`) | ported | #160 |
| `phenotype/PhenotypeBuilder.java` (`joinPhenotypes`, `mergeTwoPhenotypes`), `analysis/data/IntersectionAlignmentPlugin.java`, `UnionAlignmentPlugin.java` | `join.joinPhenotypes`, `JoinOptions`, the `join-phenotypes` spec (declared, registered by #166; `type = union` is how the translator maps `UnionAlignmentPlugin`), and `dataset.DataSet` (ADR-0024). Genotype-phenotype joins are #163 | redesign | hand-written fixtures in `testdata/phenotype/joins/` with expected files, row order included (`JoinPhenotypesTest`); `golden/phenotype-join/` (`-intersect`, `-union`) at `phenotype rows-unordered` | ported | #161, #162 |
| `phenotype/GenotypePhenotype.java`, `GenotypePhenotypeBuilder.java` | `join.GenotypePhenotype` (phenotype rows plus a `genotypeRow` map, genotypes held by reference), `join.alignGenotypePhenotype`, and the `join-genotype-phenotype` spec (declared, registered by #166). Intersect drops ungenotyped observations; union keeps them (`-1`) and appends all-missing rows for genotype-only taxa (ADR-0024) | redesign | hand-written cases, including unsorted genotype taxa, replicates, and a counting `GenotypeTable` stub (`GenotypePhenotypeTest`) | ported | #163 |
| `stats/linearmodels/LinearModelUtils.java` (`Ftest`) | `stats.fTestPValue` | port | F = 0, F = inf, non-integer df, and the throw on F < 0 (`SpecialFunctionsTest`) | parity | #164 |
| Colt 1.2.0 `cern.jet.stat.Gamma` (`gamma`, `logGamma`, `incompleteBeta`, `incompleteGamma`, `incompleteGammaComplement`; dependency replaced) | `stats.SpecialFunctions` | port | 3,358-row grid of Colt values from the pinned standalone's jar, rel 1e-13 (`golden/stats/special-functions.txt`, `SpecialFunctionsTest`) | parity | #164 |
| commons-math3 3.4.1 `FDistribution`, `TDistribution`, `ChiSquaredDistribution`; Colt `Probability` | `stats.FDistribution`, `StudentTDistribution`, `ChiSquaredDistribution` | port | Colt tails and commons-math quantiles from the pinned standalone's jars | planned | #165 |

## Later phases

To be inventoried in Phase 1: the remaining TASSEL 5 plugins, each assigned
`port`, `bridge`, or `retire`, ordered by usage (rTASSEL calls, tutorials,
citations). Known high-priority ports for the teaching release: GLM, MLM
(EMMA/P3D), LD, NJ trees, Manhattan and QQ plots.

| TASSEL 5 source | TASSEL 6 target | Strategy | Golden test | Status | Issue |
|-----------------|-----------------|----------|-------------|--------|-------|
| `plugindef/GenerateRCode.java` | rTASSEL 2.0 generated wrappers (ADR-0007) | redesign | rTASSEL 1.0 API compatibility tests | planned | |

## Local stand-ins to replace

Temporary TASSEL 6 types that must be swapped for upstream equivalents.

| Type | Replace with | Blocked on |
|------|--------------|------------|
| `net.maizegenetics.tassel.position.Chromosome` | `biokotlin-core` contig type, keeping TASSEL 5 name normalization and ordering | `biokotlin-core` published (ADR-0008) |
| `net.maizegenetics.tassel.position.Position` | `biokotlin-core` `Position` (contig, position), with TASSEL 6 site fields (SNP name, strand, variants) kept alongside | `biokotlin-core` published (ADR-0008) |
