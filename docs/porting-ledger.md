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
| `dna/snp/io/BuilderFromVCF.java` (`HeaderPositions`, `ProcessVCFBlock`, the `GT` field only) | `io.Vcf.read`, byte-level `VcfHeader` and `VcfRecordParser` | port | loads `mdp_genotype_w_depth.vcf` as 281 x 3,093 and writes back to `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down, exact (`VcfTest`, `VcfEndToEndTest`) | parity | #14, #83 |
| `util/Utils.java` (`getBufferedReader`, the `.gz` wrap) | gzip `RawSource` (BGZF members one block at a time, plain gzip through a 32 KiB window) | port | `htslib/tutorial/mdp_genotype_w_depth.vcf.gz` and a `GZIPOutputStream` of the tutorial VCF versus `golden/vcf/mdp_genotype_w_depth.hmp.txt` from the `rs#` row down, exact (`GzipSourceTest`, `GzipStreamJvmTest`, `VcfEndToEndTest`); `htslib/variants/plain.vcf.gz` and `concatenated.vcf.gz` versus `source.vcf` (`GzipSourceTest`) | parity | #84, #85 |
| `dna/snp/io/BuilderFromVCF.java` depth, quality, and likelihood fields (`keepDepth`, `AD`, `DP`, `GQ`, `PL`) and `dna/snp/score/AlleleDepthBuilder.java` | none; `io.Vcf.read` reads past every `FORMAT` field after `GT` | retire | n/a; the HapMap golden carries no depth | retired | #14 |
| `taxa/TaxaListIOUtils.parseVCFHeadersIntoMap`, `taxa/Taxon` annotations (the VCF's `##SAMPLE=` headers) | none; read and dropped by `io.Vcf.read` | retire | TASSEL 6's VCF-to-HapMap output is missing the golden's 281 `##SAMPLE=` lines, so the golden is compared from the `rs#` row down (`VcfTest`, maintainer decision 2026-11-04) | retired | #14 |
| `dna/map/GeneralPosition` `WHICH_ALLELE.Reference` and `Alternate`, and its `INFO` annotations | none; `Position` has no reference, alternate, or annotation slot, so only `knownVariants` survives | retire | a VCF-backed site summary prints `N` for `Ref` and `Alt`, as #43 hardcodes (`GenotypeSummaryTest`) | retired | #14, #43 |
| `analysis/data/GenotypeSummaryPlugin.java` | `analysis.genotypeSummary`, `ops.summaryOperation` | port | all eight `summary/*` goldens, rounded and full precision (`GenotypeSummaryTest`) | parity | #15 |
| `analysis/filter/FilterSiteBuilderPlugin.java`, `dna/snp/FilterSite.java`, `analysis/filter/FilterBySites.java` | `filter.filterSites`, `filter.SiteFilter`, `filter-sites` operation | port | all six `golden/filter-sites/*.hmp.txt`, exact, through the kernel and (for `position-range-5-6`) through the spec (`FilterSitesTest`) | parity | #16, #45 |
| `analysis/filter/FilterTaxaBuilderPlugin.java`, `dna/snp/FilterTaxa.java`, `analysis/filter/FilterByTaxa.java` | `filter.filterTaxa`, `filter.TaxaFilter`, `filter-taxa` operation | port | all four `golden/filter-taxa/*.hmp.txt`, exact, through the kernel and (for `include-list` and `exclude-list`) through the spec (`FilterTaxaTest`) | parity | #16, #45 |
| `matrixalgebra/decomposition/SymmetricEigenvalueDecomposition.java` | `math.SymmetricEigenDecomposition` | port | R `prcomp` values from `PrinCompTest`, through the Gram matrix (`PrinCompOracleTest`); the PCA goldens (`PrincipalComponentsTest`) | parity | #17 |
| `analysis/distance/KinshipPlugin.java`, `analysis/distance/EndelmanDistanceMatrix.java` | `kinship` operation (Centered IBS). Default `distance.centeredIbs` is a blocked dosage product (ADR-0020); `distance.centeredIbsTassel5Order` keeps the bit-exact port | port | `golden/kinship/centered_ibs.txt` within `abs 8e-6`, `SumPk` within `abs 1e-4` (`CenteredIbsTest`); byte-exact text from the retained kernel | parity | #18, #62 |
| `taxa/distance/DistanceMatrix.java`, `DistanceMatrixBuilder.java`, `WriteDistanceMatrix.java` (`saveDelimitedDistanceMatrix`) | `distance.DistanceMatrix` and its `writeDelimited` | port | `golden/kinship/centered_ibs.txt` written from its values (`CenteredIbsTest`) | parity | #18 |
| `analysis/data/PrincipalComponentsPlugin.java`, `stats/PCA/PrinComp.java` | `pca` operation (Gram-matrix eigen), `analysis.principalComponents` | port | `golden/pca/eigenvectors*` and `golden/pca/pcs*` whole; `golden/pca/eigenvalues*` at their `rows-unordered` MANIFEST tolerance, because TASSEL 5's SVD does not sort (`PrincipalComponentsTest`, #47) | parity | #19 |
| `analysis/numericaltransform/NumericalGenotypePlugin.java` (`setAlternateMinorAllelesToMinor`), `dna/snp/GenotypeTableUtils.convertGenotypeToFloatProbability`, `analysis/numericaltransform/ImputationPlugin.java` and `ImputationByMean.java` (the `by_mean` path) | `analysis.ReferenceProbability` | port | the PCA goldens, which the taxon-mean imputation moves by 7.6% if it is wrong (`PrincipalComponentsTest`) | parity | #19 |
| `dna/snp/score/SiteScoreUtil.java`, `dna/snp/score/ReferenceProbabilityBuilder.java` | `analysis.SiteScore` | port | the 255-step round trip (`SiteScoreTest`); the PCA goldens, which need it to reach `rel 1e-9` | parity | #19 |
| `phenotype/Phenotype.java`, `phenotype/PhenotypeUtils.java` (`write`) | `table.PhenotypeTable` and its `writeDelimited` | redesign | `golden/pca/pcs.txt` layout and values (`PhenotypeTableTest`, `PrincipalComponentsTest`) | parity | #19 |
| `util/TableReport.java`, `util/SimpleTableReport.java` | `table.ResultTable`, `table.ResultTableBuilder` | redesign | via the summary and PCA table goldens (#15, #19) | in progress | #38 |
| `util/TableReportUtils.java` (`saveDelimitedTableReport`) | `table.ResultTable.writeDelimited` | port | `summary/alleles.txt` written from its values (`ResultTableTest`) | parity | #38 |
| `util/DoubleFormat.java` | `table.DoubleFormat` | port | every numeric cell of `summary/*.full.txt` and `pca/*.full.txt` against its rounded twin (`DoubleFormatGoldenTest`) | parity | #38 |
| `java.lang.Float.toString`, `Double.toString` (how TASSEL writes float matrix cells and quantitative annotations) | `table.floatToString`, `table.doubleToString` | port | the kinship golden's cells and `SumPk` (`CenteredIbsTest`); 400,000 values against the JVM (`ShortestDecimalJvmTest`) | parity | #18 |
| `plugindef/PluginParameter.java`, `plugindef/AbstractPlugin.java` | operation spec DSL (ADR-0002, ADR-0018) | redesign | n/a | in progress | #11, #45 |
| `pipeline/TasselPipeline.java` (`-h`, `-export`, `-exportType`, and the per-plugin flags it builds from `PluginParameter`) | `tassel-cli`: one subcommand per registered operation, generated from the specs (ADR-0002, ADR-0006) | redesign | `golden/kinship/centered_ibs.txt` written by `tassel kinship`, exact (`KinshipEndToEndTest`); `golden/vcf/mdp_genotype_w_depth.hmp.txt` written by `tassel filter-sites` from a `.vcf` input (`VcfEndToEndTest`) | parity (JVM, the CLI's only target) | #20, #14 |
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
