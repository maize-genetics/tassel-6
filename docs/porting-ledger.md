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
| `dna/snp/NucleotideAlignmentConstants.java` | `tassel-core` genotype encoding | port | encoding round-trip | planned | |
| `util/SuperByteMatrix.java`, `dna/snp/genotypecall/ByteGenotypeCallTable.java` | `GenotypeTable` in-memory store with site blocks (ADR-0012) | redesign | via all I/O goldens | planned | |
| `dna/snp/genotypecall/AlleleFreqCache.java` | allele statistics | port | via summary goldens | planned | |
| `taxa/TaxaList*.java`, `dna/map/PositionList*.java` | taxa and position lists (with `biokotlin-core` types, ADR-0008) | redesign | via I/O goldens | planned | |
| `dna/snp/io/BuilderFromHapMap.java` | HapMap reader | port | load `mdp_genotype.hmp.txt` | planned | |
| `dna/snp/ExportUtils.java` (HapMap writer) | HapMap writer | port | write/read round-trip vs TASSEL 5 export | planned | |
| `dna/snp/io/BuilderFromVCF.java` | minimal VCF reader (GT only) | port | load `mdp_genotype_w_depth.vcf` | planned | |
| `analysis/data/GenotypeSummaryPlugin.java` | `summary` operation | port | overall, site, and taxa summary tables | planned | |
| `analysis/filter/FilterSiteBuilderPlugin.java`, `dna/snp/FilterSite.java` | `filter-sites` operation | port | MAF, count, and het filters | planned | |
| `analysis/filter/FilterTaxaBuilderPlugin.java`, `dna/snp/FilterTaxa.java` | `filter-taxa` operation | port | coverage and het filters | planned | |
| `analysis/distance/KinshipPlugin.java`, `analysis/distance/EndelmanDistanceMatrix.java` | `kinship` operation (Centered IBS) | port | `Centered_IBS_Matrix.txt`, tol `8e-6` | planned | |
| `taxa/distance/DistanceMatrix.java`, `WriteDistanceMatrix.java` | kinship matrix type and writer | port | text output matches | planned | |
| `analysis/data/PrincipalComponentsPlugin.java`, `stats/PCA/PrinComp.java` | `pca` operation (Gram-matrix eigen) | port | eigenvalues, vectors up to sign | planned | |
| `plugindef/PluginParameter.java`, `plugindef/AbstractPlugin.java` | operation spec DSL (ADR-0002) | redesign | n/a | planned | |
| `pipeline/TasselPipeline.java` | `tassel translate` (ADR-0006) | redesign | translated commands give golden outputs | planned (stretch) | |

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
| _none yet_ | | |
