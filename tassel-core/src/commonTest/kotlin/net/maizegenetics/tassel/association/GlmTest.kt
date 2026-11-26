package net.maizegenetics.tassel.association

import kotlinx.coroutines.CancellationException
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.join.alignGenotypePhenotype
import net.maizegenetics.tassel.ops.InvalidParametersException
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * The GLM kernel, [AssociationData], and the `glm` spec on small hand-solved designs. The goldens (`GlmGoldenTest`) and
 * R (`GlmOracleTest`) say the numbers are right on real data; these say how the pieces behave: inputs, roles, skip
 * records, the order and invariance of results, progress and cancellation, and that no genotype matrix is held.
 */
class GlmTest {
    private fun genotypes(sites: List<String>, names: List<String>? = null): InMemoryGenotypeTable {
        val calls = sites.map { site -> site.split(" ").map { Nucleotide.genotype(it) } }
        val taxa = names ?: (0 until calls.first().size).map { "t$it" }
        return InMemoryGenotypeTable.of(
            TaxaList.of(taxa),
            PositionList(sites.indices.map { Position(Chromosome.of(1), (it + 1) * 100, "m$it") }),
        ) { taxon, site -> calls[site][taxon] }
    }

    private fun phenotype(
        taxa: List<String>,
        vararg attributes: net.maizegenetics.tassel.phenotype.PhenotypeAttribute,
    ) = PhenotypeTable(listOf(PhenotypeTable.taxa(taxa)) + attributes)

    private fun trait(name: String, vararg values: Float) = NumericAttribute(name, AttributeType.data, values)

    private suspend fun run(
        genotypes: GenotypeTable,
        phenotype: PhenotypeTable,
        options: GlmOptions = GlmOptions(),
        spec: AssociationSpec = AssociationSpec(),
    ): AssociationTable {
        val data = AssociationData.of(DataSet.of(genotypes, phenotype), spec)
        return glmTable(data, GenotypeMarkerSource.of(data.data), options)
    }

    private fun assertClose(expected: Double, actual: Double, relative: Double = 1e-12) {
        if (expected.isNaN() && actual.isNaN()) return
        assertTrue(abs(expected - actual) <= relative * abs(expected), "expected $expected but was $actual")
    }

    /** Six taxa, genotype classes AA AA AC AC CC CC, and a trait whose one-way ANOVA is easy by hand. */
    private val sixTaxa = listOf("t0", "t1", "t2", "t3", "t4", "t5")
    private val sixSite = "AA AA AC AC CC CC"
    private val sixTrait = trait("y", 1f, 2f, 4f, 6f, 9f, 11f)

    @Test
    fun oneWayAnovaByHand() = runTest {
        val table = run(genotypes(listOf(sixSite)), phenotype(sixTaxa, sixTrait))
        // Class means 1.5, 10, 5 (classes A, C, M) around 5.5: SS between 73, within 4.5, total 77.5.
        val main = table[0]
        assertEquals(AssociationTestKind.Genotype, main.test)
        assertEquals(2.0, main.numeratorDf)
        assertEquals(3.0, main.denominatorDf)
        assertClose(73.0 / 2.0 / 1.5, main.statistic)
        assertClose(73.0 / 77.5, main.rSquared)
        assertClose(36.5, main.markerMeanSquare)
        assertClose(1.5, main.errorMeanSquare)
        assertEquals(2.0, main.modelDf)
        assertEquals(6, main.observations)
        assertEquals(2, main.minorObs)
        // The last class (the heterozygote M) is the zero reference.
        assertEquals(listOf("A", "C", "M"), table.effects.map { it.allele }.toList())
        val estimates = table.effects.map { it.estimate }.toList()
        assertClose(-3.5, estimates[0])
        assertClose(5.0, estimates[1])
        assertEquals(0.0, estimates[2])
        // add: copies of A (6 of 12 alleles, tied with C and first by code): SS 17^2/4 against 5.25 on 4 df.
        val add = table[1]
        assertEquals(AssociationTestKind.Additive, add.test)
        assertClose(72.25 / ((77.5 - 72.25) / 4.0), add.statistic)
        val dom = table[2]
        assertEquals(AssociationTestKind.Dominance, dom.test)
        assertClose((73.0 - 72.25) / 1.5, dom.statistic)
    }

    @Test
    fun theAdditiveEncodingNamesItsEffectAllele() = runTest {
        val options = GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.additive))
        val table = run(genotypes(listOf(sixSite)), phenotype(sixTaxa, sixTrait), options)
        val record = table[0]
        assertEquals(1, table.size)
        assertEquals(AssociationTestKind.Additive, record.test)
        assertEquals(MarkerEncoding.additive, record.encoding)
        // Minor allele is C (tied with A, which sorts first and is the major allele); slope of y on C copies is 4.25.
        assertEquals("C", record.effectAllele)
        assertEquals("A", record.comparisonAllele)
        assertClose(4.25, record.effect)
        assertEquals(6, record.alleleCount)
        assertClose(0.5, record.alleleFrequency)
        assertEquals(6, record.minorAlleleCount)
        assertEquals(1.0, record.numeratorDf)
        assertEquals(4.0, record.denominatorDf)
        // F = t^2 and the interval is centered on the effect.
        val t = record.effect / record.standardError
        assertClose(t * t, record.statistic, 1e-10)
        assertClose(record.effect, (record.ciLower + record.ciUpper) / 2, 1e-12)
        assertEquals(0, table.effectCount, "the 1-df test has no per-class effects")
    }

    @Test
    fun additiveDominanceJointTestIsTheGenotypeTestAtABiallelicSite() = runTest {
        val genotype = run(genotypes(listOf(sixSite)), phenotype(sixTaxa, sixTrait))[0]
        val table = run(
            genotypes(listOf(sixSite)),
            phenotype(sixTaxa, sixTrait),
            GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.additiveDominance)),
        )
        assertEquals(
            listOf(AssociationTestKind.Genotype, AssociationTestKind.Additive, AssociationTestKind.Dominance),
            table.records.map { it.test }.toList(),
        )
        assertClose(genotype.statistic, table[0].statistic, 1e-10)
        assertEquals(2.0, table[0].numeratorDf)
        assertTrue(table.records.all { it.encoding == MarkerEncoding.additiveDominance })
        assertEquals("dom", table[2].effectAllele)
    }

    @Test
    fun filteredAndMonomorphicMarkersAreSkippedWithReasonsAndCounted() = runTest {
        val sites = listOf("AA AA AA AA AA AA", sixSite, "AA AA AC CG CC CC", "NN NN NN NN NN NN")
        val plain = run(genotypes(sites), phenotype(sixTaxa, sixTrait))
        // A monomorphic marker is reported, as TASSEL 5 does: no marker df and NaN statistics.
        val monomorphic = plain[0]
        assertEquals(AssociationStatus.Skipped, monomorphic.status)
        assertEquals(SkipReason.Monomorphic, monomorphic.skipReason)
        assertEquals(0.0, monomorphic.numeratorDf)
        assertEquals(5.0, monomorphic.denominatorDf)
        assertTrue(monomorphic.statistic.isNaN() && monomorphic.p.isNaN())
        assertEquals(1L, plain.summary.skipped.getValue(SkipReason.Monomorphic))
        assertEquals(1L, plain.summary.skipped.getValue(SkipReason.NoObservations))
        assertEquals(6L, plain.summary.tested, "omnibus, add, and dom for each of the two markers tested")

        val filtered = run(
            genotypes(sites),
            phenotype(sixTaxa, sixTrait),
            GlmOptions(marker = MarkerOptions(biallelicOnly = true)),
        )
        // TASSEL 5 drops a marker biallelicOnly rejects, a monomorphic one included.
        assertEquals(
            listOf(SkipReason.NotBiallelic, null, null, null, SkipReason.NotBiallelic, SkipReason.NoObservations),
            filtered.records.map { it.skipReason }.toList(),
        )
        assertEquals(3L, filtered.summary.tested)
        assertEquals(3L, filtered.summary.skippedTotal)
    }

    @Test
    fun resultsDoNotDependOnObservationOrTraitOrder() = runTest {
        val sites = listOf(sixSite, "AA AC AC CC CC AA", "CC AC AA AC CC AA")
        val y2 = trait("z", 5f, 3f, 8f, 1f, 4f, 2f)
        val forward = run(genotypes(sites), phenotype(sixTaxa, sixTrait, y2))
        val order = listOf(4, 2, 5, 0, 3, 1)
        fun <T> List<T>.permute() = order.map { this[it] }
        val shuffled = run(
            genotypes(sites),
            phenotype(
                sixTaxa.permute(),
                trait("z", *floatArrayOf(5f, 3f, 8f, 1f, 4f, 2f).toList().permute().toFloatArray()),
                trait("y", *floatArrayOf(1f, 2f, 4f, 6f, 9f, 11f).toList().permute().toFloatArray()),
            ),
        )
        val key = { r: AssociationRecord -> Triple(r.trait, r.markerIndex, r.test) }
        val expected = forward.records.associateBy(key)
        assertEquals(expected.size, shuffled.size)
        for (record in shuffled.records) {
            val other = expected.getValue(key(record))
            assertClose(other.statistic, record.statistic, 1e-10)
            assertClose(other.rSquared, record.rSquared, 1e-10)
            assertEquals(other.denominatorDf, record.denominatorDf)
        }
        // Traits are listed in the analysis order, which is the phenotype's.
        assertEquals(listOf("z", "y"), shuffled.traits)
        assertEquals(listOf("y", "z"), forward.traits)
    }

    // --- inputs ----------------------------------------------------------------------------------------------

    private val roleTable = phenotype(
        listOf("t0", "t1", "t2", "t3", "t4", "t5"),
        trait("y", 1f, 2f, Float.NaN, 6f, 9f, 11f),
        trait("w", 3f, 1f, 2f, 5f, 2f, 7f),
        NumericAttribute("q", AttributeType.covariate, floatArrayOf(1f, 2f, 3f, Float.NaN, 5f, 6f)),
        FactorAttribute.of("env", listOf("a", "b", "a", "b", "a", null)),
    )

    @Test
    fun rolesFollowTassel5AndMissingValuesGiveEachTraitItsOwnObservations() {
        val data = AssociationData.of(alignGenotypePhenotype(genotypes(listOf(sixSite)), roleTable))
        assertEquals(listOf("y", "w"), data.traits)
        assertEquals(listOf("q"), data.covariates)
        assertEquals(listOf("env"), data.factors)
        // y: row 2 (trait), 3 (covariate), 5 (factor) are out; w: rows 3 and 5.
        assertEquals(listOf(0, 1, 4), data.observations(0).toList())
        assertEquals(listOf(0, 1, 2, 4), data.observations(1).toList())
        assertEquals(mapOf("trait" to 1, "covariate" to 1, "factor" to 1), data.excludedObservations[0])
        assertEquals(mapOf("covariate" to 1, "factor" to 1), data.excludedObservations[1])
        assertEquals(
            listOf("(Intercept)", "env", "q"),
            data.baseSpec(0).terms.map { it.label },
            "the mean, then factors, then covariates, as TASSEL 5 orders them",
        )
    }

    @Test
    fun namedSelectionsOverrideTheDefaultsAndAreChecked() {
        val aligned = alignGenotypePhenotype(genotypes(listOf(sixSite)), roleTable)
        val chosen = AssociationData.of(
            aligned,
            AssociationSpec(traits = listOf("w"), covariates = emptyList(), factors = emptyList()),
        )
        assertEquals(listOf("w"), chosen.traits)
        assertEquals(listOf(0, 1, 2, 3, 4, 5), chosen.observations(0).toList())
        fun message(spec: AssociationSpec) =
            assertFailsWith<IllegalArgumentException> { AssociationData.of(aligned, spec) }.message.orEmpty()
        assertContains(message(AssociationSpec(traits = listOf("nope"))), "'nope' is not an attribute")
        assertContains(message(AssociationSpec(traits = listOf("y", "y"))), "listed twice")
        assertContains(message(AssociationSpec(covariates = listOf("y"))), "both a trait and a covariate")
        assertContains(message(AssociationSpec(factors = listOf("q"))), "must be a factor")
        assertContains(message(AssociationSpec(covariates = listOf("env"))), "must be numeric")
        assertContains(message(AssociationSpec(traits = emptyList())), "no trait")
    }

    @Test
    fun aDataSetIsOneGenotypePhenotypeOrOneGenotypeTableAndPhenotypes() = runTest {
        val table = genotypes(listOf(sixSite))
        val pairing = alignGenotypePhenotype(table, roleTable)
        assertEquals(pairing, AssociationData.of(DataSet.of(pairing)).data)
        val two = AssociationData.of(
            DataSet.of(
                table,
                phenotype(sixTaxa, trait("y", 1f, 2f, 3f, 4f, 5f, 6f)),
                phenotype(sixTaxa, NumericAttribute("q", AttributeType.covariate, FloatArray(6) { it.toFloat() })),
            ),
        )
        assertEquals(listOf("y"), two.traits)
        assertEquals(listOf("q"), two.covariates)
        suspend fun message(set: DataSet) =
            assertFailsWith<IllegalArgumentException> { AssociationData.of(set) }.message.orEmpty()
        assertContains(message(DataSet.of(table)), "0 pairings, 1 genotype tables, 0 phenotype tables")
        assertContains(message(DataSet.of(roleTable)), "0 genotype tables, 1 phenotype tables")
        assertContains(message(DataSet.of(table, roleTable, "kinship")), "1 other items")
    }

    // --- progress, cancellation, and memory --------------------------------------------------------------------

    private fun largeTable(): Pair<InMemoryGenotypeTable, PhenotypeTable> {
        var state = 12345L
        fun next(): Int {
            state = state * 6364136223846793005L + 1442695040888963407L
            return ((state ushr 33) % 3).toInt()
        }
        val taxa = 60
        val calls = Array(40) { Array(taxa) { arrayOf("AA", "AC", "CC")[next()] } }
        val table = genotypes(calls.map { it.joinToString(" ") })
        val names = (0 until taxa).map { "t$it" }
        val values = FloatArray(taxa) { (it * 7 % 13).toFloat() }
        val second = FloatArray(taxa) { (it * 5 % 11).toFloat() }
        return table to phenotype(names, trait("y1", *values), trait("y2", *second))
    }

    @Test
    fun reportsProgressAfterEverySiteBlock() = runTest {
        val (table, pheno) = largeTable()
        val data = AssociationData.of(alignGenotypePhenotype(table, pheno))
        val reports = ArrayList<Progress>()
        glmTable(
            data,
            GenotypeMarkerSource.of(data.data),
            GlmOptions(blockSize = 16),
            progress = ProgressReporter { reports += it },
        )
        assertEquals(listOf(0L, 16L, 32L, 40L), reports.map { it.completed })
        assertTrue(reports.all { it.total == 40L })
    }

    @Test
    fun cancellationStopsWithinOneSiteBlock() = runTest {
        val (table, pheno) = largeTable()
        val data = AssociationData.of(alignGenotypePhenotype(table, pheno))
        val chunks = ArrayList<AssociationChunk>()
        var reports = 0
        assertFailsWith<CancellationException> {
            glm(
                data,
                GenotypeMarkerSource.of(data.data),
                GlmOptions(blockSize = 10),
                sink = object : AssociationSink {
                    override suspend fun accept(chunk: AssociationChunk) {
                        chunks += chunk
                    }
                },
                progress = ProgressReporter {
                    // Report 0 is the start; report 2 follows the second block, where the run is cancelled.
                    if (++reports == 3) throw CancellationException("stop")
                },
            )
        }
        assertEquals(2, chunks.size, "no block is analyzed after the report that cancels")
        assertEquals(listOf(0, 1), chunks.map { it.sequence })
    }

    @Test
    fun noFullGenotypeMatrixIsMaterialized() = runTest {
        val (table, pheno) = largeTable()
        val counting = Counting(table)
        val data = AssociationData.of(alignGenotypePhenotype(counting, pheno))
        assertEquals(0, counting.cellsRead, "aligning reads no genotype")
        val result = glmTable(data, GenotypeMarkerSource.of(data.data), GlmOptions(blockSize = 7))
        assertTrue(result.size > 0)
        // Each genotype is read once for the whole run, however many traits share it, and always by site block.
        assertEquals(60 * 40, counting.cellsRead)
        assertEquals(0, counting.wholeSiteOrTaxonReads, "no whole-site or whole-taxon vector is ever requested")
    }

    /** A table that counts how often genotypes are read, and does not override `siteBlocks`. */
    private class Counting(private val inner: GenotypeTable) : GenotypeTable {
        var cellsRead = 0
        var wholeSiteOrTaxonReads = 0
        override val taxa get() = inner.taxa
        override val positions get() = inner.positions

        override fun genotype(taxon: Int, site: Int): Byte {
            cellsRead++
            return inner.genotype(taxon, site)
        }

        override fun siteGenotypes(site: Int): ByteArray {
            wholeSiteOrTaxonReads++
            return super.siteGenotypes(site)
        }

        override fun taxonGenotypes(taxon: Int): ByteArray {
            wholeSiteOrTaxonReads++
            return super.taxonGenotypes(taxon)
        }
    }

    // --- the spec ------------------------------------------------------------------------------------------------

    @Test
    fun theSpecKeepsTassel5sParameterNamesAndDefaults() {
        assertEquals("glm", glmOperation.id)
        assertEquals("FixedEffectLMPlugin", glmOperation.legacyPlugin)
        assertEquals(DataSet::class, glmOperation.inputType)
        val legacy = mapOf(
            "genotypeComponent" to "genotypeComponent",
            "minClassSize" to "minClassSize",
            "biallelicOnly" to "biallelicOnly",
            "appendAddDom" to "appendAddDom",
            "maxP" to "maxP",
            "permute" to "permute",
            "numberOfPermutations" to "nperm",
            "phenotypeOnly" to "phenoOnly",
            "saveToFile" to "saveToFile",
            "siteFile" to "siteFile",
            "alleleFile" to "alleleFile",
            "bluesFile" to "bluesFile",
            "anovaFile" to "anovaFile",
        )
        legacy.forEach { (name, tassel5) ->
            assertEquals(name, glmOperation.parameterByLegacyName(tassel5)?.name, tassel5)
        }
        assertNull(glmOperation.parameterByLegacyName("siteStatsOut"), "dead in TASSEL 5, not ported")
        assertEquals(MarkerEncoding.genotype, glmOperation.parameter("encoding")?.default)
        assertEquals(1.0, glmOperation.parameter("maxP")?.default)
        assertEquals(0, glmOperation.parameter("minClassSize")?.default)
        assertEquals(false, glmOperation.parameter("biallelicOnly")?.default)
        assertEquals(MissingGenotypes.drop, glmOperation.parameter("missingGenotypes")?.default)
    }

    @Test
    fun runningTheSpecIsTheSameAsCallingTheKernel() = runTest {
        val set = DataSet.of(genotypes(listOf(sixSite, "AA AC AC CC CC AA")), phenotype(sixTaxa, sixTrait))
        val viaSpec = glmOperation.run(set, glmOperation.bind(mapOf("maxP" to 1.0, "appendAddDom" to true)))
        val data = AssociationData.of(set)
        val direct = glmTable(
            data,
            GenotypeMarkerSource.of(data.data),
            GlmOptions(appendAddDom = true),
        )
        assertEquals(direct.size, viaSpec.size)
        assertEquals(direct.records.map { it.statistic }.toList(), viaSpec.records.map { it.statistic }.toList())
        assertEquals("glm", viaSpec.provenance.operation)
        assertEquals("true", viaSpec.provenance.parameters["appendAddDom"])
    }

    @Test
    fun theSpecAppliesTheFiltersAndSelections() = runTest {
        val set = DataSet.of(genotypes(listOf(sixSite, "AA AC AC CC CC AA")), roleTable)
        val table = glmOperation.run(
            set,
            glmOperation.bind(
                mapOf(
                    "traits" to listOf("w"),
                    "covariates" to emptyList<String>(),
                    "factors" to emptyList<String>(),
                    "minMac" to 7,
                ),
            ),
        )
        // minMac 7 skips both markers (6 minor copies each) for the only trait.
        assertEquals(listOf("w"), table.traits)
        assertTrue(table.records.all { it.skipReason == SkipReason.MinorAlleleCount })
    }

    @Test
    fun theSpecRefusesWhatLaterTasksPort() = runTest {
        val set = DataSet.of(genotypes(listOf(sixSite)), phenotype(sixTaxa, sixTrait))
        suspend fun message(vararg values: Pair<String, Any?>) = assertFailsWith<IllegalArgumentException> {
            glmOperation.run(set, glmOperation.bind(mapOf(*values)))
        }.message.orEmpty()
        assertContains(message("permute" to true), "#194")
        assertContains(message("phenotypeOnly" to true), "#193")
        assertContains(message("siteFile" to "stats.txt"), "#209")
        assertContains(message("missingGenotypes" to MissingGenotypes.imputeMean), "not ported")
        assertContains(message("encoding" to MarkerEncoding.dosage), "ReferenceProbability")
        assertContains(
            message(
                "genotypeComponent" to GlmGenotypeComponent.ReferenceProbability,
                "encoding" to MarkerEncoding.additive,
            ),
            "tests a dosage",
        )
        assertFailsWith<InvalidParametersException> { glmOperation.bind(mapOf("maxP" to 2.0)) }
    }
}
