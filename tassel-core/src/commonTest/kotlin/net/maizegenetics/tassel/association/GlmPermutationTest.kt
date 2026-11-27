package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.stats.model.DesignMatrix
import net.maizegenetics.tassel.stats.model.LinearModelFit
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.sqrt
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue
import kotlin.time.Duration.Companion.minutes

/**
 * TASSEL 5's GLM permutation test (tassel-6#194). TASSEL 5's generator is unseeded, so there is no golden to match:
 * the null distribution is checked against an independent refit of explicitly permuted data, the seed against itself on
 * every worker count and (through pinned numbers the three targets must all reproduce) on every platform, and the
 * thresholds by simulation (departure D6, ADR-0029).
 */
class GlmPermutationTest {
    // --- the generator ---------------------------------------------------------------------------------------------

    @Test
    fun theGeneratorRepeatsSplitMix64sPublishedVectors() {
        // The first outputs of SplitMix64 from state 0 (Vigna's reference implementation).
        val random = PermutationRandom(0L)
        assertEquals(0xE220A8397B1DCDAFUL.toLong(), random.nextLong())
        assertEquals(0x6E789E6AA1B965F4UL.toLong(), random.nextLong())
        assertEquals(0x06C45D188009454FUL.toLong(), random.nextLong())
    }

    @Test
    fun boundedIntegersAreInRangeAndEvenlySpread() {
        val random = PermutationRandom(99L)
        val counts = IntArray(7)
        repeat(70_000) { counts[random.nextInt(7)]++ }
        // Each count is binomial(70000, 1/7): sd about 92. Six sd either side.
        counts.forEach { assertTrue(abs(it - 10_000) < 560, "count $it") }
        assertFailsWith<IllegalArgumentException> { random.nextInt(0) }
    }

    @Test
    fun aShuffleIsAPermutationAndEveryTraitHasItsOwnStream() {
        val a = IntArray(50) { it }.also { PermutationRandom.forTrait(5L, 0).shuffle(it) }
        val again = IntArray(50) { it }.also { PermutationRandom.forTrait(5L, 0).shuffle(it) }
        val other = IntArray(50) { it }.also { PermutationRandom.forTrait(5L, 1).shuffle(it) }
        assertContentEquals(a, again)
        assertNotEquals(a.toList(), other.toList())
        assertEquals((0 until 50).toList(), a.sorted())
        assertNotEquals((0 until 50).toList(), a.toList())
    }

    // --- perm_p and thresholds from a null distribution ------------------------------------------------------------

    private fun hand() = PermutationNull(
        PermutationOptions(5, 1L),
        listOf("t"),
        arrayOf(doubleArrayOf(0.01, 0.02, 0.02, 0.5, 0.9)),
    )

    @Test
    fun permPIsTheInsertionIndexPlusOneOverNperm() {
        val null5 = hand()
        assertEquals(1.0 / 5, null5.permutationP(0, 0.001))
        assertEquals(2.0 / 5, null5.permutationP(0, 0.015))
        // Equal to a minimum: only the minimums strictly below count, so the first of the tied 0.02s is the insertion point.
        assertEquals(2.0 / 5, null5.permutationP(0, 0.02))
        assertEquals(4.0 / 5, null5.permutationP(0, 0.3))
        // At or above the largest minimum the index is 4 or nperm, so (index + 1) / nperm is 1 or 6/5, capped at 1.
        assertEquals(1.0, null5.permutationP(0, 0.7))
        assertEquals(1.0, null5.permutationP(0, 0.95))
        assertTrue(null5.permutationP(0, Double.NaN).isNaN())
    }

    @Test
    fun aThresholdIsTheLargestPWhosePermPIsWithinAlpha() {
        val twenty = PermutationNull(
            PermutationOptions(20, 1L),
            listOf("t"),
            arrayOf(DoubleArray(20) { (it + 1) / 100.0 }),
        )
        // floor(0.05 * 20) = 1: the smallest minimum; floor(0.25 * 20) = 5.
        assertEquals(0.01, twenty.threshold(0, 0.05))
        assertEquals(0.05, twenty.threshold(0, 0.25))
        // 0.29 * 20 is 5.8, and floating point must not turn 0.15 * 20 = 3.0000000000000004 or 2.9999999999999996 around.
        assertEquals(0.05, twenty.threshold(0, 0.29))
        assertEquals(0.03, twenty.threshold(0, 0.15))
        // Fewer than one permutation per alpha: unresolved.
        assertTrue(twenty.threshold(0, 0.04).isNaN())
        for (alpha in listOf(0.05, 0.1, 0.25, 0.5)) {
            val t = twenty.threshold(0, alpha)
            for (p in listOf(t * 0.5, t, t + 0.0001, t * 3)) {
                assertEquals(twenty.permutationP(0, p) <= alpha + 1e-12, p <= t, "alpha $alpha p $p")
            }
        }
    }

    // --- the null against an independent refit of permuted data ----------------------------------------------------

    private fun permutedDataSet(
        input: DataSet,
        data: AssociationData,
        permutation: PermutationOptions,
        k: Int,
    ): DataSet {
        val phenotype = input.phenotypes().single()
        val replaced = HashMap<String, FloatArray>()
        data.traits.forEachIndexed { trait, name ->
            val rows = data.observations(trait)
            val design = DesignMatrix.build(data.baseSpec(trait), data.data.phenotype, rows, GLM_RANK_TOLERANCE)
            val fit = LinearModelFit.of(design)
            val random = PermutationRandom.forTrait(permutation.seed, trait)
            val arrangement = IntArray(rows.size) { it }
            repeat(k + 1) { random.shuffle(arrangement) }
            val values = FloatArray(phenotype.numberOfObservations) { Float.NaN }
            // TASSEL 5's fitted + shuffled residual. The trait mean is taken off (the base model has an intercept, so the
            // marker test cannot see it) to keep the values small, since the phenotype stores floats.
            val mean = fit.fitted.average()
            rows.forEachIndexed { i, row ->
                values[row] = (fit.fitted[i] - mean + fit.residuals[arrangement[i]]).toFloat()
            }
            replaced[name] = values
        }
        val attributes = phenotype.attributes.map { attribute ->
            replaced[attribute.name]?.let {
                NumericAttribute(attribute.name, AttributeType.data, it) as PhenotypeAttribute
            }
                ?: attribute
        }
        return DataSet.of(input.genotypes().single(), PhenotypeTable(attributes))
    }

    @Test
    fun theNullEqualsTheMinimumPOfRefittingExplicitlyPermutedData() = runTest(timeout = 5.minutes) {
        val input = randomGlmDesign(seed = 17, taxa = 45, markers = 30)
        val data = AssociationData.of(input)
        val permutation = PermutationOptions(permutations = 12, seed = 424242L)
        val fast = glmPermutationNull(data, GenotypeMarkerSource.of(data.data), GlmOptions(), permutation, workers = 2)

        val expected = Array(data.traits.size) { DoubleArray(permutation.permutations) { 1.0 } }
        for (k in 0 until permutation.permutations) {
            val permuted = permutedDataSet(input, data, permutation, k)
            val table = glmTable(
                AssociationData.of(permuted),
                GenotypeMarkerSource.of(AssociationData.of(permuted).data),
                kernel = GlmKernel.Reference,
            )
            for (record in table.records) {
                if (record.test != AssociationTestKind.Genotype || record.p.isNaN()) continue
                val trait = data.traits.indexOf(record.trait)
                if (record.p < expected[trait][k]) expected[trait][k] = record.p
            }
        }
        expected.forEach { it.sort() }
        data.traits.indices.forEach { trait ->
            // The refit reads the permuted values back from floats, which is good for a few 1e-7 of a p-value.
            assertClose(expected[trait], fast.minimumP[trait], 1e-5, "trait ${data.traits[trait]}")
        }
        // The test must be able to fail: the minimums are not all 1, and they differ across permutations.
        assertTrue(fast.minimumP.all { it.first() < 0.5 && it.first() != it.last() })
    }

    private fun assertClose(expected: DoubleArray, actual: DoubleArray, relative: Double, label: String) {
        assertEquals(expected.size, actual.size, label)
        for (i in expected.indices) {
            val bound = relative * max(abs(expected[i]), 1e-300)
            assertTrue(abs(expected[i] - actual[i]) <= bound, "$label [$i]: expected ${expected[i]}, got ${actual[i]}")
        }
    }

    // --- the seed --------------------------------------------------------------------------------------------------

    private suspend fun nullOf(
        seed: Long,
        workers: Int,
        permutations: Int = 25,
        markerBlock: Int = 16,
    ): PermutationNull {
        val data = AssociationData.of(randomGlmDesign(seed = 3, taxa = 50, markers = 40))
        return glmPermutationNull(
            data,
            GenotypeMarkerSource.of(data.data),
            GlmOptions(blockSize = markerBlock),
            PermutationOptions(permutations, seed),
            workers = workers,
        )
    }

    @Test
    fun theSameSeedGivesTheSameNullOnAnyNumberOfWorkers() = runTest(timeout = 5.minutes) {
        val one = nullOf(11L, workers = 1)
        for (workers in listOf(2, 5, 40)) {
            val many = nullOf(11L, workers = workers, markerBlock = 7)
            one.minimumP.indices.forEach {
                assertContentEquals(one.minimumP[it], many.minimumP[it], "workers $workers")
            }
        }
        val other = nullOf(12L, workers = 1)
        assertNotEquals(one.minimumP[0].toList(), other.minimumP[0].toList())
    }

    @Test
    fun theSeedReproducesPinnedNumbersOnEveryTarget() = runTest(timeout = 5.minutes) {
        // Numbers from the JVM. jvmTest, wasmJsNodeTest and linuxX64Test all run this: the same seed gives the same
        // output everywhere (the generator is integer arithmetic; the p-values agree to rounding).
        val result = nullOf(2024L, workers = 1, permutations = 8)
        assertClose(PINNED_Y1, result.minimumP[0], 1e-9, "y1")
    }

    // --- the table -------------------------------------------------------------------------------------------------

    @Test
    fun permPAndThresholdsLandInTheTableAndItsProvenance() = runTest(timeout = 5.minutes) {
        val data = AssociationData.of(randomGlmDesign(seed = 5, taxa = 60, markers = 25))
        val table = glmTable(
            data,
            GenotypeMarkerSource.of(data.data),
            permutation = PermutationOptions(100, 7L),
        )
        val summary = checkNotNull(table.provenance.permutations)
        assertEquals(100, summary.permutations)
        assertEquals(7L, table.provenance.seed)
        assertEquals(data.traits.toSet(), summary.thresholds.keys)
        assertContains(table.provenance.headerLines(), "##seed=7")
        assertContains(table.provenance.headerLines(), "##permutations=100")
        assertTrue(table.provenance.headerLines().any { it.startsWith("##permutation_threshold.y1.0.05=") })

        val omnibus = table.records.filter { it.test == AssociationTestKind.Genotype }
        assertTrue(omnibus.any { !it.p.isNaN() })
        for (record in omnibus) {
            if (record.p.isNaN()) {
                assertTrue(record.permutationP.isNaN())
                continue
            }
            assertTrue(record.permutationP in 0.01..1.0, "perm_p ${record.permutationP}")
            for ((alpha, threshold) in summary.thresholds.getValue(record.trait)) {
                assertEquals(record.permutationP <= alpha + 1e-12, record.p <= threshold, "${record.marker} at $alpha")
            }
        }
        // Only the omnibus record carries perm_p.
        assertTrue(table.records.filter { it.test == AssociationTestKind.Additive }.all { it.permutationP.isNaN() })
        // A smaller p never has a larger perm_p, trait by trait.
        omnibus.filter { !it.p.isNaN() }.groupBy { it.trait }.values.forEach { records ->
            val sorted = records.sortedBy { it.p }
            sorted.zipWithNext().forEach { (a, b) -> assertTrue(a.permutationP <= b.permutationP) }
        }
    }

    @Test
    fun withoutPermutationsThereIsNoPermPAndNoSeed() = runTest {
        val data = AssociationData.of(randomGlmDesign(seed = 5, taxa = 40, markers = 12))
        val table = glmTable(data, GenotypeMarkerSource.of(data.data))
        assertTrue(table.records.all { it.permutationP.isNaN() })
        assertEquals(null, table.provenance.seed)
        assertEquals(null, table.provenance.permutations)
    }

    @Test
    fun theOperationRunsPermutationsWithARecordedSeed() = runTest(timeout = 5.minutes) {
        val set = randomGlmDesign(seed = 5, taxa = 40, markers = 12)
        fun run(vararg values: Pair<String, Any?>) = glmOperation.bind(mapOf(*values))
        val a = glmOperation.run(set, run("permute" to true, "numberOfPermutations" to 20, "seed" to 9))
        val b = glmOperation.run(set, run("permute" to true, "numberOfPermutations" to 20, "seed" to 9))
        assertEquals(9L, a.provenance.seed)
        assertEquals(a.records.map { it.permutationP }.toList(), b.records.map { it.permutationP }.toList())
        // No seed given: one is drawn and recorded, so the run can be repeated with it.
        val drawn = glmOperation.run(set, run("permute" to true, "numberOfPermutations" to 20))
        val seed = checkNotNull(drawn.provenance.seed)
        val repeated = glmOperation.run(
            set,
            run("permute" to true, "numberOfPermutations" to 20, "seed" to seed.toInt()),
        )
        if (seed == seed.toInt().toLong()) {
            assertEquals(
                drawn.records.map {
                    it.permutationP
                }.toList(),
                repeated.records.map { it.permutationP }.toList(),
            )
        }
        val message = assertFailsWith<IllegalArgumentException> {
            glmOperation.run(set, run("permute" to true))
        }.message
        assertContains(message.orEmpty(), "nperm")
    }

    @Test
    fun replicatedTaxaAreRefused() = runTest {
        val base = randomGlmDesign(seed = 5, taxa = 20, markers = 12)
        val phenotype = base.phenotypes().single()
        val names = phenotype.taxaAttribute.names.toMutableList().also { it[1] = it[0] }
        val attributes = phenotype.attributes.map {
            if (it.type == AttributeType.taxa) PhenotypeTable.taxa(names) else it
        }
        val set = DataSet.of(base.genotypes().single(), PhenotypeTable(attributes))
        val message = assertFailsWith<IllegalArgumentException> {
            glmOperation.run(set, glmOperation.bind(mapOf("permute" to true, "numberOfPermutations" to 5, "seed" to 1)))
        }.message.orEmpty()
        assertContains(message, "replicated")
    }

    // --- the family-wise error rate --------------------------------------------------------------------------------

    /** A null data set: independent trait and markers (three classes each), one covariate, no effect anywhere. */
    private fun nullDataSet(random: Random, taxa: Int, markers: Int): DataSet {
        val names = (0 until taxa).map { "t$it" }
        val calls =
            Array(markers) {
                Array(taxa) { "ACG"[random.nextInt(3)].let { a -> "AC"[random.nextInt(2)].let { b -> "$a$b" } } }
            }
        val genotypes = InMemoryGenotypeTable.of(
            TaxaList.of(names),
            PositionList((0 until markers).map { Position(Chromosome.of(1), (it + 1) * 100, "m$it") }),
        ) { taxon, site -> Nucleotide.genotype(calls[site][taxon]) }
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(names),
                NumericAttribute("y", AttributeType.data, FloatArray(taxa) { random.nextGaussian().toFloat() }),
                NumericAttribute("c", AttributeType.covariate, FloatArray(taxa) { random.nextGaussian().toFloat() }),
            ),
        )
        return DataSet.of(genotypes, phenotype)
    }

    @Test
    fun theThresholdControlsTheFamilyWiseErrorRateWithinBinomialBounds() = runTest(timeout = 20.minutes) {
        val replicates = 1000
        val permutations = 100
        val random = Random(20261127)
        var at05 = 0
        var at01 = 0
        repeat(replicates) { replicate ->
            val data = AssociationData.of(nullDataSet(random, taxa = 36, markers = 12))
            val table = glmTable(
                data,
                GenotypeMarkerSource.of(data.data),
                permutation = PermutationOptions(permutations, seed = 1000L + replicate),
                workers = 1,
            )
            val thresholds = checkNotNull(table.provenance.permutations).thresholds.getValue("y")
            val minimum = table.records.filter { it.test == AssociationTestKind.Genotype && !it.p.isNaN() }
                .minOf { it.p }
            if (minimum <= thresholds.getValue(0.05)) at05++
            if (minimum <= thresholds.getValue(0.01)) at01++
        }
        // With the observed minimum exchangeable with the permuted ones, the chance that it falls within the k smallest of
        // nperm + 1 is k / (nperm + 1): 5/101 and 1/101. Four binomial standard deviations either side.
        assertWithin(at05, replicates, 5.0 / (permutations + 1), "alpha 0.05")
        assertWithin(at01, replicates, 1.0 / (permutations + 1), "alpha 0.01")
    }

    private fun assertWithin(count: Int, trials: Int, probability: Double, label: String) {
        val sd = sqrt(trials * probability * (1 - probability))
        val expected = trials * probability
        assertTrue(
            abs(count - expected) <= 4 * sd,
            "$label: $count of $trials rejected, expected $expected +- ${4 * sd}",
        )
    }

    private companion object {
        /** Pinned from the JVM; see `theSeedReproducesPinnedNumbersOnEveryTarget`. */
        val PINNED_Y1 = doubleArrayOf(
            0.0017515996616660836,
            0.002627675293300723,
            0.006509580251676624,
            0.016861893420541427,
            0.019801002258677147,
            0.022599030559264703,
            0.04341544151680434,
            0.0780248876174707,
        )
    }
}
