package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import kotlin.math.abs
import kotlin.math.max
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue
import kotlin.test.fail
import kotlin.time.Duration.Companion.minutes

/**
 * The fast GLM kernel (`glm`, tassel-6#192) against the reference kernel (`glmReference`, #191) it replaces by default:
 * the golden inputs and random designs with missing genotypes, through both ways of handling a marker's dropped rows,
 * every encoding, and any number of workers.
 *
 * The two kernels share the record assembly, so these tests compare numerics only. Counts, degrees of freedom, labels,
 * order, and skip reasons must be identical. Statistics, mean squares, estimates, and standard errors agree to
 * [RELATIVE], ADR-0029's GLM tolerance (decision Q7), and a p-value below `1e-300` is compared as `log10 p` to
 * [LOG10_P_ABSOLUTE].
 *
 * A number that is zero in exact arithmetic is rounding noise in both kernels, where a relative tolerance means
 * nothing. As for departure D4 in `testdata/golden/MANIFEST`, a sum of squares may then differ by a multiple of eps of
 * the error sum of squares: [cancellationEpsMultiple] (twice D4's multiple, read from that MANIFEST, since both sides carry it), twice that for
 * `dom_F`, which TASSEL 5 defines as a difference of two sums of squares. A p-value is a function of F that both
 * kernels share, so it is compared only where F agrees to [RELATIVE]; where F agrees only within that bound, its p-value
 * is a function of rounding noise (and a slightly negative F has none). An estimate that is zero in exact arithmetic is
 * held to [NOISE_ABSOLUTE] of the trait's scale. The goldens themselves are checked
 * through the fast kernel in `GlmGoldenTest`.
 */
class GlmFastTest {
    private companion object {
        /** ADR-0029's starting GLM tolerance for statistics, decision Q7. */
        const val RELATIVE = 1e-8

        /** ADR-0029's tolerance for `log10 p` past underflow. */
        const val LOG10_P_ABSOLUTE = 1e-6

        /** For a value that is zero in exact arithmetic, a bound relative to the response's own scale. */
        const val NOISE_ABSOLUTE = 1e-11

        /**
         * Twice D4's bound (`bound.cancellation_eps_multiple` in the `[glm-departures]` section of
         * `testdata/golden/MANIFEST`), since the fast kernel and the reference both round.
         */
        val cancellationEpsMultiple = 2 * GlmDepartures.cancellationEpsMultiple

        const val EPSILON = 2.220446049250313E-16

        val HEAVY = 10.minutes

        val chr9and10 = SuspendLazy { hapmap("tutorial/mdp_genotype_chr9_10.hmp.txt") }
        val synthetic = SuspendLazy { hapmap("association/glm_synthetic.hmp.txt") }
        val traits = SuspendLazy { phenotype("tutorial/mdp_traits.txt") }
        val structure = SuspendLazy { phenotype("tutorial/mdp_population_structure.txt") }
        val syntheticTraits = SuspendLazy { phenotype("association/glm_synthetic_traits.txt") }
        val syntheticFactor = SuspendLazy { phenotype("association/glm_synthetic_factor.txt") }
        val syntheticCollinear = SuspendLazy { phenotype("association/glm_synthetic_collinear.txt") }

        suspend fun hapmap(path: String): GenotypeTable =
            HapMap.read(Buffer().also { it.write(TestData.readBytes(path)) })

        suspend fun phenotype(path: String): PhenotypeTable =
            PhenotypeFiles.read(Buffer().also { it.write(TestData.readBytes(path)) }).table
    }

    // --- comparison -----------------------------------------------------------------------------------------------

    private suspend fun tables(
        input: DataSet,
        spec: AssociationSpec = AssociationSpec(),
        options: GlmOptions = GlmOptions(),
        tuning: GlmFastTuning = GlmFastTuning(),
        workers: Int = 1,
        stats: GlmFastStats? = null,
    ): Pair<AssociationTable, AssociationTable> {
        val data = AssociationData.of(input, spec)
        suspend fun source() = if (options.marker.encoding == MarkerEncoding.dosage) {
            GenotypeMarkerSource.withDosage(data.data)
        } else {
            GenotypeMarkerSource.of(data.data)
        }
        val reference = glmTable(data, source(), options, kernel = GlmKernel.Reference)
        val collector = AssociationCollector(options.maxP, MaxPRule.Glm, data.traits)
        glm(data, source(), options, collector, workers = workers, tuning = tuning, stats = stats)
        val fast = collector.finish(reference.provenance)
        return reference to fast
    }

    private suspend fun assertFastMatchesReference(
        input: DataSet,
        spec: AssociationSpec = AssociationSpec(),
        options: GlmOptions = GlmOptions(),
        tuning: GlmFastTuning = GlmFastTuning(),
        stats: GlmFastStats? = null,
    ) {
        val (reference, fast) = tables(input, spec, options, tuning, stats = stats)
        assertSameResults(reference, fast)
    }

    /** Scale of each trait's numbers, for values that are zero in exact arithmetic: its error mean square's root. */
    private fun scales(table: AssociationTable): Map<String, Double> = table.records
        .filter { it.test != AssociationTestKind.Dominance && it.errorMeanSquare.isFinite() }
        .groupBy { it.trait }
        .mapValues { (_, records) -> records.maxOf { kotlin.math.sqrt(it.errorMeanSquare) }.coerceAtLeast(1.0) }

    private fun assertSameResults(reference: AssociationTable, fast: AssociationTable) {
        assertEquals(reference.size, fast.size, "record count")
        assertEquals(reference.effectCount, fast.effectCount, "effect count")
        assertEquals(reference.traits, fast.traits)
        assertEquals(reference.summary, fast.summary)
        val scale = scales(reference)
        val problems = ArrayList<String>()
        for (i in 0 until reference.size) {
            val r = reference[i]
            val f = fast[i]
            val where = "${r.trait} ${r.marker} ${r.test}"
            val traitScale = scale[r.trait] ?: 1.0
            fun same(name: String, a: Any?, b: Any?) {
                if (a != b) problems += "$where $name: $a vs $b"
            }
            fun close(name: String, a: Double, b: Double, noise: Double = 0.0) {
                if (!near(a, b, noise)) problems += "$where $name: $a vs $b"
            }
            same(
                "identity",
                listOf(r.trait, r.markerIndex, r.marker, r.chromosome, r.position, r.test, r.encoding),
                listOf(f.trait, f.markerIndex, f.marker, f.chromosome, f.position, f.test, f.encoding),
            )
            same("alleles", r.effectAllele to r.comparisonAllele, f.effectAllele to f.comparisonAllele)
            same(
                "counts",
                listOf(r.observations, r.taxa, r.alleleCount, r.minorAlleleCount, r.minorObs),
                listOf(f.observations, f.taxa, f.alleleCount, f.minorAlleleCount, f.minorObs),
            )
            same(
                "df",
                listOf(r.numeratorDf, r.denominatorDf, r.modelDf),
                listOf(f.numeratorDf, f.denominatorDf, f.modelDf),
            )
            same("status", Triple(r.status, r.skipReason, r.message), Triple(f.status, f.skipReason, f.message))
            close("alleleFrequency", r.alleleFrequency, f.alleleFrequency)
            close("effect", r.effect, f.effect, NOISE_ABSOLUTE * traitScale)
            close("standardError", r.standardError, f.standardError, NOISE_ABSOLUTE * traitScale)
            close("ciLower", r.ciLower, f.ciLower, NOISE_ABSOLUTE * traitScale)
            close("ciUpper", r.ciUpper, f.ciUpper, NOISE_ABSOLUTE * traitScale)
            // D4's noise: a multiple of eps of the error sum of squares, in each field's units.
            val bound = cancellationEpsMultiple * EPSILON * if (r.test == AssociationTestKind.Dominance) 2 else 1
            val errorSs = r.errorMeanSquare * r.denominatorDf
            val noiseF = bound * r.denominatorDf
            close("statistic", r.statistic, f.statistic, noiseF)
            close("rSquared", r.rSquared, f.rSquared, bound)
            close("markerMeanSquare", r.markerMeanSquare, f.markerMeanSquare, bound * errorSs / r.numeratorDf)
            close("modelMeanSquare", r.modelMeanSquare, f.modelMeanSquare, bound * errorSs / r.modelDf)
            close(
                "errorMeanSquare",
                r.errorMeanSquare,
                f.errorMeanSquare,
                bound * r.modelMeanSquare * r.modelDf / r.denominatorDf,
            )
            if (!near(r.statistic, f.statistic, 0.0)) {
                // F agrees within the noise bound only (checked above); its p-values are those of rounding noise.
            } else if (r.p < 1e-300 || f.p < 1e-300) {
                if (!(abs(r.log10P - f.log10P) <= LOG10_P_ABSOLUTE)) {
                    problems +=
                        "$where log10P: ${r.log10P} vs ${f.log10P}"
                }
            } else {
                close("p", r.p, f.p)
                close("log10P", r.log10P, f.log10P, NOISE_ABSOLUTE)
            }
        }
        for (i in 0 until reference.effectCount) {
            val r = reference.effect(i)
            val f = fast.effect(i)
            val where = "${r.trait} ${r.marker} effect ${r.allele}"
            val traitScale = scale[r.trait] ?: 1.0
            if (listOf(r.trait, r.marker, r.allele, r.observations) !=
                listOf(f.trait, f.marker, f.allele, f.observations)
            ) {
                problems += "$where: identity $r vs $f"
            }
            if (!near(r.estimate, f.estimate, NOISE_ABSOLUTE * traitScale)) {
                problems +=
                    "$where estimate: ${r.estimate} vs ${f.estimate}"
            }
            if (!near(r.standardError, f.standardError, NOISE_ABSOLUTE * traitScale)) {
                problems += "$where standardError: ${r.standardError} vs ${f.standardError}"
            }
        }
        if (problems.isNotEmpty()) {
            fail("${problems.size} differences from the reference:\n" + problems.take(20).joinToString("\n"))
        }
    }

    private fun near(a: Double, b: Double, absolute: Double): Boolean {
        if (a.isNaN() || b.isNaN()) return a.isNaN() && b.isNaN()
        if (a.isInfinite() || b.isInfinite()) return a == b
        return abs(a - b) <= RELATIVE * max(abs(a), abs(b)) || abs(a - b) <= absolute
    }

    private fun AssociationTable.text(): List<String> =
        records.map { it.toString() }.toList() + effects.map { it.toString() }

    // --- golden inputs ----------------------------------------------------------------------------------------------

    private val tutorialSpec = AssociationSpec(covariates = listOf("Q1", "Q2"))

    @Test
    fun matchesTheReferenceOnTheTutorialChromosomesWithPopulationStructure() = runTest(timeout = HEAVY) {
        val input = DataSet.of(chr9and10(), traits(), structure())
        for (dropped in DroppedRows.entries) {
            val stats = GlmFastStats()
            assertFastMatchesReference(
                input,
                tutorialSpec,
                tuning = GlmFastTuning(droppedRows = dropped),
                stats = stats,
            )
            assertTrue(stats.fullSet > 0, "some markers have a call for every observation")
            assertTrue(stats.indicators + stats.patternHits + stats.patternMisses > 0, "and some drop observations")
        }
    }

    @Test
    fun matchesTheReferenceWithEveryGoldenOption() = runTest(timeout = HEAVY) {
        val input = DataSet.of(chr9and10(), traits(), structure())
        for (options in listOf(
            GlmOptions(marker = MarkerOptions(minClassSize = 5)),
            GlmOptions(marker = MarkerOptions(biallelicOnly = true)),
            GlmOptions(maxP = 0.01),
            GlmOptions(appendAddDom = true),
        )) {
            assertFastMatchesReference(input, tutorialSpec, options)
        }
    }

    @Test
    fun matchesTheReferenceWithoutCovariatesAndOnReferenceProbabilities() = runTest(timeout = HEAVY) {
        val input = DataSet.of(chr9and10(), traits())
        assertFastMatchesReference(input)
        assertFastMatchesReference(
            input,
            options = GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.dosage)),
        )
    }

    @Test
    fun matchesTheReferenceOnTheSyntheticGoldens() = runTest {
        for (phenotype in listOf(syntheticTraits(), syntheticFactor(), syntheticCollinear())) {
            for (options in listOf(
                GlmOptions(),
                GlmOptions(marker = MarkerOptions(minClassSize = 3)),
                GlmOptions(marker = MarkerOptions(biallelicOnly = true)),
                GlmOptions(appendAddDom = true),
            )) {
                assertFastMatchesReference(DataSet.of(synthetic(), phenotype), options = options)
            }
        }
    }

    // --- random designs -------------------------------------------------------------------------------------------

    private val encodings = listOf(
        GlmOptions(),
        GlmOptions(appendAddDom = true),
        GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.additive)),
        GlmOptions(
            marker = MarkerOptions(encoding = MarkerEncoding.additive, multiallelic = MultiallelicAdditive.Expand),
        ),
        GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.additiveDominance)),
        GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.dosage)),
    )

    @Test
    fun matchesTheReferenceOnRandomDesignsWithMissingGenotypes() = runTest(timeout = HEAVY) {
        for (seed in 1..3) {
            val input = randomGlmDesign(seed)
            for (options in encodings) {
                for (dropped in DroppedRows.entries) {
                    val stats = GlmFastStats()
                    assertFastMatchesReference(input, options = options, tuning = GlmFastTuning(dropped), stats = stats)
                    when (dropped) {
                        DroppedRows.Indicators -> {
                            assertTrue(stats.indicators > 0, "dropped rows in closed form: $seed $options")
                            // Markers missing every member of the rare factor level lose a base column: pattern base.
                            assertTrue(stats.patternMisses > 0, "the rare level's markers take a pattern base")
                        }
                        DroppedRows.PatternBase -> {
                            assertEquals(0L, stats.indicators)
                            assertTrue(
                                stats.patternHits > 0,
                                "the shared pattern is found in the cache: $seed $options",
                            )
                        }
                    }
                }
            }
        }
    }

    @Test
    fun matchesTheReferenceWithAnInterceptOnlyBase() = runTest {
        val spec = AssociationSpec(covariates = emptyList(), factors = emptyList())
        for (options in encodings) {
            for (dropped in DroppedRows.entries) {
                assertFastMatchesReference(randomGlmDesign(4), spec, options, GlmFastTuning(dropped))
            }
        }
    }

    @Test
    fun matchesTheReferenceWhenMostCallsAreMissing() = runTest {
        // Markers that keep only a handful of observations: more marker columns than residual rows, and no error df.
        for (dropped in DroppedRows.entries) {
            assertFastMatchesReference(
                randomGlmDesign(5, taxa = 24, markers = 40, missingRate = 0.25),
                tuning = GlmFastTuning(dropped),
            )
        }
    }

    // --- workers, batching, replicates -----------------------------------------------------------------------------

    @Test
    fun anyNumberOfWorkersGivesTheSameBytes() = runTest {
        val input = randomGlmDesign(6, taxa = 60, markers = 90)
        val data = AssociationData.of(input)
        suspend fun run(workers: Int, blockSize: Int): List<String> {
            val collector = AssociationCollector(1.0, MaxPRule.Glm, data.traits)
            glm(
                data,
                GenotypeMarkerSource.of(data.data),
                GlmOptions(blockSize = blockSize),
                collector,
                workers = workers,
                tuning = GlmFastTuning(patternCacheDoubles = 1_000L * workers),
            )
            return collector.finish(AnalysisProvenance("glm")).text()
        }
        val one = run(1, 25)
        for ((workers, blockSize) in listOf(2 to 25, 3 to 7, 8 to 1, 16 to 30, 1 to 3, 40 to 90)) {
            assertEquals(one, run(workers, blockSize), "workers $workers, block size $blockSize")
        }
    }

    @Test
    fun traitsThatShareAnObservationSetShareTheirMarkerWork() = runTest {
        val input = randomGlmDesign(7)
        val data = AssociationData.of(input)
        assertTrue(data.observations(0).contentEquals(data.observations(1)))
        assertFalse(data.observations(0).contentEquals(data.observations(2)))
        val stats = GlmFastStats()
        val (_, together) = tables(input, stats = stats)
        // Two observation sets, so two fits per marker, not three.
        assertEquals(
            2L * 70,
            stats.fullSet + stats.indicators + stats.patternHits + stats.patternMisses + stats.skipped,
        )
        for (trait in listOf("y1", "y2", "y3")) {
            val (_, alone) = tables(input, AssociationSpec(traits = listOf(trait)))
            assertEquals(
                together.records.filter { it.trait == trait }.map { it.toString() }.toList(),
                alone.text().filter { !it.startsWith("AssociationEffect") },
            )
        }
    }

    @Test
    fun replicatedTraitsAreFittedByTheReferencePath() = runTest {
        val random = Random(8)
        val base = (0 until 30).map { "t$it" }
        val names = base + base.take(12)
        val calls = Array(12) { m ->
            Array(base.size) { if (random.nextDouble() < 0.1) "NN" else listOf("AA", "AC", "CC")[(it + m) % 3] }
        }
        val genotypes = InMemoryGenotypeTable.of(
            TaxaList.of(base),
            PositionList((0 until 12).map { Position(Chromosome.of(1), (it + 1) * 10, "r$it") }),
        ) { taxon, site -> Nucleotide.genotype(calls[site][taxon]) }
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(names),
                NumericAttribute(
                    "y",
                    AttributeType.data,
                    FloatArray(names.size) {
                        random.nextDouble(0.0, 9.0).toFloat()
                    },
                ),
                NumericAttribute(
                    "w",
                    AttributeType.data,
                    FloatArray(names.size) {
                        if (it <
                            30
                        ) {
                            it.toFloat()
                        } else {
                            Float.NaN
                        }
                    },
                ),
            ),
        )
        val stats = GlmFastStats()
        val (reference, fast) = tables(DataSet.of(genotypes, phenotype), stats = stats)
        assertSameResults(reference, fast)
        assertEquals(
            reference.text().filter { "trait=y," in it },
            fast.text().filter { "trait=y," in it },
            "the replicated trait's records are the reference's own",
        )
        assertEquals(12L, stats.reference, "y is replicated and goes to the reference path")
        assertTrue(
            fast.records.filter { it.trait == "y" && it.status == AssociationStatus.Ok }.all {
                it.message ==
                    TAXA_WITHIN_MARKER
            },
        )
        assertTrue(fast.records.none { it.trait == "w" && it.message.isNotEmpty() }, "w has one observation per taxon")
    }

    // --- the pieces -----------------------------------------------------------------------------------------------

    @Test
    fun thePatternCacheVerifiesKeysAndEvictsTheLeastRecentlyUsed() {
        val cache = PatternCache<String>(budget = 3)
        var built = 0
        fun get(pattern: IntArray, set: Int = 0) =
            cache.getOrPut(set, pattern, { 1L }) { "${pattern.toList()}@$set#${built++}" }
        // [0, 31] and [1, 0] have the same content hash; the second must not be served the first's value.
        assertEquals(intArrayOf(0, 31).contentHashCode(), intArrayOf(1, 0).contentHashCode())
        val first = get(intArrayOf(0, 31))
        val second = get(intArrayOf(1, 0))
        assertTrue(first != second)
        assertEquals(first, get(intArrayOf(0, 31)))
        assertEquals(1L, cache.hits)
        // The same pattern in another observation set is another entry.
        get(intArrayOf(0, 31), set = 1)
        assertEquals(3, cache.size)
        // Full: [1, 0] is least recently used and goes first.
        get(intArrayOf(2))
        assertEquals(3, cache.size)
        val misses = cache.misses
        get(intArrayOf(0, 31))
        assertEquals(misses, cache.misses, "[0, 31] was used recently and stayed")
        get(intArrayOf(1, 0))
        assertEquals(misses + 1, cache.misses, "[1, 0] was evicted")
        // A value bigger than the whole budget is built and returned but not kept.
        assertEquals(3, cache.size)
        cache.getOrPut(9, intArrayOf(9), { 4L }) { "big" }
        cache.getOrPut(9, intArrayOf(9), { 4L }) { "big" }
        assertEquals(3, cache.size)
        assertEquals(2L, cache.hits, "only the two hits above")
    }
}
