package net.maizegenetics.tassel.distance

import kotlinx.coroutines.cancel
import kotlinx.coroutines.launch
import kotlinx.coroutines.test.runTest
import kotlinx.coroutines.yield
import kotlinx.io.Buffer
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.SquareMatrix
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGoldenMatrix
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.coroutines.coroutineContext
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * Centered IBS kinship on the tutorial HapMap, against TASSEL 5's `KinshipPlugin` output in
 * `testdata/golden/kinship/centered_ibs.txt` (ADR-0009, tolerance `abs 8e-6`).
 *
 * @see net.maizegenetics.analysis.distance.CenteredIBSTest
 */
class CenteredIbsTest {
    @Test
    fun matchesGoldenMatrix() = runTest {
        assertMatchesGoldenMatrix(GOLDEN, SquareMatrix.parse(kinship().toDelimitedString()))
    }

    @Test
    fun matchesGoldenSumPk() = runTest {
        val expected = SquareMatrix.parse(goldenText(GOLDEN)).properties
            .single { it.startsWith("##$SUMPK=") }
            .substringAfter('=')
            .toDouble()
        val actual = kinship().annotations.getValue(SUMPK).toDouble()
        assertTrue(
            abs(expected - actual) <= SUMPK_TOLERANCE,
            "SumPk $actual differs from $expected by more than $SUMPK_TOLERANCE (ADR-0020)",
        )
    }

    /**
     * The retained TASSEL 5-order kernel is still byte-identical, including `SumPk`. The default
     * kernel is held to the MANIFEST tolerance ([matchesGoldenMatrix], ADR-0020).
     */
    @Test
    fun writesTheGoldenByteForByte() = runTest {
        val written = centeredIbsTassel5Order(tutorial(), workers = 1).toDelimitedString()
        assertEquals(normalizeLineEndings(goldenText(GOLDEN)), written)
    }

    @Test
    fun writesTassel5MatrixFormat() = runTest {
        val written = kinship().toDelimitedString().lines()
        val expected = normalizeLineEndings(goldenText(GOLDEN)).lines()
        assertEquals(expected.size, written.size, "line count")
        assertEquals(expected[0].substringBefore('='), written[0].substringBefore('='), "first annotation key")
        assertEquals(expected[1], written[1], "the Matrix_Type annotation line")
        assertEquals(expected[2], written[2], "the taxa count line")
        assertEquals(
            expected.drop(3).map { it.substringBefore('\t') },
            written.drop(3).map { it.substringBefore('\t') },
            "the row labels",
        )
        assertTrue(written.last().isEmpty(), "the file ends with a newline")
    }

    @Test
    fun writesBothAnnotationsSortedByKey() = runTest {
        assertEquals(listOf(SUMPK, DistanceMatrix.MATRIX_TYPE), kinship().annotations.keys.toList())
        assertEquals(KinshipMethod.Centered_IBS.name, kinship().annotations[DistanceMatrix.MATRIX_TYPE])
    }

    @Test
    fun isSymmetricAndSurvivesTheRoundTripThroughText() = runTest {
        val written = SquareMatrix.parse(kinship().toDelimitedString())
        assertEquals(281, kinship().numberOfTaxa)
        assertEquals(tutorial().taxa.names, written.labels)
        for (row in 0 until kinship().numberOfTaxa) {
            assertEquals(kinship().distance(row, 0), kinship().distance(0, row))
            assertEquals(kinship().distance(row, row), written.values[row][row].toFloat())
        }
    }

    @Test
    fun reportsProgressOverEverySite() = runTest {
        val reports = mutableListOf<Progress>()
        centeredIbs(tutorial(), progress = ProgressReporter { reports += it })
        val steps = listOf(reports.first().completed) + reports.zipWithNext { previous, next ->
            next.completed - previous.completed
        }
        // 3,093 sites, flushed about every 100 sites: 31 reports. The last step is the remainder.
        assertEquals(31, reports.size, "steps: $steps")
        assertTrue(steps.all { it in 1..150 }, "a step left the ~100 site cadence: $steps")
        assertTrue(reports.all { it.total == tutorial().numberOfSites.toLong() })
        assertEquals(tutorial().numberOfSites.toLong(), reports.last().completed)
        assertEquals(reports.map { it.completed }.sorted(), reports.map { it.completed })
    }

    @Test
    fun rejectsMaxAllelesOutsideTassel5Range() = runTest {
        assertFailsWith<IllegalArgumentException> { centeredIbs(opposingHomozygotes(), maxAlleles = 1) }
        assertFailsWith<IllegalArgumentException> { centeredIbs(opposingHomozygotes(), maxAlleles = 7) }
    }

    /**
     * Two taxa, one site, `AA` against `CC`. One allele is evaluated with frequency 0.5, so
     * `SumPk` is `2 * 0.5 * 0.5` and the terms are `(2 - 1)^2`, `(0 - 1)(2 - 1)`, and `(0 - 1)^2`.
     */
    @Test
    fun matchesHandComputedValuesOnTwoTaxa() = runTest {
        val matrix = centeredIbs(opposingHomozygotes())
        assertEquals("0.5", matrix.annotations[SUMPK])
        assertEquals(2.0f, matrix.distance(0, 0))
        assertEquals(2.0f, matrix.distance(1, 1))
        assertEquals(-2.0f, matrix.distance(0, 1))
        assertEquals(
            "##$SUMPK=0.5\n##Matrix_Type=Centered_IBS\n2\none\t2.0\t-2.0\ntwo\t-2.0\t2.0\n",
            matrix.toDelimitedString(),
        )
    }

    /** A site with one allele left is skipped: `alleles.size - 1` evaluates to zero pseudo-sites. */
    @Test
    fun monomorphicSitesContributeNothing() = runTest {
        val padded = centeredIbs(hapMap("AA" to "CC", "GG" to "GG", "TT" to "TT"))
        val plain = centeredIbs(opposingHomozygotes())
        assertEquals(plain.annotations, padded.annotations)
        assertEquals(plain.toDelimitedString(), padded.toDelimitedString())
    }

    @Test
    fun isRegisteredWithTassel5LegacyNames() = runTest {
        assertEquals(kinshipOperation, TasselOperations.registry.get("kinship"))
        assertEquals(kinshipOperation, TasselOperations.registry.byLegacyPlugin("KinshipPlugin"))
        assertEquals(KinshipMethod.Centered_IBS, kinshipOperation.parameter("method")?.default)
        assertEquals(DEFAULT_MAX_ALLELES, kinshipOperation.parameter("maxAlleles")?.default)
        assertEquals("method", kinshipOperation.parameterByLegacyName("method")?.name)
        assertEquals("maxAlleles", kinshipOperation.parameterByLegacyName("maxAlleles")?.name)
    }

    @Test
    fun parallelBandsMatchTheSequentialKernel() = runTest {
        val table = variedTable(taxa = 16, sites = 1_200)
        val sequential = centeredIbsTassel5Order(table, workers = 1)
        val parallel = centeredIbsTassel5Order(table, workers = 4)
        assertEquals(sequential.toDelimitedString(), parallel.toDelimitedString())
    }

    @Test
    fun parallelBandsMatchTheSequentialKernelOnTheTutorial() = runTest {
        val sequential = centeredIbsTassel5Order(tutorial(), workers = 1)
        val parallel = centeredIbsTassel5Order(tutorial(), workers = 4)
        assertEquals(sequential.toDelimitedString(), parallel.toDelimitedString())
    }

    @Test
    fun cancelsWhenTheCallerIsCancelled() = runTest {
        var seen = 0
        val job = launch {
            centeredIbs(
                variedTable(taxa = 16, sites = 250),
                progress = ProgressReporter {
                    seen++
                    coroutineContext.cancel()
                    yield()
                },
            )
        }
        job.join()
        assertTrue(job.isCancelled, "the run finished instead of cancelling")
        assertEquals(1, seen, "the kernel kept going after cancellation")
    }

    @Test
    fun matchesTheTassel5OrderOnVariedSites() = runTest {
        val table = variedTable(taxa = 12, sites = 80)
        for (maxAlleles in intArrayOf(2, DEFAULT_MAX_ALLELES)) {
            val expected = centeredIbsTassel5Order(table, maxAlleles = maxAlleles, workers = 1)
            val actual = centeredIbs(table, maxAlleles = maxAlleles)
            assertWithinTolerance(expected, actual)
        }
    }

    @Test
    fun pairBandsCoverEveryRowOnce() {
        for (taxa in listOf(1, 2, 5, 10, 281, 1_000)) {
            for (workers in listOf(1, 2, 3, 4, 8, 16)) {
                val bands = pairBands(taxa, workers)
                assertEquals(0, bands.first().first, "taxa=$taxa workers=$workers")
                var cursor = 0
                var pairs = 0
                for (band in bands) {
                    assertEquals(cursor, band.first, "taxa=$taxa workers=$workers bands=$bands")
                    assertFalse(band.isEmpty(), "empty band $band for taxa=$taxa workers=$workers")
                    pairs += pairCount(taxa, band)
                    cursor = band.last + 1
                }
                assertEquals(taxa, cursor)
                assertEquals(taxa * (taxa + 1) / 2, pairs)
            }
        }
    }

    private fun assertWithinTolerance(expected: DistanceMatrix, actual: DistanceMatrix) {
        assertEquals(expected.numberOfTaxa, actual.numberOfTaxa)
        val sumExpected = expected.annotations.getValue(SUMPK).toDouble()
        val sumActual = actual.annotations.getValue(SUMPK).toDouble()
        assertTrue(abs(sumExpected - sumActual) <= SUMPK_TOLERANCE, "SumPk $sumActual vs $sumExpected")
        for (row in 0 until expected.numberOfTaxa) {
            for (column in row until expected.numberOfTaxa) {
                val diff = abs(expected.distance(row, column) - actual.distance(row, column)).toDouble()
                assertTrue(diff <= TOLERANCE, "pair $row,$column differs by $diff")
            }
        }
    }

    private companion object {
        const val GOLDEN = "golden/kinship/centered_ibs.txt"
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val TOLERANCE = 8e-6
        const val SUMPK_TOLERANCE = 1e-4
        const val SUMPK = DistanceMatrix.CENTERED_IBS_SUMPK

        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }

        val kinship = SuspendLazy { centeredIbs(tutorial()) }

        val opposingHomozygotes = SuspendLazy { hapMap("AA" to "CC") }

        /**
         * Sites that exercise missing calls, a monomorphic column, and more than two alleles.
         * 1,200 sites is several TASSEL 5 blocks (`max(sites / 4, 1000)`).
         */
        fun variedTable(taxa: Int, sites: Int): GenotypeTable {
            val names = List(taxa) { "t$it" }
            val chromosome = Chromosome.of(1)
            val positions = List(sites) { Position(chromosome, it + 1) }
            val alleles = byteArrayOf(Allele.A, Allele.C, Allele.G, Allele.T)
            return InMemoryGenotypeTable.of(TaxaList.of(names), PositionList(positions)) { taxon, site ->
                when {
                    site % 11 == 0 -> Genotype.of(Allele.G, Allele.G)
                    taxon % 13 == 0 -> Genotype.UNKNOWN
                    site % 5 == 0 -> Genotype.of(alleles[(taxon + site) % 4], alleles[(taxon * 3 + site) % 4])
                    else -> {
                        val first = if ((taxon + site) % 3 == 0) Allele.C else Allele.A
                        val second = if ((taxon * 2 + site) % 3 == 0) Allele.C else Allele.A
                        Genotype.of(first, second)
                    }
                }
            }
        }

        /** One HapMap site per pair of genotypes, with two taxa named `one` and `two`. */
        suspend fun hapMap(vararg sites: Pair<String, String>): GenotypeTable {
            val header = (HapMap.FIXED_COLUMNS + listOf("one", "two")).joinToString("\t")
            val lines = sites.mapIndexed { site, (one, two) ->
                listOf("s$site", "A/C", "1", "${site + 1}", "+", "NA", "NA", "NA", "NA", "NA", "NA", one, two)
                    .joinToString("\t")
            }
            val text = (listOf(header) + lines).joinToString("\n", postfix = "\n")
            return HapMap.read(Buffer().also { it.write(text.encodeToByteArray()) })
        }
    }
}
