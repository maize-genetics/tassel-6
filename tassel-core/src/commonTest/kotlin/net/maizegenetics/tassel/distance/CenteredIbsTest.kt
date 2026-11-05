package net.maizegenetics.tassel.distance

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.testing.SquareMatrix
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGoldenMatrix
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
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
        assertEquals(expected, kinship().annotations[SUMPK])
    }

    /**
     * The matrix is in fact identical to the golden, not merely inside its tolerance: the `float`
     * sums land on the same bits and `floatToString` spells them the same way. The MANIFEST
     * tolerance stands as the contract ([matchesGoldenMatrix]); this test says how much room is
     * left under it, and will fail first if the summation order ever drifts.
     */
    @Test
    fun writesTheGoldenByteForByte() = runTest {
        assertEquals(normalizeLineEndings(goldenText(GOLDEN)), kinship().toDelimitedString())
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
        assertTrue(reports.size > 1, "expected a report per run of sites, got ${reports.size}")
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

    private companion object {
        const val GOLDEN = "golden/kinship/centered_ibs.txt"
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val TOLERANCE = 8e-6
        const val SUMPK = DistanceMatrix.CENTERED_IBS_SUMPK

        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }

        val kinship = SuspendLazy { centeredIbs(tutorial()) }

        val opposingHomozygotes = SuspendLazy { hapMap("AA" to "CC") }

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
