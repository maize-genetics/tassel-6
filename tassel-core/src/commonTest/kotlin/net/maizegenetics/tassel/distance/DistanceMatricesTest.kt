package net.maizegenetics.tassel.distance

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.writeString
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** [DistanceMatrices.read] against `golden/kinship/centered_ibs.txt` and hand-written files. */
class DistanceMatricesTest {
    @Test
    fun readsTheKinshipGoldenAndWritesItBackByteForByte() = runTest {
        val golden = normalizeLineEndings(goldenText("golden/kinship/centered_ibs.txt"))
        val matrix = read(golden)
        assertEquals(281, matrix.numberOfTaxa)
        assertEquals(golden, matrix.toDelimitedString())
    }

    @Test
    fun readsTheGoldensAnnotationsAndFloats() = runTest {
        val matrix = read(normalizeLineEndings(goldenText("golden/kinship/centered_ibs.txt")))
        assertEquals("Centered_IBS", matrix.annotations[DistanceMatrix.MATRIX_TYPE])
        assertEquals("918.7297856754799", matrix.annotations[DistanceMatrix.CENTERED_IBS_SUMPK])
        assertEquals("33-16", matrix.taxa[0].name)
        assertEquals(1.7978103f, matrix.distance(0, 0))
        assertEquals(0.03960289f, matrix.distance(0, 1))
        assertEquals(0.03960289f, matrix.distance(1, 0))
    }

    @Test
    fun acceptsTabsSpacesBlankLinesAndAByteOrderMark() = runTest {
        val matrix = read("\uFEFF##a = b \n\n2\n  A 1.5   0.25 \n\nB\t0.25\t2\n\n")
        assertEquals(listOf("A", "B"), matrix.taxa.names)
        assertEquals(mapOf("a" to "b"), matrix.annotations)
        assertEquals(0.25f, matrix.distance(1, 0))
        assertEquals(2f, matrix.distance(1, 1))
    }

    @Test
    fun parsesValuesAsJavaFloatsDo() = runTest {
        // 16777217 is halfway between two floats; Float.parseFloat rounds to even, 16777216.
        val matrix = read("1\nA\t16777217\n")
        assertEquals(16777216f, matrix.distance(0, 0))
        assertEquals(1f, read("1\nA\t1f\n").distance(0, 0))
    }

    @Test
    fun readsAZeroTaxonMatrix() = runTest {
        assertEquals(0, read("##k=v\n0\n").numberOfTaxa)
    }

    @Test
    fun reportsOneUnitOfProgressPerRow() = runTest {
        val reports = mutableListOf<Progress>()
        DistanceMatrices.read(source("2\nA 1 0\nB 0 1\n"), progress = { reports += it })
        assertEquals(listOf(1L, 2L), reports.map { it.completed })
        assertTrue(reports.all { it.total == 2L })
    }

    @Test
    fun rejectsDuplicateLabelsNamingBothLines() = runTest {
        assertError("Line 4", "'A'", "repeats the one on line 3") { read("##k=v\n2\nA 1 0\nA 0 1\n") }
    }

    @Test
    fun rejectsRaggedRows() = runTest {
        assertError("Line 4", "expected a name and 3 values, found 2") { read("##k=v\n3\nA 1 0 0\nB 0 1\nC 0 0 1\n") }
        assertError("Line 2", "found 3") { read("2\nA 1 0 0\nB 0 1\n") }
    }

    @Test
    fun rejectsALowerTriangleOnlyFile() = runTest {
        assertError("Line 2", "found 1 values") { read("3\nA 1\nB 0 1\nC 0 0 1\n") }
    }

    @Test
    fun rejectsNonNumbersNamingTheLine() = runTest {
        assertError("Line 4", "'x' is not a number") { read("##k=v\n2\nA 1 0\nB 0 x\n") }
        assertError("Line 2", "'0x1p3'") { read("1\nA 0x1p3\n") }
    }

    @Test
    fun rejectsABadCountAndWrongRowCounts() = runTest {
        assertError("Line 2", "not a number") { read("##k=v\ntwo\nA 1\n") }
        assertError("Line 1", "negative") { read("-1\n") }
        assertError("Line 4", "expected 3 rows") { read("3\nA 1 0 0\nB 0 1 0\n") }
        assertError("Line 4", "more rows") { read("2\nA 1 0\nB 0 1\nC 0 0\n") }
        assertError("no matrix data") { read("##k=v\n\n") }
        assertError("no matrix data") { read("") }
        assertError("Line 1", "not ##key=value") { read("##nokey\n1\nA 1\n") }
    }

    @Test
    fun rejectsAnAsymmetricMatrixUnlessToldNotTo() = runTest {
        val text = "2\nA 1 0.5\nB 0.4 1\n"
        assertError("Line 3", "not symmetric") { read(text) }
        // The lower triangle wins, as in TASSEL 5.
        val lenient = DistanceMatrices.read(source(text), DistanceMatrixReadOptions(symmetryTolerance = null))
        assertEquals(0.4f, lenient.distance(0, 1))
        val loose = DistanceMatrices.read(source(text), DistanceMatrixReadOptions(symmetryTolerance = 0.2))
        assertEquals(0.4f, loose.distance(1, 0))
    }

    @Test
    fun keepsNaNForValidationToReport() = runTest {
        val matrix = read("2\nA 1 NaN\nB NaN 1\n")
        assertTrue(matrix.distance(0, 1).isNaN())
    }

    private suspend fun read(text: String): DistanceMatrix = DistanceMatrices.read(source(text))

    private fun source(text: String) = Buffer().also { it.writeString(text) }

    private suspend fun assertError(vararg parts: String, block: suspend () -> Unit) {
        val error = assertFailsWith<IllegalArgumentException> { block() }
        parts.forEach { assertContains(error.message.orEmpty(), it) }
    }
}
