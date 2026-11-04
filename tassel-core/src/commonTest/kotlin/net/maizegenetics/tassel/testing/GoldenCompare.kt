// Comparisons against the TASSEL 5 goldens in testdata/golden/ (ADR-0009), driven by the tolerance
// that testdata/golden/MANIFEST records for each file.
//
// Every ported analysis debugs against these failures, so a failure names the file, the line, the
// row label, and the column, and shows the first few differing cells with their difference and the
// tolerance that was applied.
//
//     assertMatchesGolden("golden/summary/alleles.txt", table.toDelimitedString())
//     assertMatchesGoldenMatrix("golden/kinship/centered_ibs.txt", SquareMatrix(taxa, kinship))

package net.maizegenetics.tassel.testing

import kotlin.math.abs
import kotlin.test.fail

/** `testdata/golden/MANIFEST`, parsed once per test run. */
val goldenManifest: GoldenManifest by lazy { GoldenManifest.load() }

/** The tolerance the MANIFEST records for [goldenPath], which is relative to `testdata/`. */
fun goldenTolerance(goldenPath: String): Tolerance = goldenManifest.tolerance(goldenPath)

/** The text of the golden at [goldenPath], which is relative to `testdata/`. */
fun goldenText(goldenPath: String): String = TestData.readText(goldenPath)

/**
 * Asserts that delimited [actual] matches the golden at [goldenPath] (relative to `testdata/`),
 * under the tolerance the MANIFEST records for it. See [assertDelimitedMatches].
 */
fun assertMatchesGolden(
    goldenPath: String,
    actual: String,
    tolerance: Tolerance = goldenTolerance(goldenPath),
    delimiter: String = "\t",
) = assertDelimitedMatches(goldenText(goldenPath), actual, tolerance, delimiter, goldenPath)

/**
 * Asserts that delimited [actual] matches [expected] under [tolerance]. [source] names the expected
 * text in failures, and is the golden path when there is one.
 *
 * [Tolerance.Exact] compares the whole text after normalizing line endings, and nothing else. The
 * others compare cell by cell: the header line, the row count, and the column count must be
 * identical, a cell that is not a number must be identical, and a numeric cell must agree within
 * the tolerance. [Tolerance.Format] checks everything except the numbers. `up-to-sign` first
 * negates any wholly numeric column whose largest-magnitude entry is negative, independently on
 * each side, which is how PCA eigenvectors are compared.
 */
fun assertDelimitedMatches(
    expected: String,
    actual: String,
    tolerance: Tolerance,
    delimiter: String = "\t",
    source: String = "expected output",
) {
    val expectedText = normalizeLineEndings(expected)
    val actualText = normalizeLineEndings(actual)
    if (tolerance is Tolerance.Exact) {
        assertTextMatches(expectedText, actualText, source)
        return
    }

    val expectedRows = splitRows(expectedText, delimiter)
    val actualRows = splitRows(actualText, delimiter)
    if (expectedRows.isEmpty() || actualRows.isEmpty()) {
        fail("$source has ${expectedRows.size} lines and actual has ${actualRows.size}; both need a header")
    }
    assertHeaderMatches(expectedRows[0], actualRows[0], source)
    if (expectedRows.size != actualRows.size) {
        fail("$source has ${expectedRows.size - 1} data rows, actual has ${actualRows.size - 1}")
    }
    val columns = expectedRows[0].size
    for (row in 1 until actualRows.size) {
        val cells = actualRows[row].size
        if (cells != columns) {
            fail("$source ${rowName(actualRows, row)} has $cells cells, expected $columns")
        }
    }

    val expectedSigns = columnSigns(expectedRows, tolerance)
    val actualSigns = columnSigns(actualRows, tolerance)
    val mismatches = Mismatches(source, tolerance)
    for (row in 1 until expectedRows.size) {
        for (column in 0 until columns) {
            val where = "${rowName(expectedRows, row)}, column ${column + 1} (\"${expectedRows[0][column]}\")"
            val problem = cellProblem(
                expected = expectedRows[row][column],
                actual = actualRows[row][column],
                expectedSign = expectedSigns[column],
                actualSign = actualSigns[column],
                tolerance = tolerance,
            )
            mismatches.check(where, problem)
        }
    }
    mismatches.assertNone()
}

/**
 * A square matrix of doubles with one label per row and column, as TASSEL 5 writes a distance or
 * kinship matrix: `#` property lines, a line holding the row count, then one `label<tab>value...`
 * line per row.
 */
class SquareMatrix(
    val labels: List<String>,
    val values: List<DoubleArray>,
    /** The leading `#` lines as written. [parse] fills these in; comparison ignores them. */
    val properties: List<String> = emptyList(),
) {
    init {
        require(values.size == labels.size) {
            "square matrix has ${labels.size} labels and ${values.size} rows"
        }
        values.forEachIndexed { row, cells ->
            require(cells.size == labels.size) {
                "square matrix row $row has ${cells.size} values, expected ${labels.size}"
            }
        }
    }

    val size: Int get() = labels.size

    companion object {
        /** Parses TASSEL 5's square-matrix text. */
        fun parse(text: String): SquareMatrix {
            val lines = normalizeLineEndings(text).lines().filter { it.isNotBlank() }
            val properties = lines.takeWhile { it.startsWith("#") }
            val body = lines.drop(properties.size)
            val size = body.firstOrNull()?.trim()?.toIntOrNull()
                ?: fail("square matrix: expected a row count, found \"${body.firstOrNull()}\"")
            val rows = body.drop(1)
            if (rows.size != size) {
                fail("square matrix: the count line says $size rows, found ${rows.size}")
            }
            val labels = mutableListOf<String>()
            val values = mutableListOf<DoubleArray>()
            rows.forEachIndexed { row, line ->
                val cells = line.split('\t').filter { it.isNotEmpty() }
                if (cells.size != size + 1) {
                    fail("square matrix row ${row + 1} has ${cells.size - 1} values, expected $size")
                }
                labels += cells[0]
                values += DoubleArray(size) { cells[it + 1].toDouble() }
            }
            return SquareMatrix(labels, values, properties)
        }
    }
}

/**
 * Asserts that [actual] matches the golden square matrix at [goldenPath] (relative to `testdata/`),
 * under the tolerance the MANIFEST records for it. See [assertMatrixMatches].
 */
fun assertMatchesGoldenMatrix(
    goldenPath: String,
    actual: SquareMatrix,
    tolerance: Tolerance = goldenTolerance(goldenPath),
) = assertMatrixMatches(SquareMatrix.parse(goldenText(goldenPath)), actual, tolerance, goldenPath)

/**
 * Asserts that square matrix [actual] matches [expected] under [tolerance]: the same labels in the
 * same order, then every value within the tolerance. [Tolerance.Format] checks the labels only, and
 * `up-to-sign` negates a column whose largest-magnitude entry is negative.
 */
fun assertMatrixMatches(
    expected: SquareMatrix,
    actual: SquareMatrix,
    tolerance: Tolerance,
    source: String = "expected matrix",
) {
    if (expected.size != actual.size) {
        fail("$source is ${expected.size} by ${expected.size}, actual is ${actual.size} by ${actual.size}")
    }
    val label = expected.labels.indices.firstOrNull { expected.labels[it] != actual.labels[it] }
    if (label != null) {
        fail(
            "$source row ${label + 1} is labelled \"${expected.labels[label]}\", " +
                "actual \"${actual.labels[label]}\"",
        )
    }
    if (tolerance is Tolerance.Format) return

    val expectedSigns = matrixColumnSigns(expected, tolerance)
    val actualSigns = matrixColumnSigns(actual, tolerance)
    val mismatches = Mismatches(source, tolerance)
    for (row in 0 until expected.size) {
        for (column in 0 until expected.size) {
            val where =
                "row ${row + 1} (\"${expected.labels[row]}\"), " +
                    "column ${column + 1} (\"${expected.labels[column]}\")"
            val problem = numericProblem(
                expected.values[row][column] * expectedSigns[column],
                actual.values[row][column] * actualSigns[column],
                tolerance,
            )
            mismatches.check(where, problem)
        }
    }
    mismatches.assertNone()
}

/** [text] with `\r\n` and `\r` endings rewritten as `\n`; nothing else changes. */
fun normalizeLineEndings(text: String): String = text.replace("\r\n", "\n").replace('\r', '\n')

// --- internals --------------------------------------------------------------------------------

/** Collects differing cells, so a failure can show the first few and the total. */
private class Mismatches(private val source: String, private val tolerance: Tolerance) {
    private val problems = mutableListOf<String>()
    private var compared = 0

    fun check(where: String, problem: String?) {
        compared++
        if (problem != null) problems += "$where: $problem"
    }

    fun assertNone() {
        if (problems.isEmpty()) return
        val shown = problems.take(SHOWN)
        val hidden = problems.size - shown.size
        fail(
            buildString {
                append(source)
                append(": ${problems.size} of $compared cells differ (tolerance ${describe(tolerance)})")
                shown.forEach { append("\n  ").append(it) }
                if (hidden > 0) append("\n  ...and $hidden more")
            },
        )
    }

    private companion object {
        const val SHOWN = 5
    }
}

private fun assertTextMatches(expected: String, actual: String, source: String) {
    if (expected == actual) return
    val expectedLines = expected.lines()
    val actualLines = actual.lines()
    val shared = minOf(expectedLines.size, actualLines.size)
    val line = (0 until shared).firstOrNull { expectedLines[it] != actualLines[it] }
    val detail = if (line == null) {
        val longer = if (expectedLines.size > shared) "expected" else "actual"
        "the first $shared lines match, then only the $longer text continues"
    } else {
        "line ${line + 1} differs:" +
            "\n  expected: ${excerpt(expectedLines[line])}" +
            "\n  actual:   ${excerpt(actualLines[line])}"
    }
    fail(
        "$source does not match exactly (${expectedLines.size} expected lines, " +
            "${actualLines.size} actual): $detail",
    )
}

private fun assertHeaderMatches(expected: List<String>, actual: List<String>, source: String) {
    if (expected == actual) return
    val column = (0 until minOf(expected.size, actual.size)).firstOrNull { expected[it] != actual[it] }
    fail(
        if (column != null) {
            "$source header column ${column + 1} is \"${expected[column]}\", actual \"${actual[column]}\""
        } else {
            "$source header has ${expected.size} columns, actual has ${actual.size}: " +
                "expected [${expected.joinToString()}], actual [${actual.joinToString()}]"
        },
    )
}

private fun splitRows(text: String, delimiter: String): List<List<String>> =
    text.lines().dropLastWhile { it.isEmpty() }.map { it.split(delimiter) }

/** `line 12 (row "PZB00859.1")`: the file line number, which is one more than the data row index. */
private fun rowName(rows: List<List<String>>, row: Int): String =
    "line ${row + 1} (row \"${rows[row].firstOrNull().orEmpty()}\")"

/**
 * The sign each column is multiplied by before comparing: `-1.0` for a wholly numeric column whose
 * largest-magnitude entry is negative when [tolerance] is `up-to-sign`, and `1.0` otherwise.
 */
private fun columnSigns(rows: List<List<String>>, tolerance: Tolerance): DoubleArray {
    val signs = DoubleArray(rows[0].size) { 1.0 }
    if (!tolerance.upToSign) return signs
    for (column in signs.indices) {
        var largest = 0.0
        var sign = 1.0
        var numeric = true
        for (row in 1 until rows.size) {
            val value = rows[row][column].toDoubleOrNull()
            if (value == null || value.isNaN()) {
                numeric = false
                break
            }
            if (abs(value) > largest) {
                largest = abs(value)
                sign = if (value < 0.0) -1.0 else 1.0
            }
        }
        if (numeric) signs[column] = sign
    }
    return signs
}

private fun matrixColumnSigns(matrix: SquareMatrix, tolerance: Tolerance): DoubleArray {
    val signs = DoubleArray(matrix.size) { 1.0 }
    if (!tolerance.upToSign) return signs
    for (column in signs.indices) {
        var largest = 0.0
        for (row in 0 until matrix.size) {
            val value = matrix.values[row][column]
            if (abs(value) > largest) {
                largest = abs(value)
                signs[column] = if (value < 0.0) -1.0 else 1.0
            }
        }
    }
    return signs
}

private fun cellProblem(
    expected: String,
    actual: String,
    expectedSign: Double,
    actualSign: Double,
    tolerance: Tolerance,
): String? {
    val expectedValue = expected.toDoubleOrNull()
    val actualValue = actual.toDoubleOrNull()
    return when {
        expectedValue == null && actualValue == null ->
            if (expected == actual) null else "expected the text \"$expected\", actual \"$actual\""
        expectedValue == null || actualValue == null ->
            "expected \"$expected\", actual \"$actual\"; one is a number and the other is not"
        tolerance is Tolerance.Format -> null
        else -> numericProblem(expectedValue * expectedSign, actualValue * actualSign, tolerance)
    }
}

private fun numericProblem(expected: Double, actual: Double, tolerance: Tolerance): String? {
    if (expected == actual) return null
    if (expected.isNaN() || actual.isNaN()) {
        return if (expected.isNaN() && actual.isNaN()) null else "expected $expected, actual $actual"
    }
    val difference = abs(expected - actual)
    return when (tolerance) {
        is Tolerance.Absolute ->
            if (difference <= tolerance.value) {
                null
            } else {
                "expected $expected, actual $actual (off by $difference, allowed ${tolerance.value})"
            }
        is Tolerance.Relative -> {
            val allowed = tolerance.value * abs(expected)
            if (expected != 0.0 && difference <= allowed) {
                null
            } else {
                "expected $expected, actual $actual (off by $difference, allowed $allowed)"
            }
        }
        else -> null
    }
}

/** The MANIFEST spelling of [tolerance], for failure messages. */
private fun describe(tolerance: Tolerance): String {
    val text = when (tolerance) {
        Tolerance.Exact -> "exact"
        Tolerance.Format -> "format"
        is Tolerance.Absolute -> "abs ${tolerance.value}"
        is Tolerance.Relative -> "rel ${tolerance.value}"
    }
    return if (tolerance.upToSign) "$text up-to-sign" else text
}

private val Tolerance.upToSign: Boolean
    get() = when (this) {
        is Tolerance.Absolute -> upToSign
        is Tolerance.Relative -> upToSign
        else -> false
    }

private fun excerpt(line: String): String = if (line.length <= EXCERPT) line else line.take(EXCERPT) + "..."

private const val EXCERPT = 160
