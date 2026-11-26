package net.maizegenetics.tassel.distance

import kotlinx.io.Source
import kotlinx.io.readLine
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.JavaFloat
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.math.abs

/**
 * How [DistanceMatrices.read] treats a file.
 *
 * @property symmetryTolerance The largest absolute difference allowed between the two copies of an
 *   off-diagonal value (row `i`, column `j` and row `j`, column `i`). A larger difference is an
 *   error naming both lines. `null` accepts any file and keeps the lower triangle, which is what
 *   TASSEL 5 does without saying so.
 */
data class DistanceMatrixReadOptions(val symmetryTolerance: Double? = DEFAULT_SYMMETRY_TOLERANCE) {
    init {
        require(symmetryTolerance == null || symmetryTolerance >= 0.0) {
            "symmetryTolerance must not be negative, not $symmetryTolerance"
        }
    }

    companion object {
        /** Absolute tolerance on the two copies of an off-diagonal value. */
        const val DEFAULT_SYMMETRY_TOLERANCE: Double = 1e-6
    }
}

/**
 * Reads TASSEL 5's square-matrix text, the format [DistanceMatrix.writeDelimited] writes:
 *
 * ```
 * ##Matrix_Type=Centered_IBS
 * 3
 * B73     1.2   0.1   -0.2
 * Mo17    0.1   0.9   0.0
 * W22    -0.2   0.0   1.1
 * ```
 *
 * Zero or more `##key=value` annotation lines, the number of taxa, then one line per taxon holding
 * its name and its whole row. Cells are separated by runs of whitespace, so tabs and spaces both
 * work, and a taxon name cannot contain whitespace. Values are parsed as Java's `Float.parseFloat`
 * parses them, so they are the same [Float]s TASSEL 5 reads. (TASSEL 5 parses a `Double` and
 * narrows it. The two differ only for text with more digits than a `Float` holds.)
 *
 * TASSEL 5 reads only full square rows, so a file holding only the lower triangle is an error
 * here, naming the line, as it is there.
 *
 * Where this departs from TASSEL 5 (ADR-0026):
 * - Duplicate taxon names are an error. TASSEL 5 would build a matrix it cannot index by name.
 * - A matrix whose two triangles disagree by more than [DistanceMatrixReadOptions.symmetryTolerance]
 *   is an error. TASSEL 5 silently keeps the lower triangle.
 * - A byte-order mark and blank lines are skipped, and every error names its line.
 * - An annotation without `=` is an error. TASSEL 5 fails with an array-index exception.
 *
 * @see net.maizegenetics.taxa.distance.ReadDistanceMatrix.readDistanceMatrix
 */
object DistanceMatrices {
    /**
     * Reads one matrix from [source], reporting one unit of [progress] per row.
     *
     * @throws IllegalArgumentException with the line number if the input is empty, has no count
     *   line, has the wrong number of rows or cells, repeats a taxon name, holds a cell that is not
     *   a number, or is asymmetric beyond [DistanceMatrixReadOptions.symmetryTolerance].
     */
    suspend fun read(
        source: Source,
        options: DistanceMatrixReadOptions = DistanceMatrixReadOptions(),
        progress: ProgressReporter = ProgressReporter.None,
    ): DistanceMatrix {
        var lineNumber = 0

        fun nextLine(): String? {
            while (true) {
                var text = source.readLine() ?: return null
                lineNumber++
                if (lineNumber == 1) text = text.removePrefix(BYTE_ORDER_MARK)
                text = text.trim { it <= ' ' }
                if (text.isNotEmpty()) return text
            }
        }

        val annotations = LinkedHashMap<String, String>()
        var countLine: String
        while (true) {
            countLine = nextLine() ?: throw IllegalArgumentException("Distance matrix has no matrix data")
            if (!countLine.startsWith("##")) break
            val body = countLine.substring(2)
            val equals = body.indexOf('=')
            if (equals < 0) {
                throw IllegalArgumentException("Line $lineNumber: annotation '$countLine' is not ##key=value")
            }
            annotations[body.substring(0, equals).trim()] = body.substring(equals + 1).trim()
        }
        val count = countLine.toIntOrNull()
            ?: throw IllegalArgumentException("Line $lineNumber: the number of taxa is not a number: '$countLine'")
        require(count >= 0) { "Line $lineNumber: the number of taxa is negative: $count" }
        val countLineNumber = lineNumber

        val names = ArrayList<String>(count)
        val lineOfRow = IntArray(count)
        val firstRowOf = HashMap<String, Int>(count)
        val values = FloatArray(DistanceMatrix.triangleSize(count))
        val tolerance = options.symmetryTolerance
        val tokens = ArrayList<String>(count + 1)

        for (row in 0 until count) {
            val line = nextLine() ?: throw IllegalArgumentException(
                "Line ${lineNumber + 1}: expected $count rows after the count on line $countLineNumber, found $row",
            )
            lineOfRow[row] = lineNumber
            splitWhitespace(line, tokens)
            if (tokens.size != count + 1) {
                throw IllegalArgumentException(
                    "Line $lineNumber: expected a name and $count values, found ${tokens.size - 1} values " +
                        "after '${tokens[0]}'",
                )
            }
            val name = tokens[0]
            firstRowOf.put(name, row)?.let { first ->
                throw IllegalArgumentException(
                    "Line $lineNumber: taxon '$name' repeats the one on line ${lineOfRow[first]}",
                )
            }
            names += name

            val base = row * (row + 1) / 2
            for (column in 0 until count) {
                val cell = tokens[column + 1]
                val value = JavaFloat.parse(cell)
                    ?: throw IllegalArgumentException("Line $lineNumber: '$cell' is not a number")
                when {
                    // Above the diagonal: parked in the mirrored slot, which row `column` checks and overwrites.
                    column > row -> values[column * (column + 1) / 2 + row] = value
                    column == row -> values[base + column] = value
                    else -> {
                        if (tolerance != null) {
                            checkSymmetry(
                                values[base + column],
                                value,
                                tolerance,
                                row,
                                column,
                                names,
                                lineOfRow,
                                lineNumber,
                            )
                        }
                        values[base + column] = value
                    }
                }
            }
            progress.report(row + 1, count, "Reading distance matrix")
        }
        nextLine()?.let { extra ->
            throw IllegalArgumentException(
                "Line $lineNumber: more rows than the $count on the count line ($countLineNumber): '${extra.take(
                    EXCERPT,
                )}'",
            )
        }
        return DistanceMatrix(TaxaList.of(names), values, annotations)
    }

    private fun checkSymmetry(
        upperCopy: Float,
        lowerCopy: Float,
        tolerance: Double,
        row: Int,
        column: Int,
        names: List<String>,
        lineOfRow: IntArray,
        lineNumber: Int,
    ) {
        if (upperCopy.isNaN() && lowerCopy.isNaN()) return
        if (upperCopy == lowerCopy) return
        if (abs(upperCopy.toDouble() - lowerCopy.toDouble()) <= tolerance) return
        throw IllegalArgumentException(
            "Line $lineNumber: the matrix is not symmetric; '${names[row]}' by '${names[column]}' is $lowerCopy " +
                "but '${names[column]}' by '${names[row]}' on line ${lineOfRow[column]} is $upperCopy",
        )
    }

    private fun splitWhitespace(line: String, into: MutableList<String>) {
        into.clear()
        var start = -1
        for (i in line.indices) {
            if (line[i] <= ' ') {
                if (start >= 0) into += line.substring(start, i)
                start = -1
            } else if (start < 0) {
                start = i
            }
        }
        if (start >= 0) into += line.substring(start)
    }

    private const val BYTE_ORDER_MARK = "\uFEFF"
    private const val EXCERPT = 40
}
