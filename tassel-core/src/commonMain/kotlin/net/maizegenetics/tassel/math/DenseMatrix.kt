package net.maizegenetics.tassel.math

/**
 * A dense real matrix in column-major order: element `(row, column)` is `data[column * rows + row]`, so each column
 * is a contiguous run of [rows] values. Householder reflections, triangular solves, and `X'X` all walk down columns,
 * which makes that the cache-friendly layout (ADR-0025). [SymmetricEigenDecomposition] stores its eigenvectors the
 * same way.
 *
 * The matrix wraps [data] rather than copying it, so a caller that owns a workspace array can view it as a matrix
 * for free, and a hot loop can reuse one matrix across many fits with [copyFrom]. Nothing here allocates per
 * element access.
 *
 * Values are not checked for finiteness; the decompositions that need finite input say so.
 */
class DenseMatrix(
    /** Number of rows. */
    val rows: Int,
    /** Number of columns. */
    val columns: Int,
    /** The values, column-major. Shared, not copied. */
    val data: DoubleArray,
) {
    init {
        require(rows >= 0 && columns >= 0) { "dimensions must not be negative, got $rows x $columns" }
        require(rows.toLong() * columns <= Int.MAX_VALUE) { "a $rows x $columns matrix does not fit in one array" }
        require(data.size == rows * columns) {
            "expected ${rows * columns} values for $rows x $columns, got ${data.size}"
        }
    }

    /** A zero-filled `rows x columns` matrix. */
    constructor(rows: Int, columns: Int) : this(rows, columns, DoubleArray(checkedSize(rows, columns)))

    /** Index in [data] of the first element of [column]. */
    fun columnOffset(column: Int): Int {
        requireColumn(column)
        return column * rows
    }

    operator fun get(row: Int, column: Int): Double {
        requireRow(row)
        requireColumn(column)
        return data[column * rows + row]
    }

    operator fun set(row: Int, column: Int, value: Double) {
        requireRow(row)
        requireColumn(column)
        data[column * rows + row] = value
    }

    /** A copy of [column]. */
    fun columnCopy(column: Int): DoubleArray {
        val from = columnOffset(column)
        return data.copyOfRange(from, from + rows)
    }

    /** A read-write window onto [column] that shares [data]; see [DenseColumn]. */
    fun column(column: Int): DenseColumn = DenseColumn(data, columnOffset(column), rows)

    /** An independent copy of this matrix. */
    fun copy(): DenseMatrix = DenseMatrix(rows, columns, data.copyOf())

    /** Overwrites this matrix with [other], which must have the same shape. Reuses this matrix's storage. */
    fun copyFrom(other: DenseMatrix) {
        require(other.rows == rows && other.columns == columns) {
            "cannot copy a ${other.rows} x ${other.columns} matrix into a $rows x $columns one"
        }
        other.data.copyInto(data)
    }

    /** The transpose, as a new matrix. */
    fun transpose(): DenseMatrix {
        val result = DenseMatrix(columns, rows)
        for (c in 0 until columns) {
            for (r in 0 until rows) result.data[r * columns + c] = data[c * rows + r]
        }
        return result
    }

    override fun toString(): String = "DenseMatrix($rows x $columns)"

    private fun requireRow(row: Int) = require(row in 0 until rows) { "row $row is outside 0 until $rows" }

    private fun requireColumn(column: Int) =
        require(column in 0 until columns) { "column $column is outside 0 until $columns" }

    companion object {
        /** Builds a matrix from rows written the way a person reads them. Every row must have the same length. */
        fun ofRows(vararg rows: DoubleArray): DenseMatrix {
            val columns = if (rows.isEmpty()) 0 else rows[0].size
            rows.forEachIndexed { index, row ->
                require(row.size == columns) { "row $index has ${row.size} values, expected $columns" }
            }
            val matrix = DenseMatrix(rows.size, columns)
            for (r in rows.indices) {
                for (c in 0 until columns) matrix.data[c * rows.size + r] = rows[r][c]
            }
            return matrix
        }

        /** The `size x size` identity. */
        fun identity(size: Int): DenseMatrix {
            val matrix = DenseMatrix(size, size)
            for (i in 0 until size) matrix.data[i * size + i] = 1.0
            return matrix
        }

        private fun checkedSize(rows: Int, columns: Int): Int {
            require(rows >= 0 && columns >= 0) { "dimensions must not be negative, got $rows x $columns" }
            require(rows.toLong() * columns <= Int.MAX_VALUE) { "a $rows x $columns matrix does not fit in one array" }
            return rows * columns
        }
    }
}

/**
 * A window of [length] values starting at [offset] in [array], used to name one matrix column without copying it.
 * Writes go through to the matrix.
 */
class DenseColumn internal constructor(
    /** The backing array, shared with the matrix. */
    val array: DoubleArray,
    /** Index of this column's first value in [array]. */
    val offset: Int,
    /** Number of values. */
    val length: Int,
) {
    operator fun get(index: Int): Double {
        require(index in 0 until length) { "index $index is outside 0 until $length" }
        return array[offset + index]
    }

    operator fun set(index: Int, value: Double) {
        require(index in 0 until length) { "index $index is outside 0 until $length" }
        array[offset + index] = value
    }

    /** A copy of the window. */
    fun toDoubleArray(): DoubleArray = array.copyOfRange(offset, offset + length)
}
