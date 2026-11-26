package net.maizegenetics.tassel.math

import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertSame

class DenseMatrixTest {
    @Test
    fun storageIsColumnMajor() {
        val matrix = DenseMatrix.ofRows(doubleArrayOf(1.0, 2.0, 3.0), doubleArrayOf(4.0, 5.0, 6.0))
        assertEquals(2, matrix.rows)
        assertEquals(3, matrix.columns)
        assertContentEquals(doubleArrayOf(1.0, 4.0, 2.0, 5.0, 3.0, 6.0), matrix.data)
        assertEquals(6.0, matrix[1, 2])
        assertEquals(4, matrix.columnOffset(2))
    }

    @Test
    fun wrapsTheArrayItIsGiven() {
        val data = DoubleArray(6)
        val matrix = DenseMatrix(3, 2, data)
        matrix[2, 1] = 7.0
        assertSame(data, matrix.data)
        assertEquals(7.0, data[5])
    }

    @Test
    fun columnViewsShareStorage() {
        val matrix = DenseMatrix.ofRows(doubleArrayOf(1.0, 2.0), doubleArrayOf(3.0, 4.0), doubleArrayOf(5.0, 6.0))
        val column = matrix.column(1)
        assertEquals(3, column.length)
        assertEquals(3, column.offset)
        assertEquals(4.0, column[1])
        column[1] = 40.0
        assertEquals(40.0, matrix[1, 1])
        assertContentEquals(doubleArrayOf(2.0, 40.0, 6.0), matrix.columnCopy(1))
        assertContentEquals(doubleArrayOf(2.0, 40.0, 6.0), column.toDoubleArray())
    }

    @Test
    fun copyFromReusesStorageAndCopyDoesNot() {
        val source = DenseMatrix.ofRows(doubleArrayOf(1.0, 2.0), doubleArrayOf(3.0, 4.0))
        val target = DenseMatrix(2, 2)
        val storage = target.data
        target.copyFrom(source)
        assertSame(storage, target.data)
        assertContentEquals(source.data, target.data)
        val independent = source.copy()
        independent[0, 0] = 9.0
        assertEquals(1.0, source[0, 0])
        assertFailsWith<IllegalArgumentException> { target.copyFrom(DenseMatrix(3, 2)) }
    }

    @Test
    fun transposeAndIdentity() {
        val matrix = DenseMatrix.ofRows(doubleArrayOf(1.0, 2.0, 3.0), doubleArrayOf(4.0, 5.0, 6.0))
        val t = matrix.transpose()
        assertEquals(3, t.rows)
        assertEquals(2, t.columns)
        assertEquals(matrix[1, 2], t[2, 1])
        assertContentEquals(doubleArrayOf(1.0, 0.0, 0.0, 1.0), DenseMatrix.identity(2).data)
    }

    @Test
    fun rejectsBadShapesAndIndices() {
        assertFailsWith<IllegalArgumentException> { DenseMatrix(2, 2, DoubleArray(3)) }
        assertFailsWith<IllegalArgumentException> { DenseMatrix(-1, 2) }
        assertFailsWith<IllegalArgumentException> { DenseMatrix.ofRows(doubleArrayOf(1.0), doubleArrayOf(1.0, 2.0)) }
        val matrix = DenseMatrix(2, 2)
        assertFailsWith<IllegalArgumentException> { matrix[2, 0] }
        assertFailsWith<IllegalArgumentException> { matrix[0, 2] = 1.0 }
        assertFailsWith<IllegalArgumentException> { matrix.column(2) }
    }

    @Test
    fun emptyMatricesAreAllowed() {
        assertEquals(0, DenseMatrix(0, 3).data.size)
        assertEquals(0, DenseMatrix(3, 0).data.size)
    }
}
