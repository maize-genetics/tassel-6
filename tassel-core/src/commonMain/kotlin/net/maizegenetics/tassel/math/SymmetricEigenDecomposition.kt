package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.sqrt

/**
 * Eigenvalues and eigenvectors of a dense real symmetric matrix.
 *
 * The decomposition is `A = V * diag(eigenvalues) * V'`, computed by Householder reduction to tridiagonal form
 * (`tred2`) followed by the implicit QL algorithm with shifts (`tql2`). Both come from the Algol procedures of
 * Bowdler, Martin, Reinsch, and Wilkinson (*Handbook for Automatic Computation*, Vol. ii, Linear Algebra) and the
 * EISPACK Fortran that followed. This is a port of the public-domain JAMA `EigenvalueDecomposition` symmetric path
 * (NIST/MathWorks, JAMA 1.0.3), which is the same algorithm TASSEL 5 runs through its Colt-derived
 * `SymmetricEigenvalueDecomposition`. No part of JAMA or EJML is linked; the algorithm is reimplemented.
 *
 * ## Output contract
 *
 * - [eigenvalues] are sorted descending. Ties keep the order the QL iteration produced; the result is a
 *   deterministic function of the input on every target.
 * - Eigenvectors are the columns of `V`: column `j`, read with [eigenvector] or `get(row, j)`, belongs to
 *   `eigenvalue(j)`.
 * - Every column has unit Euclidean length.
 * - Sign convention: in each column the entry with the largest absolute value is made positive. If several entries
 *   tie exactly in absolute value, the one with the smallest row index decides. Columns of a degenerate
 *   (repeated-eigenvalue) subspace are orthonormal but otherwise arbitrary, as any eigensolver's are.
 *
 * Only the lower triangle of the input, including the diagonal, is read; the upper triangle is ignored. A caller
 * that accumulates a Gram matrix can therefore fill the lower triangle only. The input must be finite.
 *
 * Working memory is one `size * size` array plus two vectors of length `size`, so a 281 x 281 taxa Gram matrix
 * costs about 630 kB and nothing is allocated per site.
 *
 * Round-off can push eigenvalues of a positive-semidefinite matrix slightly below zero; callers that need
 * non-negative values (variance proportions, for example) should clamp them.
 *
 * @see net.maizegenetics.matrixalgebra.decomposition.SymmetricEigenvalueDecomposition
 */
class SymmetricEigenDecomposition private constructor(
    /** Number of rows (and columns) of the decomposed matrix. */
    val size: Int,
    private val values: DoubleArray,
    // Column-major: element (row, column) of V is at vectors[column * size + row]. The QL sweeps and the
    // Householder updates both run down a column, so columns are contiguous.
    private val vectors: DoubleArray,
) {
    /** The eigenvalue at [index], counting from the largest. */
    fun eigenvalue(index: Int): Double {
        requireIndex(index)
        return values[index]
    }

    /** A copy of all eigenvalues, largest first. */
    fun eigenvalues(): DoubleArray = values.copyOf()

    /** A copy of the unit-length eigenvector belonging to `eigenvalue(index)`. */
    fun eigenvector(index: Int): DoubleArray {
        requireIndex(index)
        val from = index * size
        return vectors.copyOfRange(from, from + size)
    }

    /** Element [row] of the eigenvector belonging to `eigenvalue(column)`; that is, `V[row, column]`. */
    operator fun get(row: Int, column: Int): Double {
        requireIndex(row)
        requireIndex(column)
        return vectors[column * size + row]
    }

    private fun requireIndex(index: Int) {
        require(index in 0 until size) { "index $index is outside 0 until $size" }
    }

    companion object {
        // 2^-52, the gap between 1.0 and the next Double. A literal keeps it identical on every target.
        private const val EPSILON = 2.220446049250313E-16

        // The shifted QL iteration needs a handful of sweeps per eigenvalue; this only guards against a
        // non-converging input that slipped past the finiteness check.
        private const val MAX_ITERATIONS = 50

        /** Decomposes [matrix], a `size x size` matrix of rows, reading its lower triangle. */
        fun of(matrix: Array<DoubleArray>): SymmetricEigenDecomposition {
            val size = matrix.size
            require(size > 0) { "matrix must have at least one row" }
            matrix.forEachIndexed { row, values ->
                require(values.size == size) { "row $row has ${values.size} values, expected $size" }
            }
            return of(size) { row, column -> matrix[row][column] }
        }

        /**
         * Decomposes a `size x size` matrix held in [values] as one array, reading its lower triangle. Row-major
         * and column-major orders agree for a symmetric matrix, so either may be passed.
         */
        fun of(size: Int, values: DoubleArray): SymmetricEigenDecomposition {
            require(values.size == size * size) {
                "expected ${size * size} values for a $size x $size matrix, got ${values.size}"
            }
            return of(size) { row, column -> values[row * size + column] }
        }

        /**
         * Decomposes the `size x size` matrix whose entries [element] returns. [element] is called once per entry
         * of the lower triangle, so a caller holding only that triangle need not mirror it.
         */
        fun of(size: Int, element: (row: Int, column: Int) -> Double): SymmetricEigenDecomposition {
            require(size > 0) { "size must be positive, got $size" }
            val vectors = DoubleArray(size * size)
            for (column in 0 until size) {
                for (row in column until size) {
                    val value = element(row, column)
                    require(value.isFinite()) { "element ($row, $column) is $value; the matrix must be finite" }
                    vectors[column * size + row] = value
                    vectors[row * size + column] = value
                }
            }
            val diagonal = DoubleArray(size)
            val offDiagonal = DoubleArray(size)
            tridiagonalize(vectors, size, diagonal, offDiagonal)
            diagonalize(vectors, size, diagonal, offDiagonal)
            sortDescending(vectors, size, diagonal)
            normalize(vectors, size)
            return SymmetricEigenDecomposition(size, diagonal, vectors)
        }

        /**
         * Householder reduction of [v] to tridiagonal form, leaving the accumulated transformation in [v], the
         * diagonal in [d], and the sub-diagonal in `e[1 until n]`. JAMA's `tred2`.
         */
        private fun tridiagonalize(v: DoubleArray, n: Int, d: DoubleArray, e: DoubleArray) {
            for (j in 0 until n) d[j] = v[j * n + n - 1]

            for (i in n - 1 downTo 1) {
                // Scale to avoid under- and overflow.
                var scale = 0.0
                var h = 0.0
                for (k in 0 until i) scale += abs(d[k])
                if (scale == 0.0) {
                    e[i] = d[i - 1]
                    for (j in 0 until i) {
                        d[j] = v[j * n + i - 1]
                        v[j * n + i] = 0.0
                        v[i * n + j] = 0.0
                    }
                } else {
                    // Generate the Householder vector.
                    for (k in 0 until i) {
                        d[k] /= scale
                        h += d[k] * d[k]
                    }
                    var f = d[i - 1]
                    var g = sqrt(h)
                    if (f > 0) g = -g
                    e[i] = scale * g
                    h -= f * g
                    d[i - 1] = f - g
                    for (j in 0 until i) e[j] = 0.0

                    // Apply the similarity transformation to the remaining columns.
                    for (j in 0 until i) {
                        f = d[j]
                        v[i * n + j] = f
                        val column = j * n
                        g = e[j] + v[column + j] * f
                        for (k in j + 1 until i) {
                            g += v[column + k] * d[k]
                            e[k] += v[column + k] * f
                        }
                        e[j] = g
                    }
                    f = 0.0
                    for (j in 0 until i) {
                        e[j] /= h
                        f += e[j] * d[j]
                    }
                    val hh = f / (h + h)
                    for (j in 0 until i) e[j] -= hh * d[j]
                    for (j in 0 until i) {
                        f = d[j]
                        g = e[j]
                        val column = j * n
                        for (k in j until i) v[column + k] -= f * e[k] + g * d[k]
                        d[j] = v[column + i - 1]
                        v[column + i] = 0.0
                    }
                }
                d[i] = h
            }

            // Accumulate the transformations.
            for (i in 0 until n - 1) {
                v[i * n + n - 1] = v[i * n + i]
                v[i * n + i] = 1.0
                val h = d[i + 1]
                val next = (i + 1) * n
                if (h != 0.0) {
                    for (k in 0..i) d[k] = v[next + k] / h
                    for (j in 0..i) {
                        var g = 0.0
                        val column = j * n
                        for (k in 0..i) g += v[next + k] * v[column + k]
                        for (k in 0..i) v[column + k] -= g * d[k]
                    }
                }
                for (k in 0..i) v[next + k] = 0.0
            }
            for (j in 0 until n) {
                d[j] = v[j * n + n - 1]
                v[j * n + n - 1] = 0.0
            }
            v[(n - 1) * n + n - 1] = 1.0
            e[0] = 0.0
        }

        /**
         * Implicit QL iteration with shifts on the tridiagonal matrix in [d] and [e], rotating [v] as it goes.
         * JAMA's `tql2`, without its ascending sort (see [sortDescending]).
         */
        private fun diagonalize(v: DoubleArray, n: Int, d: DoubleArray, e: DoubleArray) {
            for (i in 1 until n) e[i - 1] = e[i]
            e[n - 1] = 0.0

            var f = 0.0
            var tst1 = 0.0
            for (l in 0 until n) {
                // Find a small sub-diagonal element.
                tst1 = max(tst1, abs(d[l]) + abs(e[l]))
                var m = l
                while (m < n) {
                    if (abs(e[m]) <= EPSILON * tst1) break
                    m++
                }

                // If m == l then d[l] is already an eigenvalue; otherwise iterate.
                if (m > l) {
                    var iteration = 0
                    do {
                        iteration++
                        check(iteration <= MAX_ITERATIONS) {
                            "QL iteration did not converge for eigenvalue $l after $MAX_ITERATIONS sweeps"
                        }

                        // Compute the implicit shift.
                        var g = d[l]
                        var p = (d[l + 1] - g) / (2.0 * e[l])
                        var r = pythag(p, 1.0)
                        if (p < 0) r = -r
                        d[l] = e[l] / (p + r)
                        d[l + 1] = e[l] * (p + r)
                        val dl1 = d[l + 1]
                        var h = g - d[l]
                        for (i in l + 2 until n) d[i] -= h
                        f += h

                        // Implicit QL transformation.
                        p = d[m]
                        var c = 1.0
                        var c2 = c
                        var c3 = c
                        val el1 = e[l + 1]
                        var s = 0.0
                        var s2 = 0.0
                        for (i in m - 1 downTo l) {
                            c3 = c2
                            c2 = c
                            s2 = s
                            g = c * e[i]
                            h = c * p
                            r = pythag(p, e[i])
                            e[i + 1] = s * r
                            s = e[i] / r
                            c = p / r
                            p = c * d[i] - s * g
                            d[i + 1] = h + s * (c * g + s * d[i])

                            // Accumulate the rotation into columns i and i + 1.
                            val column = i * n
                            val next = column + n
                            for (k in 0 until n) {
                                h = v[next + k]
                                v[next + k] = s * v[column + k] + c * h
                                v[column + k] = c * v[column + k] - s * h
                            }
                        }
                        p = -s * s2 * c3 * el1 * e[l] / dl1
                        e[l] = s * p
                        d[l] = c * p
                    } while (abs(e[l]) > EPSILON * tst1)
                }
                d[l] += f
                e[l] = 0.0
            }
        }

        /** Selection sort of [d] into descending order, carrying the columns of [v] along. */
        private fun sortDescending(v: DoubleArray, n: Int, d: DoubleArray) {
            for (i in 0 until n - 1) {
                var k = i
                var p = d[i]
                for (j in i + 1 until n) {
                    if (d[j] > p) {
                        k = j
                        p = d[j]
                    }
                }
                if (k != i) {
                    d[k] = d[i]
                    d[i] = p
                    val left = i * n
                    val right = k * n
                    for (j in 0 until n) {
                        val swap = v[left + j]
                        v[left + j] = v[right + j]
                        v[right + j] = swap
                    }
                }
            }
        }

        /** Scales each column of [v] to unit length and to the documented sign convention. */
        private fun normalize(v: DoubleArray, n: Int) {
            for (c in 0 until n) {
                val column = c * n
                var sumOfSquares = 0.0
                var pivot = 0
                var largest = -1.0
                for (r in 0 until n) {
                    val value = v[column + r]
                    sumOfSquares += value * value
                    val magnitude = abs(value)
                    if (magnitude > largest) {
                        largest = magnitude
                        pivot = r
                    }
                }
                val norm = sqrt(sumOfSquares)
                if (norm == 0.0) continue
                val scale = if (v[column + pivot] < 0.0) -1.0 / norm else 1.0 / norm
                for (r in 0 until n) v[column + r] *= scale
            }
        }

        /**
         * `sqrt(a * a + b * b)` without intermediate overflow. Ported rather than taken from `kotlin.math.hypot`
         * so that every target computes the same bits: this form uses only division, multiplication, and `sqrt`,
         * which IEEE 754 requires to be exactly rounded.
         */
        private fun pythag(a: Double, b: Double): Double = if (abs(a) > abs(b)) {
            val r = b / a
            abs(a) * sqrt(1 + r * r)
        } else if (b != 0.0) {
            val r = a / b
            abs(b) * sqrt(1 + r * r)
        } else {
            0.0
        }
    }
}
