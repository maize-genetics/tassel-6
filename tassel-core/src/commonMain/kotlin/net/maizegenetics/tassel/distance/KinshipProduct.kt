package net.maizegenetics.tassel.distance

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope

/**
 * Adds one panel's `Z Z^T` into a pairwise triangle. [panel] is row-major with [stride] columns of
 * storage and [columns] of them live. [triangle] is packed row by row over the upper triangle,
 * diagonal included: `(0,0), (0,1), ... (0,n-1), (1,1), ...`.
 *
 * The pure-Kotlin implementation is [PureKotlinSyrk]. A JVM or native fast path replaces
 * [kinshipTriangleUpdate] and leaves this contract alone (tassel-6#68, tassel-6#69).
 */
internal fun interface KinshipTriangleUpdate {
    suspend fun add(panel: DoubleArray, taxa: Int, columns: Int, stride: Int, triangle: DoubleArray)
}

/** The product the default kinship kernel calls. Off-by-default spikes assign a different one. */
internal var kinshipTriangleUpdate: KinshipTriangleUpdate = PureKotlinSyrk

/** Cache-blocked `Z Z^T` in pure Kotlin, parallel across pair-balanced row bands. */
internal object PureKotlinSyrk : KinshipTriangleUpdate {
    override suspend fun add(panel: DoubleArray, taxa: Int, columns: Int, stride: Int, triangle: DoubleArray) {
        if (taxa == 0 || columns == 0) return
        val bands = pairBands(taxa, kinshipWorkerCount())
        if (bands.size <= 1) {
            addBand(panel, taxa, columns, stride, triangle, 0 until taxa)
            return
        }
        val slices = Array(bands.size) { index -> DoubleArray(pairCount(taxa, bands[index])) }
        // Inherit the caller's dispatcher. LocalEngine puts the JVM and Native on Dispatchers.Default,
        // so the bands run on that pool; a test dispatcher and the browser's one thread stay single-threaded.
        coroutineScope {
            bands.mapIndexed { index, rows ->
                async {
                    addBand(panel, taxa, columns, stride, slices[index], rows)
                }
            }.awaitAll()
        }
        // Parent thread, band order. Scheduling does not change the sums (ADR-0020).
        var offset = 0
        for (slice in slices) {
            for (value in slice) {
                triangle[offset] += value
                offset++
            }
        }
    }
}

/**
 * Upper triangle of rows [rows], each cell the dot product of two panel rows. Four partial sums
 * keep the inner loop from waiting on a single accumulator; the panel is too small for that to
 * change a result anyone compares beyond the kinship tolerance.
 */
private fun addBand(panel: DoubleArray, taxa: Int, columns: Int, stride: Int, dest: DoubleArray, rows: IntRange) {
    var index = 0
    val limit = columns - 3
    for (first in rows) {
        val row = first * stride
        for (second in first until taxa) {
            val other = second * stride
            var k = 0
            var s0 = 0.0
            var s1 = 0.0
            var s2 = 0.0
            var s3 = 0.0
            while (k < limit) {
                s0 += panel[row + k] * panel[other + k]
                s1 += panel[row + k + 1] * panel[other + k + 1]
                s2 += panel[row + k + 2] * panel[other + k + 2]
                s3 += panel[row + k + 3] * panel[other + k + 3]
                k += 4
            }
            var sum = (s0 + s1) + (s2 + s3)
            while (k < columns) {
                sum += panel[row + k] * panel[other + k]
                k++
            }
            dest[index] += sum
            index++
        }
    }
}
