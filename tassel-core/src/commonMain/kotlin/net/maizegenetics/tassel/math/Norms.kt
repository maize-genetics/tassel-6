package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.math.sqrt

/**
 * `sqrt(a * a + b * b)` without intermediate overflow. Ported rather than taken from `kotlin.math.hypot` so that every
 * target computes the same bits: this form uses only division, multiplication, and `sqrt`, which IEEE 754 requires to
 * be exactly rounded.
 */
internal fun pythag(a: Double, b: Double): Double = if (abs(a) > abs(b)) {
    val r = b / a
    abs(a) * sqrt(1 + r * r)
} else if (b != 0.0) {
    val r = a / b
    abs(b) * sqrt(1 + r * r)
} else {
    0.0
}

/**
 * Euclidean norm of `array[offset until offset + length]`, scaled by the largest magnitude so that neither squaring
 * nor summing overflows or underflows. Division, multiplication, addition, and `sqrt` only, so it is bit-identical on
 * every target.
 */
internal fun norm2(array: DoubleArray, offset: Int, length: Int): Double {
    var largest = 0.0
    for (i in offset until offset + length) {
        val magnitude = abs(array[i])
        if (magnitude > largest) largest = magnitude
    }
    if (largest == 0.0 || !largest.isFinite()) return largest
    var sum = 0.0
    for (i in offset until offset + length) {
        val scaled = array[i] / largest
        sum += scaled * scaled
    }
    return largest * sqrt(sum)
}

/**
 * A running sum carried in two doubles (about 32 significant digits), after Ogita, Rump, and Oishi's `Dot2`: each
 * product is split exactly by Veltkamp's method and each addition keeps its rounding error, so the sum is accurate
 * as if computed in twice the working precision and then rounded once. Only `+`, `-`, and `*`, so every target
 * produces the same bits. One instance is reused across sums; call [reset] between them.
 */
internal class CompensatedSum {
    private var sum = 0.0
    private var error = 0.0

    /** The sum so far, rounded to one double. */
    val value: Double get() = sum + error

    fun reset() {
        sum = 0.0
        error = 0.0
    }

    fun add(a: Double) {
        val s = sum + a
        val z = s - sum
        error += (sum - (s - z)) + (a - z)
        sum = s
    }

    fun addProduct(a: Double, b: Double) {
        val p = a * b
        val ta = SPLITTER * a
        val ahi = ta - (ta - a)
        val alo = a - ahi
        val tb = SPLITTER * b
        val bhi = tb - (tb - b)
        val blo = b - bhi
        error += alo * blo - (((p - ahi * bhi) - alo * bhi) - ahi * blo)
        add(p)
    }

    private companion object {
        // 2^27 + 1: splits a double into two halves of 26 bits whose product is exact.
        const val SPLITTER = 134217729.0
    }
}
