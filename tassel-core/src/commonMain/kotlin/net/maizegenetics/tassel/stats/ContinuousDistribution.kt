package net.maizegenetics.tassel.stats

import kotlin.math.abs
import kotlin.math.ln
import kotlin.math.ln1p
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/**
 * A continuous probability distribution with accurate tails: [sf] is computed directly rather than as `1 - cdf`, so
 * an upper-tail p-value is exact far below machine epsilon, and [inverseSf] solves for it directly.
 *
 * Nothing here depends on phenotype storage; a distribution is just its parameters.
 */
interface ContinuousDistribution {
    /** The probability density at [x]; zero outside the support. NaN gives NaN. */
    fun pdf(x: Double): Double

    /** The natural log of [pdf], `-Infinity` outside the support. NaN gives NaN. */
    fun logPdf(x: Double): Double

    /** `P(X <= x)`. Zero below the support and one above it. NaN gives NaN. */
    fun cdf(x: Double): Double

    /** `P(X > x)`, the upper tail, accurate to a small relative error however small it is. NaN gives NaN. */
    fun sf(x: Double): Double

    /**
     * `ln P(X > x)`, finite where [sf] underflows to zero (down to about `-1e300`) and accurate in log space
     * everywhere. `-Infinity` above the support and 0 below it. NaN gives NaN.
     */
    fun logSf(x: Double): Double

    /** `log10 P(X > x)`, [logSf] divided by `ln 10`. */
    fun log10Sf(x: Double): Double = logSf(x) / LN_10

    /**
     * The `x` with `cdf(x) = p`, to about `1e-14` relative where that is representable. `p = 0` and `p = 1` give the
     * lower and upper support bounds, which may be infinite.
     *
     * @throws IllegalArgumentException if [p] is not in `[0, 1]` (NaN included).
     */
    fun quantile(p: Double): Double

    /**
     * The `x` with `sf(x) = q`, solved on the upper tail so a tiny [q] keeps its relative accuracy. `q = 0` gives the
     * upper support bound and `q = 1` the lower one.
     *
     * @throws IllegalArgumentException if [q] is not in `[0, 1]` (NaN included).
     */
    fun inverseSf(q: Double): Double
}

internal fun requireProbability(name: String, p: Double) {
    require(p in 0.0..1.0) { "$name must be in [0, 1], got $p" }
}

internal fun requireDegreesOfFreedom(name: String, df: Double) {
    require(df.isFinite() && df > 0.0) { "$name must be a finite number greater than 0, got $df" }
}

/**
 * Solves `tail(x) = target` for a positive `x`, where [tail] is increasing in `x` for a lower tail and decreasing for
 * an upper tail ([upper]). `target` must be in `(0, 0.5]` so the tail is the small side and its relative accuracy is
 * what matters.
 *
 * The root is found on `ln tail(x) - ln target`, which is smooth even when the tail spans 300 orders of magnitude,
 * with a safeguarded Newton step (the derivative is `pdf / tail`) inside a bracket, falling back to bisection (the
 * geometric mean while the bracket spans more than a factor of four) whenever the step leaves it. It stops when a
 * step is below `1e-15` relative. A root below the smallest positive double gives 0, and above the largest gives
 * `Infinity`.
 */
internal fun solveTail(
    target: Double,
    upper: Boolean,
    start: Double,
    tail: (Double) -> Double,
    logPdf: (Double) -> Double,
): Double {
    val logTarget = ln(target)

    // Increasing in x: negative below the root, positive above it.
    fun g(x: Double): Double {
        val diff = ln(tail(x)) - logTarget
        return if (upper) -diff else diff
    }

    var x = start
    var gx = g(x)
    if (gx == 0.0) return x
    var lo: Double
    var hi: Double
    var factor = 2.0
    if (gx < 0.0) {
        lo = x
        while (true) {
            val next = min(lo * factor, Double.MAX_VALUE)
            val gn = g(next)
            if (gn == 0.0) return next
            if (gn > 0.0) {
                hi = next
                break
            }
            if (next == Double.MAX_VALUE) return Double.POSITIVE_INFINITY
            lo = next
            factor *= factor
        }
    } else {
        hi = x
        while (true) {
            val next = max(hi / factor, Double.MIN_VALUE)
            val gn = g(next)
            if (gn == 0.0) return next
            if (gn < 0.0) {
                lo = next
                break
            }
            if (next == Double.MIN_VALUE) return 0.0
            hi = next
            factor *= factor
        }
    }

    x = sqrt(lo) * sqrt(hi)
    repeat(400) {
        gx = g(x)
        if (gx == 0.0) return x
        if (gx < 0.0) lo = x else hi = x
        val slope = kotlin.math.exp(logPdf(x) - ln(tail(x)))
        var next = x - gx / slope
        if (!(next > lo && next < hi)) {
            next = if (hi > 4.0 * lo) sqrt(lo) * sqrt(hi) else lo + (hi - lo) / 2.0
        }
        if (abs(next - x) <= 1e-15 * abs(next)) return next
        x = next
    }
    return x
}

/**
 * `ln B(a, b)` for the beta function. Cancellation in `lnGamma(a) + lnGamma(b) - lnGamma(a + b)` costs a digit per
 * order of magnitude of the arguments, so when the larger is at least 20 the difference of its log gamma and that of
 * the sum is taken analytically from Stirling's series (error below `1e-15`) instead.
 */
internal fun logBeta(a: Double, b: Double): Double {
    val big = max(a, b)
    val small = min(a, b)
    val sum = a + b
    if (big < 20.0) {
        return SpecialFunctions.logGamma(a) + SpecialFunctions.logGamma(b) - SpecialFunctions.logGamma(sum)
    }
    val difference =
        -(big - 0.5) * ln1p(small / big) - small * ln(sum) + small + stirlingCorrection(big) - stirlingCorrection(sum)
    return SpecialFunctions.logGamma(small) + difference
}

/** The tail of Stirling's series for `lnGamma(z)` beyond `(z - 1/2) ln z - z + ln(2 pi) / 2`, for `z >= 20`. */
private fun stirlingCorrection(z: Double): Double {
    val inverse = 1.0 / z
    val square = inverse * inverse
    return inverse * (1.0 / 12.0 - square * (1.0 / 360.0 - square * (1.0 / 1260.0 - square / 1680.0)))
}
