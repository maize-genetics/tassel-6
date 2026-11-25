package net.maizegenetics.tassel.stats

import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.ln1p

/**
 * Student's t distribution with [df] degrees of freedom, which need not be an integer.
 *
 * Both tails come from the incomplete beta identity `P(T > |t|) = I_w(df / 2, 1/2) / 2` with `w = df / (df + t^2)`,
 * and the other tail follows by symmetry, so [cdf] and [sf] are each accurate in their own small tail. For `|t| < 1`,
 * where `w` approaches 1, the complementary form `1/2 - I_z(1/2, df / 2) / 2` with `z = t^2 / (df + t^2)` avoids
 * losing a small `t`. Past `w = 1e-300` the tail comes from the leading term of the beta series, which is exact there.
 *
 * The incomplete beta is [SpecialFunctions]'s, Colt's, which takes one argument and forms `1 - x` itself, so for a
 * very large [df] (beyond about 1e6) the relative accuracy falls to about 1e-8.
 *
 * @see cern.jet.stat.Probability.studentT
 * @see org.apache.commons.math3.distribution.TDistribution
 */
class StudentTDistribution(val df: Double) : ContinuousDistribution {
    init {
        requireDegreesOfFreedom("df", df)
    }

    private val halfDf = df / 2.0
    private val logNormalizer = 0.5 * ln(df) + logBeta(0.5, halfDf)

    override fun logPdf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x.isInfinite() -> Double.NEGATIVE_INFINITY
        else -> -(halfDf + 0.5) * ln1p(x * x / df) - logNormalizer
    }

    override fun pdf(x: Double): Double = exp(logPdf(x))

    /** `P(T > a)` for `a >= 0`. */
    private fun upperTail(a: Double): Double {
        if (a == 0.0) return 0.5
        if (a == Double.POSITIVE_INFINITY) return 0.0
        if (a < 1.0) {
            // The tail is at least 0.15 here, so 1/2 - half the central mass loses at most a digit and keeps small a.
            val square = a * a
            return 0.5 - 0.5 * SpecialFunctions.regularizedIncompleteBeta(0.5, halfDf, square / (df + square))
        }
        // w = df / (df + a^2), written so that a^2 does not overflow beyond a = 1e154.
        val ratio = df / a
        val w = ratio / (ratio + a)
        if (w >= 1e-300) return 0.5 * SpecialFunctions.regularizedIncompleteBeta(halfDf, 0.5, w)
        // w has underflowed. I_w(a, b) is x^a / (a B(a, b)) to far better than double precision for w below 1e-300.
        val logW = ln(ratio) - ln(ratio + a)
        return 0.5 * exp(halfDf * logW - ln(halfDf) - logBeta(halfDf, 0.5))
    }

    override fun sf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x >= 0.0 -> upperTail(x)
        else -> 1.0 - upperTail(-x)
    }

    override fun cdf(x: Double): Double = sf(-x)

    /** `ln P(T > a)` for `a >= 0`, from the same beta identity as [upperTail] with both ends of `w` exact. */
    private fun logUpperTail(a: Double): Double {
        if (a == 0.0) return -LN_2
        if (a == Double.POSITIVE_INFINITY) return Double.NEGATIVE_INFINITY
        // w = df / (df + a^2) and 1 - w = a^2 / (df + a^2), written so that a^2 does not overflow.
        val ratio = df / a
        val w = ratio / (ratio + a)
        val wc = a / (ratio + a)
        // w can underflow beyond a = 1e154 * sqrt(df), so its log comes from the ratio rather than from w.
        val logW = if (w > 0.5) ln1p(-wc) else ln(ratio) - ln(ratio + a)
        val logWc = if (wc > 0.5) ln1p(-w) else ln(a) - ln(ratio + a)
        return -LN_2 + SpecialFunctions.logIncompleteBeta(halfDf, 0.5, w, wc, logW, logWc)
    }

    override fun logSf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x >= 0.0 -> logUpperTail(x)
        else -> SpecialFunctions.log1mExp(logUpperTail(-x))
    }

    override fun quantile(p: Double): Double {
        requireProbability("p", p)
        return when {
            p == 0.0 -> Double.NEGATIVE_INFINITY
            p == 1.0 -> Double.POSITIVE_INFINITY
            p == 0.5 -> 0.0
            p < 0.5 -> -solveTail(p, true, 1.0, ::sf, ::logPdf)
            else -> solveTail(1.0 - p, true, 1.0, ::sf, ::logPdf)
        }
    }

    override fun inverseSf(q: Double): Double {
        requireProbability("q", q)
        return when {
            q == 0.0 -> Double.POSITIVE_INFINITY
            q == 1.0 -> Double.NEGATIVE_INFINITY
            q == 0.5 -> 0.0
            q < 0.5 -> solveTail(q, true, 1.0, ::sf, ::logPdf)
            else -> -solveTail(1.0 - q, true, 1.0, ::sf, ::logPdf)
        }
    }

    private companion object {
        const val LN_2 = 0.6931471805599453
    }
}
