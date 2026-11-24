package net.maizegenetics.tassel.stats

import kotlin.math.ln
import kotlin.math.ln1p

/**
 * The F distribution with [df1] numerator and [df2] denominator degrees of freedom, neither needing to be an integer.
 *
 * [sf] is exactly TASSEL 5's F-test p-value, [fTestPValue], which is the upper tail computed directly with no
 * `1 - cdf`. [cdf] is the matching lower-tail incomplete beta.
 *
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.Ftest
 * @see org.apache.commons.math3.distribution.FDistribution
 */
class FDistribution(val df1: Double, val df2: Double) : ContinuousDistribution {
    init {
        requireDegreesOfFreedom("df1", df1)
        requireDegreesOfFreedom("df2", df2)
    }

    private val half1 = df1 / 2.0
    private val half2 = df2 / 2.0
    private val logBeta = logBeta(half1, half2)

    override fun logPdf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x < 0.0 -> Double.NEGATIVE_INFINITY
        x == 0.0 ->
            when {
                df1 < 2.0 -> Double.POSITIVE_INFINITY
                df1 == 2.0 -> 0.0
                else -> Double.NEGATIVE_INFINITY
            }
        x == Double.POSITIVE_INFINITY -> Double.NEGATIVE_INFINITY
        else -> {
            // pdf = z^(df1 / 2) w^(df2 / 2) / (x B), with z = df1 x / (df2 + df1 x) and w = 1 - z. Each log is
            // taken from the form that does not cancel.
            val dx = df1 * x
            val logZ = if (dx <= df2) ln(dx / (df2 + dx)) else -ln1p(df2 / dx)
            val logW = if (dx <= df2) -ln1p(dx / df2) else ln(df2 / (df2 + dx))
            half1 * logZ + half2 * logW - ln(x) - logBeta
        }
    }

    override fun pdf(x: Double): Double = kotlin.math.exp(logPdf(x))

    override fun cdf(x: Double): Double {
        if (x.isNaN()) return Double.NaN
        if (x <= 0.0) return 0.0
        val dx = df1 * x
        if (dx.isInfinite()) return 1.0
        return SpecialFunctions.regularizedIncompleteBeta(half1, half2, dx / (df2 + dx))
    }

    override fun sf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x <= 0.0 -> 1.0
        else -> fTestPValue(x, df1, df2)
    }

    override fun quantile(p: Double): Double {
        requireProbability("p", p)
        return when {
            p == 0.0 -> 0.0
            p == 1.0 -> Double.POSITIVE_INFINITY
            p <= 0.5 -> solveTail(p, false, 1.0, ::cdf, ::logPdf)
            else -> solveTail(1.0 - p, true, 1.0, ::sf, ::logPdf)
        }
    }

    override fun inverseSf(q: Double): Double {
        requireProbability("q", q)
        return when {
            q == 0.0 -> Double.POSITIVE_INFINITY
            q == 1.0 -> 0.0
            q <= 0.5 -> solveTail(q, true, 1.0, ::sf, ::logPdf)
            else -> solveTail(1.0 - q, false, 1.0, ::cdf, ::logPdf)
        }
    }
}
