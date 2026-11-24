package net.maizegenetics.tassel.stats

import kotlin.math.ln

/**
 * The chi-squared distribution with [df] degrees of freedom, `Gamma(df / 2, scale 2)`. [df] need not be an integer.
 *
 * [cdf] and [sf] are the regularized incomplete gamma functions `P(df / 2, x / 2)` and `Q(df / 2, x / 2)` from
 * [SpecialFunctions], so [sf] is accurate far into the tail.
 *
 * @see cern.jet.stat.Probability.chiSquare
 * @see cern.jet.stat.Probability.chiSquareComplemented
 * @see org.apache.commons.math3.distribution.ChiSquaredDistribution
 */
class ChiSquaredDistribution(val df: Double) : ContinuousDistribution {
    init {
        requireDegreesOfFreedom("df", df)
    }

    private val shape = df / 2.0
    private val logNormalizer = shape * LN_2 + SpecialFunctions.logGamma(shape)

    override fun logPdf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x < 0.0 -> Double.NEGATIVE_INFINITY
        x == 0.0 ->
            when {
                shape < 1.0 -> Double.POSITIVE_INFINITY
                shape == 1.0 -> -logNormalizer
                else -> Double.NEGATIVE_INFINITY
            }
        x == Double.POSITIVE_INFINITY -> Double.NEGATIVE_INFINITY
        else -> (shape - 1.0) * ln(x) - x / 2.0 - logNormalizer
    }

    override fun pdf(x: Double): Double = kotlin.math.exp(logPdf(x))

    override fun cdf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x <= 0.0 -> 0.0
        x == Double.POSITIVE_INFINITY -> 1.0
        else -> SpecialFunctions.regularizedGammaP(shape, x / 2.0)
    }

    override fun sf(x: Double): Double = when {
        x.isNaN() -> Double.NaN
        x <= 0.0 -> 1.0
        x == Double.POSITIVE_INFINITY -> 0.0
        else -> SpecialFunctions.regularizedGammaQ(shape, x / 2.0)
    }

    override fun quantile(p: Double): Double {
        requireProbability("p", p)
        return when {
            p == 0.0 -> 0.0
            p == 1.0 -> Double.POSITIVE_INFINITY
            p <= 0.5 -> solveTail(p, false, df, ::cdf, ::logPdf)
            else -> solveTail(1.0 - p, true, df, ::sf, ::logPdf)
        }
    }

    override fun inverseSf(q: Double): Double {
        requireProbability("q", q)
        return when {
            q == 0.0 -> Double.POSITIVE_INFINITY
            q == 1.0 -> 0.0
            q <= 0.5 -> solveTail(q, true, df, ::sf, ::logPdf)
            else -> solveTail(1.0 - q, false, df, ::cdf, ::logPdf)
        }
    }

    private companion object {
        const val LN_2 = 0.6931471805599453
    }
}
