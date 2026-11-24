package net.maizegenetics.tassel.stats

import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.floor
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sin

/**
 * Special functions behind every TASSEL 5 GLM and MLM p-value: a pure-Kotlin port of Colt 1.2.0's
 * `cern.jet.stat.Gamma`, which is itself a port of the Cephes `gamma.c`, `igam.c`, and `incbet.c`.
 *
 * The control flow, constants, series lengths, and convergence tests follow Colt exactly, including the swap of
 * `a` and `b` in [regularizedIncompleteBeta] with its `1 - MACHEP` cap and the switch to log space when a power
 * would overflow, so values agree with TASSEL 5 to rounding. Only `kotlin.math` is used.
 *
 * Copyright (c) 1999 CERN - European Organization for Nuclear Research. Permission to use, copy, modify, and
 * distribute this software and its documentation for any purpose is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that copyright notice and this permission notice appear
 * in supporting documentation. CERN makes no representations about the suitability of this software for any
 * purpose. It is provided "as is" without expressed or implied warranty. The Cephes routines are by Stephen L.
 * Moshier.
 *
 * @see cern.jet.stat.Gamma
 */
object SpecialFunctions {
    private const val MACHEP = 1.11022302462515654042E-16
    private const val MAXLOG = 7.09782712893383996732E2
    private const val MINLOG = -7.451332191019412076235E2
    private const val MAXGAM = 171.624376956302725
    private const val BIG = 4.503599627370496e15
    private const val BIGINV = 2.22044604925031308085e-16
    private const val LOGPI = 1.14472988584940017414

    private val GAMMA_P =
        doubleArrayOf(
            1.60119522476751861407E-4,
            1.19135147006586384913E-3,
            1.04213797561761569935E-2,
            4.76367800457137231464E-2,
            2.07448227648435975150E-1,
            4.94214826801497100753E-1,
            9.99999999999999996796E-1,
        )
    private val GAMMA_Q =
        doubleArrayOf(
            -2.31581873324120129819E-5,
            5.39605580493303397842E-4,
            -4.45641913851797240494E-3,
            1.18139785222060435552E-2,
            3.58236398605498653373E-2,
            -2.34591795718243348568E-1,
            7.14304917030273074085E-2,
            1.00000000000000000320E0,
        )
    private val STIR =
        doubleArrayOf(
            7.87311395793093628397E-4,
            -2.29549961613378126380E-4,
            -2.68132617805781232825E-3,
            3.47222221605458667310E-3,
            8.33333333333482257126E-2,
        )
    private const val MAXSTIR = 143.01608
    private const val SQTPI = 2.50662827463100050242E0
    private val LGAM_A =
        doubleArrayOf(
            8.11614167470508450300E-4,
            -5.95061904284301438324E-4,
            7.93650340457716943945E-4,
            -2.77777777730099687205E-3,
            8.33333333333331927722E-2,
        )
    private val LGAM_B =
        doubleArrayOf(
            -1.37825152569120859100E3,
            -3.88016315134637840924E4,
            -3.31612992738871184744E5,
            -1.16237097492762307383E6,
            -1.72173700820839662146E6,
            -8.53555664245765465627E5,
        )
    private val LGAM_C =
        doubleArrayOf(
            -3.51815701436523470549E2,
            -1.70642106651881159223E4,
            -2.20528590553854454839E5,
            -1.13933444367982507207E6,
            -2.53252307177582951285E6,
            -2.01889141433532773231E6,
        )

    /**
     * The gamma function. Throws [ArithmeticException] where Colt does: at zero and the negative integers, and
     * for a large negative argument whose sine term is zero. Overflows to infinity above about 171.62.
     *
     * Colt's reflection branch for `x < -33` returns `-z` whatever the parity of `floor(-x)`, so the sign of
     * gamma at large negative non-integers is Colt's, not the true sign.
     */
    fun gamma(x: Double): Double {
        var xv = x
        var q = abs(xv)
        if (q > 33.0) {
            if (xv < 0.0) {
                var p = floor(q)
                if (p == q) throw ArithmeticException("gamma: overflow")
                var z = q - p
                if (z > 0.5) {
                    p += 1.0
                    z = q - p
                }
                z = q * sin(PI * z)
                if (z == 0.0) throw ArithmeticException("gamma: overflow")
                z = abs(z)
                z = PI / (z * stirlingFormula(q))
                return -z
            }
            return stirlingFormula(xv)
        }
        var z = 1.0
        while (xv >= 3.0) {
            xv -= 1.0
            z *= xv
        }
        while (xv < 0.0) {
            if (xv > -1.0E-9) return z / ((1.0 + 0.5772156649015329 * xv) * xv)
            z /= xv
            xv += 1.0
        }
        while (xv < 2.0) {
            if (xv == 0.0) throw ArithmeticException("gamma: singular")
            if (xv < 1.0e-9) return z / ((1.0 + 0.5772156649015329 * xv) * xv)
            z /= xv
            xv += 1.0
        }
        if (xv == 2.0 || xv == 3.0) return z
        xv -= 2.0
        val p = polevl(xv, GAMMA_P, 6)
        q = polevl(xv, GAMMA_Q, 7)
        return z * p / q
    }

    /** The natural log of the absolute value of the gamma function. Throws [ArithmeticException] where Colt does. */
    fun logGamma(x: Double): Double {
        var xv = x
        if (xv < -34.0) {
            val q = -xv
            val w = logGamma(q)
            var p = floor(q)
            if (p == q) throw ArithmeticException("lgam: Overflow")
            var z = q - p
            if (z > 0.5) {
                p += 1.0
                z = p - q
            }
            z = q * sin(PI * z)
            if (z == 0.0) throw ArithmeticException("lgamma: Overflow")
            return LOGPI - ln(z) - w
        }
        if (xv < 13.0) {
            var z = 1.0
            while (xv >= 3.0) {
                xv -= 1.0
                z *= xv
            }
            while (xv < 2.0) {
                if (xv == 0.0) throw ArithmeticException("lgamma: Overflow")
                z /= xv
                xv += 1.0
            }
            if (z < 0.0) z = -z
            if (xv == 2.0) return ln(z)
            xv -= 2.0
            val p = xv * polevl(xv, LGAM_B, 5) / p1evl(xv, LGAM_C, 6)
            return ln(z) + p
        }
        if (xv > 2.556348e305) throw ArithmeticException("lgamma: Overflow")
        var q = (xv - 0.5) * ln(xv) - xv + 0.91893853320467274178
        if (xv > 1.0e8) return q
        val p = 1.0 / (xv * xv)
        q +=
            if (xv >= 1000.0) {
                ((7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p + 0.0833333333333333333333) / xv
            } else {
                polevl(p, LGAM_A, 4) / xv
            }
        return q
    }

    /**
     * The regularized incomplete beta function `I_x(a, b)`, the lower tail, as TASSEL 5 computes it.
     *
     * @throws IllegalArgumentException for `a <= 0`, `b <= 0`, or `x` outside `[0, 1]` (Colt's "incbet: Domain
     *   error!"). `x = 0` gives 0 and `x = 1` gives 1. A NaN argument gives NaN.
     */
    fun regularizedIncompleteBeta(a: Double, b: Double, x: Double): Double {
        if (a.isNaN() || b.isNaN() || x.isNaN()) return Double.NaN
        require(a > 0.0 && b > 0.0) { "incbet: Domain error! a=$a, b=$b" }
        if (x <= 0.0 || x >= 1.0) {
            if (x == 0.0) return 0.0
            if (x == 1.0) return 1.0
            throw IllegalArgumentException("incbet: Domain error! x=$x")
        }
        if (b * x <= 1.0 && x <= 0.95) return betaPowerSeries(a, b, x)

        var w = 1.0 - x
        val flag = x > a / (a + b)
        val aa: Double
        val bb: Double
        val xc: Double
        val xx: Double
        if (flag) {
            aa = b
            bb = a
            xc = x
            xx = w
        } else {
            aa = a
            bb = b
            xc = w
            xx = x
        }
        if (flag && bb * xx <= 1.0 && xx <= 0.95) return capped(betaPowerSeries(aa, bb, xx))

        var y = xx * (aa + bb - 2.0) - (aa - 1.0)
        w = if (y < 0.0) betaFraction1(aa, bb, xx) else betaFraction2(aa, bb, xx) / xc

        y = aa * ln(xx)
        var t = bb * ln(xc)
        if (aa + bb < MAXGAM && abs(y) < MAXLOG && abs(t) < MAXLOG) {
            t = xc.pow(bb)
            t *= xx.pow(aa)
            t /= aa
            t *= w
            t *= gamma(aa + bb) / (gamma(aa) * gamma(bb))
            return if (flag) capped(t) else t
        }
        y += t + logGamma(aa + bb) - logGamma(aa) - logGamma(bb)
        y += ln(w / aa)
        t = if (y < MINLOG) 0.0 else exp(y)
        return if (flag) capped(t) else t
    }

    /** `1 - t` for the swapped branch, held just below 1 when `t` is at most machine epsilon. */
    private fun capped(t: Double): Double = if (t <= MACHEP) 1.0 - MACHEP else 1.0 - t

    /**
     * The regularized lower incomplete gamma function `P(a, x)`. As in Colt, `x <= 0` or `a <= 0` gives 0 rather
     * than an error. A NaN argument gives NaN.
     */
    fun regularizedGammaP(a: Double, x: Double): Double {
        if (a.isNaN() || x.isNaN()) return Double.NaN
        if (x <= 0.0 || a <= 0.0) return 0.0
        if (x > 1.0 && x > a) return 1.0 - regularizedGammaQ(a, x)
        var ax = a * ln(x) - x - logGamma(a)
        if (ax < -MAXLOG) return 0.0
        ax = exp(ax)
        var r = a
        var c = 1.0
        var ans = 1.0
        do {
            r += 1.0
            c *= x / r
            ans += c
        } while (c / ans > MACHEP)
        return ans * ax / a
    }

    /**
     * The regularized upper incomplete gamma function `Q(a, x)`, the complement of [regularizedGammaP]. As in
     * Colt, `x <= 0` or `a <= 0` gives 1. A NaN argument gives NaN.
     */
    fun regularizedGammaQ(a: Double, x: Double): Double {
        if (a.isNaN() || x.isNaN()) return Double.NaN
        if (x <= 0.0 || a <= 0.0) return 1.0
        if (x < 1.0 || x < a) return 1.0 - regularizedGammaP(a, x)
        var ax = a * ln(x) - x - logGamma(a)
        if (ax < -MAXLOG) return 0.0
        ax = exp(ax)

        var y = 1.0 - a
        var z = x + y + 1.0
        var c = 0.0
        var pkm2 = 1.0
        var qkm2 = x
        var pkm1 = x + 1.0
        var qkm1 = z * x
        var ans = pkm1 / qkm1
        var t: Double
        do {
            c += 1.0
            y += 1.0
            z += 2.0
            val yc = y * c
            val pk = pkm1 * z - pkm2 * yc
            val qk = qkm1 * z - qkm2 * yc
            if (qk != 0.0) {
                val r = pk / qk
                t = abs((ans - r) / r)
                ans = r
            } else {
                t = 1.0
            }
            pkm2 = pkm1
            pkm1 = pk
            qkm2 = qkm1
            qkm1 = qk
            if (abs(pk) > BIG) {
                pkm2 *= BIGINV
                pkm1 *= BIGINV
                qkm2 *= BIGINV
                qkm1 *= BIGINV
            }
        } while (t > MACHEP)
        return ans * ax
    }

    /** Continued fraction expansion #1 for the incomplete beta integral. */
    private fun betaFraction1(a: Double, b: Double, x: Double): Double {
        var k1 = a
        var k2 = a + b
        var k3 = a
        var k4 = a + 1.0
        var k5 = 1.0
        var k6 = b - 1.0
        var k7 = k4
        var k8 = a + 2.0
        var pkm2 = 0.0
        var qkm2 = 1.0
        var pkm1 = 1.0
        var qkm1 = 1.0
        var ans = 1.0
        var r = 1.0
        val thresh = 3.0 * MACHEP
        for (n in 0 until 300) {
            var xk = -(x * k1 * k2) / (k3 * k4)
            var pk = pkm1 + pkm2 * xk
            var qk = qkm1 + qkm2 * xk
            pkm2 = pkm1
            pkm1 = pk
            qkm2 = qkm1
            qkm1 = qk

            xk = (x * k5 * k6) / (k7 * k8)
            pk = pkm1 + pkm2 * xk
            qk = qkm1 + qkm2 * xk
            pkm2 = pkm1
            pkm1 = pk
            qkm2 = qkm1
            qkm1 = qk

            if (qk != 0.0) r = pk / qk
            val t: Double
            if (r != 0.0) {
                t = abs((ans - r) / r)
                ans = r
            } else {
                t = 1.0
            }
            if (t < thresh) return ans

            k1 += 1.0
            k2 += 1.0
            k3 += 2.0
            k4 += 2.0
            k5 += 1.0
            k6 -= 1.0
            k7 += 2.0
            k8 += 2.0

            if (abs(qk) + abs(pk) > BIG) {
                pkm2 *= BIGINV
                pkm1 *= BIGINV
                qkm2 *= BIGINV
                qkm1 *= BIGINV
            }
            if (abs(qk) < BIGINV || abs(pk) < BIGINV) {
                pkm2 *= BIG
                pkm1 *= BIG
                qkm2 *= BIG
                qkm1 *= BIG
            }
        }
        return ans
    }

    /** Continued fraction expansion #2 for the incomplete beta integral. */
    private fun betaFraction2(a: Double, b: Double, x: Double): Double {
        var k1 = a
        var k2 = b - 1.0
        var k3 = a
        var k4 = a + 1.0
        var k5 = 1.0
        var k6 = a + b
        var k7 = a + 1.0
        var k8 = a + 2.0
        var pkm2 = 0.0
        var qkm2 = 1.0
        var pkm1 = 1.0
        var qkm1 = 1.0
        val z = x / (1.0 - x)
        var ans = 1.0
        var r = 1.0
        val thresh = 3.0 * MACHEP
        for (n in 0 until 300) {
            var xk = -(z * k1 * k2) / (k3 * k4)
            var pk = pkm1 + pkm2 * xk
            var qk = qkm1 + qkm2 * xk
            pkm2 = pkm1
            pkm1 = pk
            qkm2 = qkm1
            qkm1 = qk

            xk = (z * k5 * k6) / (k7 * k8)
            pk = pkm1 + pkm2 * xk
            qk = qkm1 + qkm2 * xk
            pkm2 = pkm1
            pkm1 = pk
            qkm2 = qkm1
            qkm1 = qk

            if (qk != 0.0) r = pk / qk
            val t: Double
            if (r != 0.0) {
                t = abs((ans - r) / r)
                ans = r
            } else {
                t = 1.0
            }
            if (t < thresh) return ans

            k1 += 1.0
            k2 -= 1.0
            k3 += 2.0
            k4 += 2.0
            k5 += 1.0
            k6 += 1.0
            k7 += 2.0
            k8 += 2.0

            if (abs(qk) + abs(pk) > BIG) {
                pkm2 *= BIGINV
                pkm1 *= BIGINV
                qkm2 *= BIGINV
                qkm1 *= BIGINV
            }
            if (abs(qk) < BIGINV || abs(pk) < BIGINV) {
                pkm2 *= BIG
                pkm1 *= BIG
                qkm2 *= BIG
                qkm1 *= BIG
            }
        }
        return ans
    }

    /** Power series for the incomplete beta integral, used when `b * x <= 1` and `x <= 0.95`. */
    private fun betaPowerSeries(a: Double, b: Double, x: Double): Double {
        val ai = 1.0 / a
        var u = (1.0 - b) * x
        var v = u / (a + 1.0)
        val t1 = v
        var t = u
        var n = 2.0
        var s = 0.0
        val z = MACHEP * ai
        while (abs(v) > z) {
            u = (n - b) * x / n
            t *= u
            v = t / (a + n)
            s += v
            n += 1.0
        }
        s += t1
        s += ai

        u = a * ln(x)
        if (a + b < MAXGAM && abs(u) < MAXLOG) {
            t = gamma(a + b) / (gamma(a) * gamma(b))
            s = s * t * x.pow(a)
        } else {
            t = logGamma(a + b) - logGamma(a) - logGamma(b) + u + ln(s)
            s = if (t < MINLOG) 0.0 else exp(t)
        }
        return s
    }

    /** Gamma function for large positive arguments by Stirling's formula. */
    private fun stirlingFormula(x: Double): Double {
        var w = 1.0 / x
        w = 1.0 + w * polevl(w, STIR, 4)
        var y = exp(x)
        if (x > MAXSTIR) {
            val v = x.pow(0.5 * x - 0.25)
            y = v * (v / y)
        } else {
            y = x.pow(x - 0.5) / y
        }
        return SQTPI * y * w
    }

    /** Evaluates the degree-[n] polynomial with coefficients [coef], highest power first. */
    private fun polevl(x: Double, coef: DoubleArray, n: Int): Double {
        var ans = coef[0]
        for (i in 1..n) ans = ans * x + coef[i]
        return ans
    }

    /** As [polevl] for a polynomial whose leading coefficient is 1 and not stored; [coef] holds the other [n]. */
    private fun p1evl(x: Double, coef: DoubleArray, n: Int): Double {
        var ans = x + coef[0]
        for (i in 1 until n) ans = ans * x + coef[i]
        return ans
    }
}

/**
 * The upper-tail p-value of an F statistic with [df1] numerator and [df2] denominator degrees of freedom, as
 * TASSEL 5 computes it for every GLM and MLM term: `incompleteBeta(df2 / 2, df1 / 2, df2 / (df2 + df1 * f))`,
 * which is the upper tail directly with no `1 - cdf`. Degrees of freedom need not be integers (Wilks and Rao
 * F approximations).
 *
 * `f = 0` gives 1 and `f = Infinity` gives 0. A NaN input gives NaN.
 *
 * @throws IllegalArgumentException for `f < 0` (the beta argument exceeds 1), or `df1 <= 0` or `df2 <= 0`.
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.Ftest
 */
fun fTestPValue(f: Double, df1: Double, df2: Double): Double =
    SpecialFunctions.regularizedIncompleteBeta(df2 / 2.0, df1 / 2.0, df2 / (df2 + df1 * f))
