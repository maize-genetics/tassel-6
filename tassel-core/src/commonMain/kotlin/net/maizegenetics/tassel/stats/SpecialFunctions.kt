package net.maizegenetics.tassel.stats

import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.expm1
import kotlin.math.floor
import kotlin.math.ln
import kotlin.math.ln1p
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

    // ---- Log-space tails (ADR-0025) -------------------------------------------------------------------------------

    private const val LOG_EPS = 1.1102230246251565e-16
    private const val FPMIN = 1.0e-300
    private const val LOG_MAX_ITERATIONS = 200_000

    /**
     * `ln I_x(a, b)`, the natural log of [regularizedIncompleteBeta], evaluated in log space so it stays finite where
     * `I_x(a, b)` underflows (down to a log of about `-1e300`) and accurate where `I_x(a, b)` is within rounding of 1.
     *
     * The prefactor `x^a (1 - x)^b / (a B(a, b))` is formed as a sum of logs and the Lentz continued fraction, which is
     * `O(1)`, is added as a log. The fraction converges on the side of `x < (a + 1) / (a + b + 2)`; the other side is
     * `ln(1 - I_{1-x}(b, a))`. Unlike [regularizedIncompleteBeta] this is not a branch-for-branch Cephes port: the
     * Cephes branches that stay in `p` are only taken to agree with it where `p` is representable.
     *
     * `1 - x` is formed as `1.0 - x`, which loses relative accuracy when `x` is near 1; callers that have it exactly
     * should use the four-argument internal form.
     *
     * @throws IllegalArgumentException for `a <= 0`, `b <= 0`, or `x` outside `[0, 1]`. `x = 0` gives `-Infinity`
     *   and `x = 1` gives 0. A NaN argument gives NaN.
     */
    fun logRegularizedIncompleteBeta(a: Double, b: Double, x: Double): Double {
        if (a.isNaN() || b.isNaN() || x.isNaN()) return Double.NaN
        require(a > 0.0 && b > 0.0) { "incbet: Domain error! a=$a, b=$b" }
        require(x in 0.0..1.0) { "incbet: Domain error! x=$x" }
        return logIncompleteBeta(a, b, x, 1.0 - x)
    }

    /**
     * `ln I_x(a, b)` for the pair `x`, `xc = 1 - x`, each given to full relative accuracy, with their logs when `x` or
     * `xc` is below the smallest double and so cannot carry its own log. `logX = -Infinity` means `x = 0` and
     * `logXc = -Infinity` means `x = 1`.
     */
    internal fun logIncompleteBeta(
        a: Double,
        b: Double,
        x: Double,
        xc: Double,
        logX: Double = naturalLog(x, xc),
        logXc: Double = naturalLog(xc, x),
    ): Double {
        if (logX == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY
        if (logXc == Double.NEGATIVE_INFINITY) return 0.0
        if (x < (a + 1.0) / (a + b + 2.0)) return logBetaDirect(a, b, x, xc, logX, logXc)
        return log1mExp(logBetaDirect(b, a, xc, x, logXc, logX))
    }

    /** `ln v` for `v` with exact complement [vc], taking the `ln1p` form where `v` is near 1. */
    internal fun naturalLog(v: Double, vc: Double): Double = if (v > 0.5) ln1p(-vc) else ln(v)

    /** `ln I_x(a, b)` by the continued fraction, for `x` on the side where it converges quickly. */
    private fun logBetaDirect(a: Double, b: Double, x: Double, xc: Double, logX: Double, logXc: Double): Double =
        a * logX + b * logXc - logBeta(a, b) - ln(a) + ln(betaLentz(a, b, x, xc))

    /**
     * The incomplete beta continued fraction by the modified Lentz method (Numerical Recipes `betacf`).
     *
     * For `x > 1/2` every `1 + (partial term)` that sits near `1 - x` is formed from the exact complement [xc]
     * instead (see [betaLentzNearOne]): with large `a` the plain form subtracts two numbers within `1 - x` of each other
     * at each step and loses up to `log10(1 / (1 - x))` digits.
     */
    private fun betaLentz(a: Double, b: Double, x: Double, xc: Double): Double {
        if (x > 0.5) return betaLentzNearOne(a, b, x, xc)
        val qab = a + b
        val qap = a + 1.0
        val qam = a - 1.0
        var c = 1.0
        var d = 1.0 - qab * x / qap
        if (abs(d) < FPMIN) d = FPMIN
        d = 1.0 / d
        var h = d
        for (m in 1..LOG_MAX_ITERATIONS) {
            val mm = m.toDouble()
            val m2 = 2.0 * mm
            var aa = mm * (b - mm) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            if (abs(d) < FPMIN) d = FPMIN
            c = 1.0 + aa / c
            if (abs(c) < FPMIN) c = FPMIN
            d = 1.0 / d
            h *= d * c
            aa = -(a + mm) * (qab + mm) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            if (abs(d) < FPMIN) d = FPMIN
            c = 1.0 + aa / c
            if (abs(c) < FPMIN) c = FPMIN
            d = 1.0 / d
            val delta = d * c
            h *= delta
            if (abs(delta - 1.0) < LOG_EPS) break
        }
        return h
    }

    /**
     * [betaLentz] for `x > 1/2`, the same recurrence with the odd step's `1 + aa * d` and `1 + aa / c`, where
     * `aa = -x r` and `r = (a + m)(a + b + m) / ((a + 2m)(a + 1 + 2m))` is within `O(1/a)` of 1, rewritten as
     * `(1 - 1/p) + (1/p)((1 - r) + r (1 - x))` with `p` the previous `d` or `c`. Here `1 - r` is a closed form with no
     * subtraction and `1 - x` is [xc], so no step cancels to a small difference of nearly equal numbers.
     */
    private fun betaLentzNearOne(a: Double, b: Double, x: Double, xc: Double): Double {
        val qab = a + b
        val qap = a + 1.0
        val qam = a - 1.0
        var c = 1.0
        var d = ((1.0 - b) + qab * xc) / qap
        if (abs(d) < FPMIN) d = FPMIN
        d = 1.0 / d
        var h = d
        for (m in 1..LOG_MAX_ITERATIONS) {
            val mm = m.toDouble()
            val m2 = 2.0 * mm
            // Even step: aa is small, so D = 1 + aa d and C = 1 + aa / c are taken as they are.
            val even = mm * (b - mm) * x / ((qam + m2) * (a + m2))
            val dEven = even * d
            val cEven = even / c
            var dd = 1.0 + dEven
            if (abs(dd) < FPMIN) dd = FPMIN
            var cc = 1.0 + cEven
            if (abs(cc) < FPMIN) cc = FPMIN
            d = 1.0 / dd
            c = cc
            h *= d * c
            // Odd step: aa = -x r.
            val denominator = (a + m2) * (qap + m2)
            val r = (a + mm) * (qab + mm) / denominator
            val oneMinusR = (a * (1.0 - b) + mm * (2.0 * a + 2.0 - b) + 3.0 * mm * mm) / denominator
            val oneMinusXr = oneMinusR + r * xc
            // 1 - d = (D - 1) / D = dEven * d, and likewise 1 - 1/c = cEven / c.
            dd = d * (dEven + oneMinusXr)
            if (abs(dd) < FPMIN) dd = FPMIN
            cc = (cEven + oneMinusXr) / c
            if (abs(cc) < FPMIN) cc = FPMIN
            d = 1.0 / dd
            c = cc
            val delta = d * c
            h *= delta
            if (abs(delta - 1.0) < LOG_EPS) break
        }
        return h
    }

    /**
     * `ln Q(a, x)`, the natural log of [regularizedGammaQ], finite where `Q` underflows and accurate where `Q` is
     * within rounding of 1. As in [regularizedGammaQ], `x <= 0` or `a <= 0` gives 0 (`Q = 1`), and `x = Infinity`
     * gives `-Infinity`. A NaN argument gives NaN.
     *
     * For `x < a + 1` it is `ln(1 - P)` from the power series; beyond that the Lentz continued fraction for `Q` is
     * added to the log prefactor `a ln x - x - lnGamma(a)`.
     */
    fun logRegularizedGammaQ(a: Double, x: Double): Double {
        if (a.isNaN() || x.isNaN()) return Double.NaN
        if (x <= 0.0 || a <= 0.0) return 0.0
        if (x == Double.POSITIVE_INFINITY) return Double.NEGATIVE_INFINITY
        val logPrefactor = logGammaPrefactor(a, x)
        if (x < a + 1.0) {
            var term = 1.0 / a
            var sum = term
            var ap = a
            for (n in 1..LOG_MAX_ITERATIONS) {
                ap += 1.0
                term *= x / ap
                sum += term
                if (abs(term) < abs(sum) * LOG_EPS) break
            }
            return log1mExp(logPrefactor + ln(sum))
        }
        var b = x + 1.0 - a
        var c = 1.0 / FPMIN
        var d = 1.0 / b
        var h = d
        for (i in 1..LOG_MAX_ITERATIONS) {
            val an = -i.toDouble() * (i.toDouble() - a)
            b += 2.0
            d = an * d + b
            if (abs(d) < FPMIN) d = FPMIN
            c = b + an / c
            if (abs(c) < FPMIN) c = FPMIN
            d = 1.0 / d
            val delta = d * c
            h *= delta
            if (abs(delta - 1.0) < LOG_EPS) break
        }
        return logPrefactor + ln(h)
    }

    /**
     * `ln(x^a e^-x / Gamma(a))`. Near `x = a` with large `a` the direct form cancels terms of size `a`, so from
     * `a = 15` it is taken from Loader's saddle-point form, `ln a - stirlerr(a) - bd0(a, x) - ln(2 pi a) / 2`, where
     * `stirlerr` is the tail of Stirling's series and `bd0` is `a ln(a / x) + x - a` summed as a series near `x = a`.
     */
    private fun logGammaPrefactor(a: Double, x: Double): Double {
        if (a < 15.0) return a * ln(x) - x - logGamma(a)
        return ln(a) - stirlingError(a) - deviance(a, x) - 0.5 * ln(2.0 * PI * a)
    }

    /** `lnGamma(n + 1) - ((n + 1/2) ln n - n + ln(2 pi) / 2)` for `n >= 15`, by its asymptotic series. */
    internal fun stirlingError(n: Double): Double {
        val inverse = 1.0 / n
        val square = inverse * inverse
        return inverse *
            (1.0 / 12.0 - square * (1.0 / 360.0 - square * (1.0 / 1260.0 - square * (1.0 / 1680.0 - square / 1188.0))))
    }

    /** `x ln(x / mean) + mean - x` without cancellation when `x` is close to [mean] (Loader's `bd0`). */
    internal fun deviance(x: Double, mean: Double): Double {
        if (abs(x - mean) < 0.1 * (x + mean)) {
            val v = (x - mean) / (x + mean)
            var sum = (x - mean) * v
            var term = 2.0 * x * v
            val vSquare = v * v
            var j = 1
            while (j < 1000) {
                term *= vSquare
                val next = sum + term / (2 * j + 1)
                if (next == sum) return next
                sum = next
                j++
            }
            return sum
        }
        return x * ln(x / mean) + mean - x
    }

    /** `ln(1 - exp(l))` for `l <= 0`, without cancellation at either end. */
    internal fun log1mExp(l: Double): Double = if (l > -0.6931471805599453) ln(-expm1(l)) else ln1p(-exp(l)) + 0.0

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

/**
 * `ln` of [fTestPValue], evaluated in log space so it stays finite where the p-value underflows to 0. The beta
 * argument is formed as the exact pair `q / (1 + q)` and `1 / (1 + q)` with `q = df2 / (df1 f)`, so neither end loses
 * relative accuracy. Same domain rules as [fTestPValue]: `f = 0` gives 0, `f = Infinity` gives `-Infinity`, NaN gives
 * NaN, and `f < 0` or a non-positive degree of freedom throws.
 */
internal fun fTestLogPValue(f: Double, df1: Double, df2: Double): Double {
    if (f.isNaN() || df1.isNaN() || df2.isNaN()) return Double.NaN
    require(df1 > 0.0 && df2 > 0.0) { "incbet: Domain error! a=${df2 / 2.0}, b=${df1 / 2.0}" }
    require(f >= 0.0) { "incbet: Domain error! x=${df2 / (df2 + df1 * f)}" }
    if (f == 0.0) return 0.0
    if (f == Double.POSITIVE_INFINITY) return Double.NEGATIVE_INFINITY
    val dx = df1 * f
    val q = if (dx.isInfinite() || dx == 0.0) 0.0 else df2 / dx
    if (q < 1e-300 || q > 1e300) {
        // Either end of the beta argument has left the normal range: take both from logQ = ln q.
        val logQ = ln(df2) - ln(df1) - ln(f)
        val logOnePlusQ = if (logQ > 0.0) logQ + ln1p(exp(-logQ)) else ln1p(exp(logQ))
        return SpecialFunctions.logIncompleteBeta(
            df2 / 2.0,
            df1 / 2.0,
            exp(logQ - logOnePlusQ),
            exp(-logOnePlusQ),
            logQ - logOnePlusQ,
            -logOnePlusQ,
        )
    }
    return SpecialFunctions.logIncompleteBeta(df2 / 2.0, df1 / 2.0, q / (1.0 + q), 1.0 / (1.0 + q))
}

/**
 * `log10` of [fTestPValue], finite where the p-value underflows to 0 (TASSEL 5's `0.0`), and equal to
 * `log10(fTestPValue(f, df1, df2))` to about `1e-12` wherever the p-value is at least `1e-290`. Same domain rules as
 * [fTestPValue].
 *
 * @throws IllegalArgumentException for `f < 0`, or `df1 <= 0` or `df2 <= 0`.
 * @see fTestPValue
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.Ftest
 */
fun fTestLog10PValue(f: Double, df1: Double, df2: Double): Double = fTestLogPValue(f, df1, df2) / LN_10

internal val LN_10: Double = ln(10.0)
