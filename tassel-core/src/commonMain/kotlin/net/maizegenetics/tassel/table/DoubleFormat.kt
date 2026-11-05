package net.maizegenetics.tassel.table

import kotlin.math.absoluteValue

/**
 * TASSEL 5's number format for table cells, reproduced exactly on every target.
 *
 * A value in `[0.001, 1000000)` is written with up to five decimals and no trailing zeros
 * (pattern `0.#####`); anything else, including **every** negative number, is written in scientific
 * notation with five significant digits (pattern `0.0###E0`). Both patterns are `Locale.ROOT` and
 * round half-even, which is the JDK `DecimalFormat` default.
 *
 * `DecimalFormat` is JVM only, so the rounding is a port, not a delegation: [format] expands the
 * value's exact binary-to-decimal expansion and rounds that. See the KDoc on [exactDigits] for the
 * two cases where the JDK's own conversion makes it round differently, one of which TASSEL's whole
 * number counts depend on.
 *
 * @see net.maizegenetics.util.DoubleFormat
 */
object DoubleFormat {
    /** What TASSEL 5 writes for a missing or undefined value. */
    const val NAN = "NaN"

    private const val FRACTION_DIGITS = 5
    private const val SIGNIFICANT_DIGITS = 5
    private const val SCIENTIFIC_BELOW = 0.001
    private const val SCIENTIFIC_FROM = 1_000_000.0

    /** [value] as TASSEL 5 writes it; a `null` is [NAN], as TASSEL 5's `format(Double)` returns. */
    fun format(value: Double?): String = if (value == null) NAN else format(value)

    /** [value] as TASSEL 5 writes it. */
    fun format(value: Double): String {
        if (value.isNaN()) return NAN
        // TASSEL 5 writes both infinities as "Infinity", without a sign.
        if (value.isInfinite()) return "Infinity"
        val negative = value < 0.0 || (value == 0.0 && 1.0 / value < 0.0)
        val magnitude = value.absoluteValue
        val text = when {
            value == 0.0 -> "0"
            value < SCIENTIFIC_BELOW || value >= SCIENTIFIC_FROM -> scientific(magnitude)
            else -> plain(magnitude)
        }
        return if (negative) "-$text" else text
    }

    /** [magnitude] under pattern `0.#####`: at most [FRACTION_DIGITS] decimals, no trailing zeros. */
    private fun plain(magnitude: Double): String {
        val expansion = exactDigits(magnitude)
        val digits = expansion.digits
        // The digits of magnitude * 10^FRACTION_DIGITS, rounded to a whole number.
        val keep = digits.length - expansion.scale + FRACTION_DIGITS
        val scaled = when {
            keep < 0 -> ""
            keep >= digits.length -> digits.padEnd(keep, '0')
            else -> round(digits, keep, expansion.tieRoundsUp)
        }
        val padded = scaled.padStart(FRACTION_DIGITS + 1, '0')
        val whole = padded.substring(0, padded.length - FRACTION_DIGITS)
        val fraction = padded.substring(padded.length - FRACTION_DIGITS).trimEnd('0')
        return if (fraction.isEmpty()) whole else "$whole.$fraction"
    }

    /**
     * [magnitude] under pattern `0.0###E0`: one digit before the point, one to four after, and an
     * exponent with no `+` and no padding.
     */
    private fun scientific(magnitude: Double): String {
        val expansion = exactDigits(magnitude)
        var exponent = expansion.digits.length - expansion.scale - 1
        var mantissa = round(expansion.digits, SIGNIFICANT_DIGITS, expansion.tieRoundsUp)
        if (mantissa.length > SIGNIFICANT_DIGITS) {
            // Rounding carried, as 9.99999e-4 does; "100000" is 1.0 one decade up.
            mantissa = mantissa.substring(0, SIGNIFICANT_DIGITS)
            exponent++
        }
        val fraction = mantissa.substring(1).trimEnd('0').ifEmpty { "0" }
        return "${mantissa[0]}.${fraction}E$exponent"
    }

    /**
     * [digits] rounded to its [keep] leading digits, half-even unless [tieRoundsUp]. The result is
     * [keep] digits long, or [keep] + 1 when rounding carried, as `"99"` kept to one digit gives
     * `"10"`.
     */
    private fun round(digits: String, keep: Int, tieRoundsUp: Boolean): String {
        if (keep >= digits.length) return digits.padEnd(keep, '0')
        val next = digits[keep]
        val roundUp = when {
            next > '5' -> true
            next < '5' -> false
            // digits is the exact value, so anything past the tie decides it.
            (keep + 1 until digits.length).any { digits[it] != '0' } -> true
            tieRoundsUp -> true
            else -> keep > 0 && (digits[keep - 1] - '0') % 2 != 0
        }
        if (!roundUp) return digits.substring(0, keep)
        val kept = digits.substring(0, keep).toCharArray()
        for (index in keep - 1 downTo 0) {
            if (kept[index] != '9') {
                kept[index]++
                return kept.concatToString()
            }
            kept[index] = '0'
        }
        return "1" + kept.concatToString()
    }

    /**
     * The exact decimal expansion of a finite, non-zero magnitude: the [digits], with no leading
     * zeros, and how many of them are fractional ([scale]).
     */
    private class Expansion(val digits: String, val scale: Int, val tieRoundsUp: Boolean)

    /**
     * Expands [magnitude] exactly. Every `Double` is `mantissa * 2^exponent`, so the decimal
     * expansion is finite and this loses nothing.
     *
     * `DecimalFormat` rounds the digits the JDK's binary-to-decimal conversion hands it, which are
     * the shortest decimal that round-trips, not the exact expansion. The two agree for every
     * normal `Double` at five digits, because the shortest decimal is then within half an ulp of
     * the value and so is also the nearest five-digit decimal. Two cases do differ:
     *
     * - A whole number under `2^53` takes the conversion's fast path, which does not report that
     *   its digits are exact, so `DecimalFormat` rounds a tie up instead of to even: it writes
     *   `1234450.0` as `1.2345E6`, not `1.2344E6`. [tieRoundsUp] reproduces that, because TASSEL 5
     *   holds counts as whole `Double`s. Whole numbers of `2^53` and up are not reproduced.
     * - Below roughly `5e-319` an ulp is wider than the fifth digit, so the shortest decimal has
     *   fewer than five digits: `DecimalFormat` writes `Double.MIN_VALUE` as `4.9E-324`, this
     *   writes its exact `4.9407E-324`. No analysis reaches those values.
     */
    private fun exactDigits(magnitude: Double): Expansion {
        val bits = magnitude.toRawBits()
        val biasedExponent = ((bits ushr 52) and 0x7FF).toInt()
        val significand = bits and 0x000F_FFFF_FFFF_FFFFL
        val mantissa = if (biasedExponent == 0) significand else significand or (1L shl 52)
        val exponent = if (biasedExponent == 0) -1074 else biasedExponent - 1075
        val whole = exponent >= 0 || mantissa.countTrailingZeroBits() >= -exponent
        val value = Decimal(mantissa)
        return if (exponent >= 0) {
            // mantissa * 2^exponent is a whole number.
            value.multiplyByPowerOf(2, exponent)
            Expansion(value.toString(), 0, whole && magnitude < TWO_TO_THE_53)
        } else {
            // mantissa / 2^-exponent is mantissa * 5^-exponent shifted -exponent places.
            value.multiplyByPowerOf(5, -exponent)
            Expansion(value.toString(), -exponent, whole && magnitude < TWO_TO_THE_53)
        }
    }

    private const val TWO_TO_THE_53 = 9007199254740992.0
}

/** A non-negative whole number, in base 10^9 limbs, little endian: just enough to expand a `Double`. */
private class Decimal(value: Long) {
    private var limbs = IntArray(4)
    private var size = 0

    init {
        var rest = value
        do {
            append((rest % BASE).toInt())
            rest /= BASE
        } while (rest > 0)
    }

    /** Multiplies by [base] (2 or 5) raised to [power]. */
    fun multiplyByPowerOf(base: Int, power: Int) {
        val chunk = if (base == 2) 30 else 13 // the largest power of each that fits an Int
        var remaining = power
        while (remaining > 0) {
            val step = minOf(chunk, remaining)
            var factor = 1
            repeat(step) { factor *= base }
            multiply(factor)
            remaining -= step
        }
    }

    /** The digits, with no leading zeros. */
    override fun toString(): String {
        val text = StringBuilder(size * LIMB_DIGITS)
        text.append(limbs[size - 1])
        for (index in size - 2 downTo 0) {
            val limb = limbs[index].toString()
            repeat(LIMB_DIGITS - limb.length) { text.append('0') }
            text.append(limb)
        }
        return text.toString()
    }

    private fun multiply(factor: Int) {
        var carry = 0L
        for (index in 0 until size) {
            // A limb is under 10^9 and a factor under 2^31, so this stays well inside a Long.
            val product = limbs[index].toLong() * factor + carry
            limbs[index] = (product % BASE).toInt()
            carry = product / BASE
        }
        while (carry > 0) {
            append((carry % BASE).toInt())
            carry /= BASE
        }
    }

    private fun append(limb: Int) {
        if (size == limbs.size) limbs = limbs.copyOf(limbs.size * 2)
        limbs[size++] = limb
    }

    private companion object {
        const val BASE = 1_000_000_000L
        const val LIMB_DIGITS = 9
    }
}
