package net.maizegenetics.tassel.phenotype

import net.maizegenetics.tassel.table.Decimal
import kotlin.math.max
import kotlin.math.min

/**
 * Java's `Float.parseFloat` for decimal text, reproduced bit for bit on every target.
 *
 * TASSEL 5 parses every phenotype cell with `Float.parseFloat`, so a TASSEL 6 table holds the same
 * floats only if it parses the same way. Kotlin's `String.toFloat()` is that function on the JVM
 * alone. Elsewhere it parses a `Double` and narrows it, which rounds twice, and it does not accept
 * Java's `f` and `d` suffixes. This parses the `Double` too, but when that `Double` lands exactly
 * halfway between two floats, where the second rounding can go the wrong way, it settles the tie
 * against the exact value of the text.
 *
 * The grammar is Java's for decimal numbers: an optional sign, digits with an optional point (at
 * least one digit on either side), an optional exponent, and an optional `f`, `F`, `d`, or `D`;
 * or `NaN` or `Infinity` after the sign. Surrounding characters up to and including a space are
 * trimmed, as `String.trim` trims them. Hexadecimal literals (`0x1p3`) are not accepted.
 *
 * @see java.lang.Float.parseFloat
 */
internal object JavaFloat {
    /** [text] as `Float.parseFloat` reads it, or `null` where that would throw. */
    fun parse(text: String): Float? {
        val token = text.trim { it <= ' ' }
        if (token.isEmpty()) return null
        val negative = token[0] == '-'
        val body = if (token[0] == '-' || token[0] == '+') token.substring(1) else token
        val magnitude = when (body) {
            "NaN" -> return Float.NaN
            "Infinity" -> Float.POSITIVE_INFINITY
            else -> parseMagnitude(body) ?: return null
        }
        return if (negative) -magnitude else magnitude
    }

    /** An unsigned decimal number, or `null` if [body] is not one. */
    private fun parseMagnitude(body: String): Float? {
        var index = 0
        val wholeStart = index
        while (index < body.length && body[index].isAsciiDigit()) index++
        val whole = body.substring(wholeStart, index)
        var fraction = ""
        if (index < body.length && body[index] == '.') {
            val fractionStart = ++index
            while (index < body.length && body[index].isAsciiDigit()) index++
            fraction = body.substring(fractionStart, index)
        }
        if (whole.isEmpty() && fraction.isEmpty()) return null

        var exponent = 0L
        if (index < body.length && (body[index] == 'e' || body[index] == 'E')) {
            index++
            val exponentNegative = index < body.length && body[index] == '-'
            if (index < body.length && (body[index] == '-' || body[index] == '+')) index++
            val exponentStart = index
            while (index < body.length && body[index].isAsciiDigit()) {
                // Clamped: anything past this is already far outside the float range.
                exponent = min(exponent * 10 + (body[index] - '0'), EXPONENT_CLAMP)
                index++
            }
            if (index == exponentStart) return null
            if (exponentNegative) exponent = -exponent
        }
        if (index < body.length && body[index] in SUFFIXES) index++
        if (index != body.length) return null

        // The value is digits * 10^power, with no leading or trailing zeros in digits.
        val allDigits = (whole + fraction).trimStart('0')
        if (allDigits.isEmpty()) return 0f
        val digits = allDigits.trimEnd('0')
        val power = exponent - fraction.length + (allDigits.length - digits.length)
        return nearestFloat(digits, power)
    }

    /** The float nearest digits * 10^[power], ties to even, as Java rounds. */
    private fun nearestFloat(digits: String, power: Long): Float {
        // The value lies in [10^(magnitude - 1), 10^magnitude).
        val magnitude = digits.length + power
        if (magnitude > MAX_MAGNITUDE) return Float.POSITIVE_INFINITY
        if (magnitude < MIN_MAGNITUDE) return 0f

        val nearestDouble = "${digits}e$power".toDouble()
        val float = nearestDouble.toFloat()
        val floatValue = exactValue(float)
        if (floatValue == nearestDouble) return float

        // Narrowing rounds wrongly only when the double is exactly halfway between two floats. The
        // float is non-negative, so its neighbors are one bit pattern away (infinity's lower one is
        // the largest float), which also works where kotlin.math's nextUp is JVM only.
        val neighbor = Float.fromBits(float.toRawBits() + if (floatValue > nearestDouble) -1 else 1)
        val halfway = (floatValue + exactValue(neighbor)) / 2
        if (halfway != nearestDouble) return float
        val order = compareExactly(digits, power, halfway)
        return when {
            order > 0 -> max(float, neighbor)
            order < 0 -> min(float, neighbor)
            else -> float // a true tie; narrowing already rounded it to even
        }
    }

    /** [float] as a `Double`, with infinity standing for 2^128, the next value past the largest float. */
    private fun exactValue(float: Float): Double = if (float.isInfinite()) TWO_TO_THE_128 else float.toDouble()

    /** The sign of digits * 10^[power] minus [value], computed without rounding. */
    private fun compareExactly(digits: String, power: Long, value: Double): Int {
        val bits = value.toRawBits()
        val biasedExponent = ((bits ushr 52) and 0x7FF).toInt()
        val significand = bits and 0x000F_FFFF_FFFF_FFFFL
        val mantissa = if (biasedExponent == 0) significand else significand or (1L shl 52)
        val binaryExponent = if (biasedExponent == 0) -1074 else biasedExponent - 1075
        val expansion = Decimal(mantissa)
        // value = valueDigits * 10^-valueScale
        val valueScale = if (binaryExponent >= 0) {
            expansion.multiplyByPowerOf(2, binaryExponent)
            0
        } else {
            expansion.multiplyByPowerOf(5, -binaryExponent)
            -binaryExponent
        }
        val valueDigits = expansion.toString()
        val wholeDigits = digits.length + power
        val valueWholeDigits = (valueDigits.length - valueScale).toLong()
        if (wholeDigits != valueWholeDigits) return wholeDigits.compareTo(valueWholeDigits)
        val width = max(digits.length, valueDigits.length)
        return digits.padEnd(width, '0').compareTo(valueDigits.padEnd(width, '0'))
    }

    private fun Char.isAsciiDigit(): Boolean = this in '0'..'9'

    private const val SUFFIXES = "fFdD"
    private const val EXPONENT_CLAMP = 1_000_000L

    /** 10^39 is past the largest float (about 3.4e38), so it is infinite. */
    private const val MAX_MAGNITUDE = 39L

    /** Below 10^-46 is under half the smallest float (about 1.4e-45), so it is zero. */
    private const val MIN_MAGNITUDE = -45L

    private const val TWO_TO_THE_128 = 3.4028236692093846E38
}
