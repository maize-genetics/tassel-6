// Java's Float.toString and Double.toString: the shortest decimal that round-trips, written plainly
// for 1e-3 <= |value| < 1e7 and in scientific notation outside that range.
//
// TASSEL 5 writes distance-matrix cells with String.valueOf(float) and quantitative annotations with
// Double.toString, so a port that writes those files needs the same spelling. Kotlin produces the same
// digits on every target but lays them out differently off the JVM: Kotlin/Wasm writes 0.00031861675
// and 3.4028235e+38 where Java writes 3.1861675E-4 and 3.4028235E38. These functions take the
// platform's digits and re-emit them Java's way.
//
// Table cells are a different problem: TASSEL rounds those to five significant digits, which is
// DoubleFormat in this package.

package net.maizegenetics.tassel.table

/** [value] as Java's `Float.toString` writes it. */
fun floatToString(value: Float): String = when {
    value.isNaN() -> NAN
    value == Float.POSITIVE_INFINITY -> INFINITY
    value == Float.NEGATIVE_INFINITY -> "-$INFINITY"
    else -> javaLayout(value.toString())
}

/** [value] as Java's `Double.toString` writes it. */
fun doubleToString(value: Double): String = when {
    value.isNaN() -> NAN
    value == Double.POSITIVE_INFINITY -> INFINITY
    value == Double.NEGATIVE_INFINITY -> "-$INFINITY"
    else -> javaLayout(shortestRoundTrip(value))
}

/**
 * The platform's digits for [value], shortened when they are not the shortest. Kotlin/Wasm writes `7668.1141812907445`
 * where Java writes `7668.114181290744`: both read back as the same double, and Java's is shorter. A string of 15
 * digits or fewer is always the shortest, so only longer ones are searched: the fewest digits that read back as
 * [value], the one nearest the exact value when two do, as Java does.
 *
 * "Reads back" is decided exactly, in integers, by [readsBackAs], not by `toDouble`: Kotlin/Wasm's parser is off by an
 * ulp for some 16 and 17 digit strings.
 */
private fun shortestRoundTrip(value: Double): String {
    val platform = value.toString()
    val mantissa = platform.trimStart('-').substringBefore('e').substringBefore('E')
    val significant = mantissa.replace(".", "").trim('0').length
    if (significant <= MAX_ALWAYS_SHORTEST) return platform
    val sign = if (value < 0.0) "-" else ""
    for (count in 1 until significant) {
        for (candidate in DoubleFormat.significantDigits(kotlin.math.abs(value), count)) {
            val digits = candidate.digits
            if (!readsBackAs(value, digits.toLong(), candidate.power - digits.length + 1)) continue
            return sign + digits[0] + "." + digits.drop(1).ifEmpty { "0" } + "E" + candidate.power
        }
    }
    return platform
}

/**
 * Whether the decimal `digits * 10^exponent` lies in the interval of reals that round to [value]: between the
 * midpoints to its neighbours, a midpoint itself counting only when [value]'s significand is even. Exact: both
 * sides are scaled to whole numbers and compared digit by digit.
 */
internal fun readsBackAs(value: Double, digits: Long, exponent: Int): Boolean {
    val bits = kotlin.math.abs(value).toRawBits()
    val biased = ((bits ushr 52) and 0x7FF).toInt()
    val fraction = bits and 0x000F_FFFF_FFFF_FFFFL
    val significand = if (biased == 0) fraction else fraction or (1L shl 52)
    val binaryExponent = if (biased == 0) -1074 else biased - 1075
    // Below a power of two the neighbour is twice as close, so the interval is half as wide there.
    val narrowBelow = biased > 1 && fraction == 0L
    val above = compare(digits, exponent, 2 * significand + 1, binaryExponent - 1)
    val below = if (narrowBelow) {
        compare(digits, exponent, 4 * significand - 1, binaryExponent - 2)
    } else {
        compare(digits, exponent, 2 * significand - 1, binaryExponent - 1)
    }
    val even = significand % 2 == 0L
    return (above < 0 || (above == 0 && even)) && (below > 0 || (below == 0 && even))
}

/** The sign of `a * 10^tenPower - b * 2^twoPower`, for positive [a] and [b]. */
private fun compare(a: Long, tenPower: Int, b: Long, twoPower: Int): Int {
    // Multiply both sides by 10^max(-tenPower, 0) * 2^max(-twoPower, 0) so that no power is negative.
    val left = Decimal(a)
    if (tenPower > 0) {
        left.multiplyByPowerOf(2, tenPower)
        left.multiplyByPowerOf(5, tenPower)
    }
    if (twoPower < 0) left.multiplyByPowerOf(2, -twoPower)
    val right = Decimal(b)
    if (twoPower > 0) right.multiplyByPowerOf(2, twoPower)
    if (tenPower < 0) {
        right.multiplyByPowerOf(2, -tenPower)
        right.multiplyByPowerOf(5, -tenPower)
    }
    val x = left.toString()
    val y = right.toString()
    return if (x.length != y.length) x.length.compareTo(y.length) else x.compareTo(y)
}

private const val MAX_ALWAYS_SHORTEST = 15

private const val NAN = "NaN"
private const val INFINITY = "Infinity"

/** Smallest power of ten Java writes plainly; below it `Float.toString` switches to `1.0E-4`. */
private const val MIN_PLAIN_POWER = -3

/** Smallest power of ten Java writes in scientific notation: `1.0E7` rather than `10000000.0`. */
private const val MAX_PLAIN_POWER = 7

/**
 * Rewrites a finite platform [Float.toString] or [Double.toString] in Java's layout. The digits are
 * taken as given; only where the decimal point goes and how the exponent is spelled change.
 */
private fun javaLayout(platform: String): String {
    val negative = platform.startsWith("-")
    val sign = if (negative) "-" else ""
    val body = if (negative) platform.substring(1) else platform

    val exponentAt = body.indexOfFirst { it == 'e' || it == 'E' }
    val mantissa = if (exponentAt < 0) body else body.substring(0, exponentAt)
    val exponent = if (exponentAt < 0) 0 else body.substring(exponentAt + 1).removePrefix("+").toInt()

    val point = mantissa.indexOf('.')
    val wholeLength = if (point < 0) mantissa.length else point
    val allDigits = if (point < 0) mantissa else mantissa.substring(0, point) + mantissa.substring(point + 1)

    val firstDigit = allDigits.indexOfFirst { it != '0' }
    if (firstDigit < 0) return "${sign}0.0"

    // power is the power of ten of the leading digit, so the value is digits[0].digits[1..] * 10^power.
    val digits = allDigits.substring(firstDigit).trimEnd('0')
    val power = wholeLength - firstDigit - 1 + exponent
    val laidOut = when {
        power < MIN_PLAIN_POWER || power >= MAX_PLAIN_POWER ->
            "${digits[0]}.${digits.drop(1).ifEmpty { "0" }}E$power"
        power < 0 -> "0." + "0".repeat(-power - 1) + digits
        else -> digits.take(power + 1).padEnd(power + 1, '0') + "." + digits.drop(power + 1).ifEmpty { "0" }
    }
    return sign + laidOut
}
