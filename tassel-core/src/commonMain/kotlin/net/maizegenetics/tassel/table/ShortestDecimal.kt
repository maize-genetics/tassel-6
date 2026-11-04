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
    else -> javaLayout(value.toString())
}

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
