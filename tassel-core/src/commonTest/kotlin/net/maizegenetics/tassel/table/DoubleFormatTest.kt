package net.maizegenetics.tassel.table

import kotlin.test.Test
import kotlin.test.assertEquals

/**
 * Expected values come from TASSEL 5 `v5.2.98`'s own `net.maizegenetics.util.DoubleFormat`, run on
 * a JDK 21 and printed one value per line.
 */
class DoubleFormatTest {
    @Test
    fun formatsTheExamplesFromThePcaGoldens() {
        assertEquals("0.00176", DoubleFormat.format(0.0017601161965581133))
        assertEquals("-1.5526E-2", DoubleFormat.format(-0.015526324308868099))
        assertEquals("9.1313E-4", DoubleFormat.format(9.131257175122849E-4))
        assertEquals("281", DoubleFormat.format(281.0))
    }

    @Test
    fun writesNonNumbersAsTassel5Does() {
        assertEquals("NaN", DoubleFormat.format(Double.NaN))
        assertEquals("NaN", DoubleFormat.format(null))
        // TASSEL 5 drops the sign of an infinity.
        assertEquals("Infinity", DoubleFormat.format(Double.POSITIVE_INFINITY))
        assertEquals("Infinity", DoubleFormat.format(Double.NEGATIVE_INFINITY))
        assertEquals("0", DoubleFormat.format(0.0))
        assertEquals("-0", DoubleFormat.format(-0.0))
    }

    @Test
    fun switchesToScientificBelowThousandthsAndAtAMillion() {
        assertEquals("0.001", DoubleFormat.format(0.001))
        assertEquals("1.0E-3", DoubleFormat.format(9.999999E-4))
        assertEquals("999999.5", DoubleFormat.format(999999.5))
        // Rounding inside the plain range can reach a million without switching formats.
        assertEquals("1000000", DoubleFormat.format(999999.9999999))
        assertEquals("1.0E6", DoubleFormat.format(1_000_000.0))
        assertEquals("1.7383E6", DoubleFormat.format(1738266.0))
        assertEquals("869133", DoubleFormat.format(869133.0))
    }

    @Test
    fun writesEveryNegativeNumberInScientificNotation() {
        assertEquals("-1.0E0", DoubleFormat.format(-1.0))
        assertEquals("-1.0E-10", DoubleFormat.format(-1.0E-10))
        assertEquals("-3.1088E-4", DoubleFormat.format(-3.1088287E-4))
    }

    @Test
    fun keepsFiveDecimalsAndDropsTrailingZeros() {
        assertEquals("0.28342", DoubleFormat.format(0.2834169223812696))
        assertEquals("0.0375", DoubleFormat.format(0.037495732474496314))
        assertEquals("0.001", DoubleFormat.format(0.0010015255657604791))
        assertEquals("123456.789", DoubleFormat.format(123456.789))
        assertEquals("1", DoubleFormat.format(0.9999999))
    }

    @Test
    fun keepsFiveSignificantDigitsInScientificNotation() {
        assertEquals("6.4087E-4", DoubleFormat.format(6.4087E-4))
        assertEquals("9.6533E-4", DoubleFormat.format(9.653298171856321E-4))
        assertEquals("1.25E-5", DoubleFormat.format(1.25E-5))
        assertEquals("1.0E-300", DoubleFormat.format(1.0E-300))
        assertEquals("1.0E300", DoubleFormat.format(1.0E300))
        assertEquals("1.7977E308", DoubleFormat.format(Double.MAX_VALUE))
        assertEquals("2.2251E-308", DoubleFormat.format(2.2250738585072014E-308))
    }

    @Test
    fun roundsTheExactBinaryValue() {
        // The stored double is just below the 0.123455 tie, so it rounds down even though the
        // shortest decimal that prints it ends in 5.
        assertEquals("0.12345", DoubleFormat.format(0.123455))
        // Just above the 0.123465 tie, so it rounds up even though the digit before is even.
        assertEquals("0.12347", DoubleFormat.format(0.123465))
        assertEquals("0.12346", DoubleFormat.format(0.1234550000000001))
        // 1/64 and 3/64 are exact in binary and tie at the fifth decimal, so half-even decides:
        // 0.015625 keeps the even 2 and 0.046875 rounds the odd 7 up.
        assertEquals("0.01562", DoubleFormat.format(0.015625))
        assertEquals("0.04688", DoubleFormat.format(0.046875))
        assertEquals("-2.0312E-1", DoubleFormat.format(-0.203125))
    }

    @Test
    fun roundsAWholeNumberTieUpAsTassel5Does() {
        // The JDK's binary-to-decimal fast path does not report that a whole number's digits are
        // exact, so DecimalFormat rounds these ties up rather than to even. TASSEL 5 holds counts
        // as whole Doubles, so the summary tables depend on it.
        assertEquals("1.2345E6", DoubleFormat.format(1234450.0))
        assertEquals("7.6435E7", DoubleFormat.format(76434500.0))
        assertEquals("-2.0313E5", DoubleFormat.format(-203125.0))
        // Nothing changes when the digit before the tie is already odd.
        assertEquals("1.2346E6", DoubleFormat.format(1234550.0))
    }
}
