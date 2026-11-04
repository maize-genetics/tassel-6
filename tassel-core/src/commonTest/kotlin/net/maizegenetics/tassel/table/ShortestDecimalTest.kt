package net.maizegenetics.tassel.table

import kotlin.test.Test
import kotlin.test.assertEquals

/**
 * [floatToString] and [doubleToString] against Java's `Float.toString` and `Double.toString`, which
 * is how TASSEL 5 spells distance-matrix cells and quantitative annotations. The expected strings
 * were printed by a JVM; Kotlin itself only agrees on the JVM, so this test is what keeps the other
 * targets honest.
 */
class ShortestDecimalTest {
    @Test
    fun writesFloatsAsJavaDoes() {
        val expected = mapOf(
            1.7978103f to "1.7978103",
            0.03960289f to "0.03960289",
            -0.008272289f to "-0.008272289",
            -3.1861675e-4f to "-3.1861675E-4",
            3.1348677e-6f to "3.1348677E-6",
            8.1197644e-4f to "8.1197644E-4",
            0.0f to "0.0",
            -0.0f to "-0.0",
            1.0f to "1.0",
            -1.0f to "-1.0",
            100.0f to "100.0",
            0.1f to "0.1",
            1.0e-3f to "0.001",
            9.99e-4f to "9.99E-4",
            9999999.0f to "9999999.0",
            1.0e7f to "1.0E7",
            1.23456e10f to "1.23456E10",
            Float.MAX_VALUE to "3.4028235E38",
            Float.NaN to "NaN",
            Float.POSITIVE_INFINITY to "Infinity",
            Float.NEGATIVE_INFINITY to "-Infinity",
        )
        assertEquals(expected, expected.keys.associateWith { floatToString(it) })
    }

    @Test
    fun writesDoublesAsJavaDoes() {
        val expected = mapOf(
            918.7297856754799 to "918.7297856754799",
            0.0 to "0.0",
            -0.0 to "-0.0",
            1.0 to "1.0",
            0.001 to "0.001",
            9.99e-4 to "9.99E-4",
            9999999.0 to "9999999.0",
            1.0e7 to "1.0E7",
            1.0 / 3.0 to "0.3333333333333333",
            Double.MAX_VALUE to "1.7976931348623157E308",
            Double.NaN to "NaN",
            Double.POSITIVE_INFINITY to "Infinity",
            Double.NEGATIVE_INFINITY to "-Infinity",
        )
        assertEquals(expected, expected.keys.associateWith { doubleToString(it) })
    }

    /**
     * The one place the digits are not Java's. For the smallest subnormals, Java writes a digit more
     * than it needs (`1.4E-45` for `Float.MIN_VALUE`, `4.9E-324` for `Double.MIN_VALUE`), while every
     * other target writes the genuinely shortest decimal. Both parse back to the same value, and no
     * analysis produces a number near there, so the layout is what matters and the digits are left
     * to the platform rather than reimplemented.
     */
    @Test
    fun keepsThePlatformDigitsForTheSmallestSubnormals() {
        assertEquals(Float.MIN_VALUE, floatToString(Float.MIN_VALUE).toFloat())
        assertEquals(Double.MIN_VALUE, doubleToString(Double.MIN_VALUE).toDouble())
    }

    @Test
    fun everyFormattedFloatParsesBackToItself() {
        var bits = 1
        repeat(2000) {
            bits = bits * 1664525 + 1013904223
            val value = Float.fromBits(bits)
            if (!value.isNaN() && !value.isInfinite()) {
                assertEquals(value, floatToString(value).toFloat(), "round trip of $value")
            }
        }
    }
}
