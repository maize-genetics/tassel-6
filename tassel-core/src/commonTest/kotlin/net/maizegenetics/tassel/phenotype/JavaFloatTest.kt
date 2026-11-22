package net.maizegenetics.tassel.phenotype

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull

/**
 * [JavaFloat.parse] against bit patterns `java.lang.Float.parseFloat` (Java 21) produced for the
 * same text, so every target is held to the JVM's answer. The corpus has the tutorial's values,
 * Java's suffixes and spellings, the edges of the float range, and decimals a hair either side of
 * the halfway point between two floats, where parsing a `Double` and narrowing it rounds wrongly.
 * `JavaFloatJvmTest` compares against `Float.parseFloat` itself on many more inputs.
 */
class JavaFloatTest {
    @Test
    fun matchesJavasFloatParseFloatBitForBit() {
        val wrong = corpus.mapNotNull { (text, bits) ->
            val parsed = JavaFloat.parse(text)?.toRawBits()
            if (parsed ==
                bits
            ) {
                null
            } else {
                "'$text': expected 0x${bits.toUInt().toString(16)}, got ${parsed?.toUInt()?.toString(16)}"
            }
        }
        assertEquals(emptyList(), wrong)
    }

    @Test
    fun rejectsWhatJavaRejects() {
        val notNumbers = listOf(
            "", " ", "abc", "1e", "e5", ".", "-", "+", "1.2.3", "1,5", "--1", "1f5", "1ff", "Inf", "nan", "NAN",
            "infinity", "1_000", "1e+", "1 2", "NaNf", "Infinityf",
        )
        notNumbers.forEach { assertNull(JavaFloat.parse(it), "'$it' should not parse") }
    }

    /** Java also reads hexadecimal floats; TASSEL 6 does not (ADR-0023). */
    @Test
    fun doesNotReadHexadecimalFloats() {
        assertNull(JavaFloat.parse("0x1p3"))
    }

    /** The digits of 2^-150, half the smallest float, all but the last two. */
    private val halfSmallest = "7.0064923216240853546186479164495806564013097093825788587853" +
        "41419448955413429303007433190941810607910156"

    private val corpus: List<Pair<String, Int>> = listOf(
        "1.1" to 0x3F8CCCCD.toInt(),
        "64.75" to 0x42818000.toInt(),
        "-999" to 0xC479C000.toInt(),
        "-99" to 0xC2C60000.toInt(),
        "37.897" to 0x42179687.toInt(),
        "32.21933" to 0x4200E098.toInt(),
        "0.014" to 0x3C656042.toInt(),
        "92.25" to 0x42B88000.toInt(),
        "65" to 0x42820000.toInt(),
        "3.1861675E-4" to 0x39A70C04.toInt(),
        "0" to 0x00000000.toInt(),
        "-0" to 0x80000000.toInt(),
        "-0.0" to 0x80000000.toInt(),
        "0e999" to 0x00000000.toInt(),
        "1e-999" to 0x00000000.toInt(),
        "1e39" to 0x7F800000.toInt(),
        "-1e39" to 0xFF800000.toInt(),
        "3.4028235e38" to 0x7F7FFFFF.toInt(),
        "3.4028236e38" to 0x7F800000.toInt(),
        "3.40282356779733661637539395458142568448e38" to 0x7F800000.toInt(),
        "3.4028235677973366163753939545814256844e38" to 0x7F7FFFFF.toInt(),
        "3.40282356779733661637539395458142568449e38" to 0x7F800000.toInt(),
        "1.401298464324817e-45" to 0x00000001.toInt(),
        "1e-45" to 0x00000001.toInt(),
        "7.006492321624085e-46" to 0x00000000.toInt(),
        "7.0064923216240854e-46" to 0x00000001.toInt(),
        halfSmallest + "25e-46" to
            0x00000000.toInt(),
        halfSmallest + "26e-46" to
            0x00000001.toInt(),
        "1.000000059604644775390625" to 0x3F800000.toInt(),
        "1.0000000596046447753906250001" to 0x3F800001.toInt(),
        "1.0000000596046447753906249999" to 0x3F800000.toInt(),
        "1.00000017881393432617187499" to 0x3F800001.toInt(),
        "1.000000178813934326171875" to 0x3F800002.toInt(),
        "1.00000017881393432617187501" to 0x3F800002.toInt(),
        "8388609.5" to 0x4B000002.toInt(),
        "8388610.5" to 0x4B000002.toInt(),
        "16777217" to 0x4B800000.toInt(),
        "16777219" to 0x4B800002.toInt(),
        "33554435" to 0x4C000001.toInt(),
        "0.1" to 0x3DCCCCCD.toInt(),
        "0.2" to 0x3E4CCCCD.toInt(),
        "0.3" to 0x3E99999A.toInt(),
        "1e10" to 0x501502F9.toInt(),
        "1.17549435E-38" to 0x00800000.toInt(),
        "1.1754942E-38" to 0x007FFFFF.toInt(),
        "2.5e-44" to 0x00000012.toInt(),
        "123456789012345678901234567890" to 0x6FC77488.toInt(),
        "0.000000000000000000000000000000001" to 0x08A6274C.toInt(),
        "1.5f" to 0x3FC00000.toInt(),
        "2D" to 0x40000000.toInt(),
        "7F" to 0x40E00000.toInt(),
        "3d" to 0x40400000.toInt(),
        ".5" to 0x3F000000.toInt(),
        "5." to 0x40A00000.toInt(),
        "+.5e+1" to 0x40A00000.toInt(),
        "-1.e-2" to 0xBC23D70A.toInt(),
        "  12 " to 0x41400000.toInt(),
        "\t42\n" to 0x42280000.toInt(),
        "NaN" to 0x7FC00000.toInt(),
        "-Infinity" to 0xFF800000.toInt(),
        "+Infinity" to 0x7F800000.toInt(),
        "Infinity" to 0x7F800000.toInt(),
        "00000.000001234500000" to 0x35A5B119.toInt(),
        "1E+0" to 0x3F800000.toInt(),
        "9.999999e-39" to 0x006CE3EE.toInt(),
        "4.2e-45" to 0x00000003.toInt(),
    )
}
