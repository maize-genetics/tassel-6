package net.maizegenetics.tassel.phenotype

import java.math.BigDecimal
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals

/**
 * [JavaFloat.parse] against `java.lang.Float.parseFloat` itself, on tens of thousands of inputs:
 * every kind of float written by `Float.toString`, random decimals across the float range, and
 * decimals exactly at, just above, and just below the halfway point between two adjacent floats.
 */
class JavaFloatJvmTest {
    private val random = Random(20261122)

    @Test
    fun agreesWithFloatParseFloatOnFloatToStringOutput() {
        val texts = List(20_000) { java.lang.Float.toString(Float.fromBits(random.nextInt())) }
        assertAgrees(texts)
    }

    @Test
    fun agreesWithFloatParseFloatOnRandomDecimals() {
        val texts = List(20_000) {
            val digits = (1..random.nextInt(1, 20)).joinToString("") { random.nextInt(10).toString() }
            val point = random.nextInt(digits.length + 1)
            val sign = if (random.nextBoolean()) "-" else ""
            "$sign${digits.substring(0, point)}.${digits.substring(point)}e${random.nextInt(-60, 45)}"
        }
        assertAgrees(texts)
    }

    @Test
    fun agreesWithFloatParseFloatAtHalfwayPoints() {
        val texts = ArrayList<String>()
        repeat(5_000) {
            // Positive, finite, below the largest float, normal or subnormal.
            val float = Float.fromBits(random.nextInt(0, 0x7F7F_FFFF))
            val halfway = (BigDecimal(float.toDouble()) + BigDecimal(Math.nextUp(float).toDouble()))
                .divide(BigDecimal(2))
            val nudge = BigDecimal.ONE.movePointLeft(halfway.scale() + 3)
            texts += halfway.toString()
            texts += (halfway + nudge).toString()
            texts += (halfway - nudge).toString()
        }
        assertAgrees(texts)
    }

    @Test
    fun rejectsExactlyWhatFloatParseFloatRejects() {
        val texts = listOf("", "abc", "1e", ".", "1.2.3", "NaNf", "Infinityf", "nan", "1 2", "1ff", "--1", "+")
        assertAgrees(texts)
    }

    private fun assertAgrees(texts: List<String>) {
        val disagreements = texts.mapNotNull { text ->
            val expected = runCatching { java.lang.Float.parseFloat(text) }.getOrNull()?.toRawBits()
            val actual = JavaFloat.parse(text)?.toRawBits()
            if (expected == actual) null else "'$text': Java $expected, JavaFloat $actual"
        }
        assertEquals(emptyList(), disagreements.take(10), "${disagreements.size} of ${texts.size} disagree")
    }
}
