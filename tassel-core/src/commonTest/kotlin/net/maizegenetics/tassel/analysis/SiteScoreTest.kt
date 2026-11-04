package net.maizegenetics.tassel.analysis

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * The byte round trip TASSEL 5 puts every reference probability through. Expected values come from
 * `SiteScoreUtil`: `Math.round(254f * value)` on the way in, and `value / 254f` rounded half up to
 * three decimals on the way out.
 *
 * @see net.maizegenetics.dna.snp.score.SiteScoreUtil
 */
class SiteScoreTest {
    /** The only values an observed genotype produces, and the reason they survive the round trip. */
    @Test
    fun observedProbabilitiesAreFixedPoints() {
        assertEquals(0, SiteScore.encode(0.0f))
        assertEquals(127, SiteScore.encode(0.5f))
        assertEquals(254, SiteScore.encode(1.0f))
        assertEquals(0.0f, SiteScore.quantize(0.0f))
        assertEquals(0.5f, SiteScore.quantize(0.5f))
        assertEquals(1.0f, SiteScore.quantize(1.0f))
    }

    /** `decode` rounds to three decimals, so an imputed mean loses everything beyond them. */
    @Test
    fun imputedValuesAreRoundedToThreeDecimals() {
        assertEquals(0.787f, SiteScore.quantize(0.7873754f))
        assertEquals(0.996f, SiteScore.quantize(0.9955f))
        assertEquals(0.004f, SiteScore.quantize(0.004f))
        assertEquals(0.0f, SiteScore.quantize(0.001f))
    }

    @Test
    fun everyStepIsThreeDecimalsApartAndIncreasing() {
        val steps = (0..254).map { SiteScore.decode(it) }
        assertEquals(0.0f, steps.first())
        assertEquals(1.0f, steps.last())
        steps.zipWithNext { low, high -> assertTrue(high > low, "$high should follow $low") }
        steps.forEach { assertEquals(it, (it * 1000.0f).toInt() / 1000.0f, "not a three-decimal value: $it") }
    }

    @Test
    fun missingStaysMissing() {
        assertEquals(SiteScore.NAN, SiteScore.encode(Float.NaN))
        assertTrue(SiteScore.decode(SiteScore.NAN).isNaN())
        assertTrue(SiteScore.quantize(Float.NaN).isNaN())
    }

    @Test
    fun rejectsScoresOutsideZeroToOne() {
        assertFailsWith<IllegalArgumentException> { SiteScore.encode(-0.001f) }
        assertFailsWith<IllegalArgumentException> { SiteScore.encode(1.001f) }
    }
}
