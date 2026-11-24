package net.maizegenetics.tassel.io

import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

/**
 * [bcfGenotype] against [bcfTextGenotype] for every integer width the genotype block can use.
 *
 * The text function is the reader that shipped before the direct path, so a match is the same call byte or the
 * same exception. A short payload must fail as a truncated record before a genotype error is possible.
 */
class BcfGenotypeTest {
    @Test
    fun directDecodeMatchesTheTextPath() {
        val alleleCounts = intArrayOf(1, 2, 4, 10, 11)
        for (type in intArrayOf(TYPE_INT8, TYPE_INT16, TYPE_INT32)) {
            val grid = grid(type)
            for (perSample in 0..3) {
                for (alleleCount in alleleCounts) {
                    val alleles = ByteArray(alleleCount) { it.toByte() }
                    val values = IntArray(perSample)
                    each(grid, values, 0) { sample ->
                        val block = BcfCursor(encode(type, perSample, sample)).intBlock()
                        val expected = outcome { bcfTextGenotype(block, 0, perSample, alleles, SITE) }
                        val actual = outcome { bcfGenotype(block, 0, perSample, alleles, SITE) }
                        if (expected != actual) {
                            assertEquals(
                                expected,
                                actual,
                                "type=$type alleles=$alleleCount values=${sample.contentToString()}",
                            )
                        }
                    }
                }
            }
        }
    }

    @Test
    fun randomBlockUsesTheSampleStride() {
        val perSample = 2
        val samples = 1_000
        val alleles = ByteArray(4) { it.toByte() }
        val pool = intArrayOf(-128, -127, 0, 1, 2, 3, 4, 5, 20, 22, 23, 127, -1, -5)
        val random = Random(SEED)
        val values = IntArray(samples * perSample) { pool[random.nextInt(pool.size)] }
        values[0] = 2
        values[1] = 2
        values[values.lastIndex - 1] = 2
        values[values.lastIndex] = 4
        val block = BcfCursor(encode(TYPE_INT8, perSample, values)).intBlock(samples)
        assertEquals(samples * perSample, block.count)
        for (sample in 0 until samples) {
            val from = sample * perSample
            val expected = outcome { bcfTextGenotype(block, from, perSample, alleles, SITE) }
            val actual = outcome { bcfGenotype(block, from, perSample, alleles, SITE) }
            if (expected != actual) {
                assertEquals(expected, actual, "sample=$sample from=$from")
            }
        }
    }

    @Test
    fun truncatedPayloadThrowsBeforeACallError() {
        // Two int8 values were declared and one byte is present. That byte is end-of-vector, which would be an
        // empty genotype if it were decoded before the length check.
        assertTruncated(byteArrayOf(0x21, 0x81.toByte()), times = 1)
        // One int16 was declared and only the first byte is present.
        assertTruncated(byteArrayOf(0x12, 0x00), times = 1)
        // One int32 was declared and three of the four bytes are present.
        assertTruncated(byteArrayOf(0x13, 0x00, 0x00, 0x00), times = 1)
        // Diploid int8 for two samples needs four bytes. One byte must not become a call error.
        assertTruncated(byteArrayOf(0x21, 0x81.toByte()), times = 2)
    }

    private fun assertTruncated(bytes: ByteArray, times: Int) {
        val failure = assertFailsWith<IllegalArgumentException> {
            val block = BcfCursor(bytes).intBlock(times)
            bcfGenotype(block, 0, 1, byteArrayOf(0), SITE)
        }
        assertEquals(TRUNCATED_BCF_RECORD, failure.message)
    }

    private fun each(grid: IntArray, values: IntArray, slot: Int, body: (IntArray) -> Unit) {
        if (slot == values.size) {
            body(values)
            return
        }
        for (value in grid) {
            values[slot] = value
            each(grid, values, slot + 1, body)
        }
    }

    private fun grid(type: Int): IntArray {
        val values = ArrayList<Int>()
        when (type) {
            TYPE_INT8 -> {
                values += INT8_MISSING
                values += INT8_VECTOR_END
            }
            TYPE_INT16 -> {
                values += INT16_MISSING
                values += INT16_VECTOR_END
            }
            else -> {
                values += INT32_MISSING
                values += INT32_VECTOR_END
            }
        }
        values += 0
        values += 1
        for (allele in 0..10) {
            values += (allele + 1) shl 1
            values += ((allele + 1) shl 1) or 1
        }
        // Allele 100 is a three-digit token. It does not fit in int8.
        if (type != TYPE_INT8) {
            values += (100 + 1) shl 1
            values += ((100 + 1) shl 1) or 1
        }
        values += -1
        values += -2
        values += -5
        when (type) {
            TYPE_INT8 -> {
                values += 127
                values += -126
            }
            TYPE_INT16 -> {
                values += 32767
                values += -32766
            }
            else -> {
                values += Int.MAX_VALUE
                values += Int.MIN_VALUE + 2
            }
        }
        return values.distinct().toIntArray()
    }

    private fun encode(type: Int, perSample: Int, values: IntArray): ByteArray {
        val width = when (type) {
            TYPE_INT8 -> 1
            TYPE_INT16 -> 2
            else -> 4
        }
        val out = ByteArray(1 + values.size * width)
        out[0] = ((perSample shl 4) or type).toByte()
        var at = 1
        for (value in values) {
            out[at++] = value.toByte()
            if (width > 1) out[at++] = (value ushr 8).toByte()
            if (width > 2) {
                out[at++] = (value ushr 16).toByte()
                out[at++] = (value ushr 24).toByte()
            }
        }
        return out
    }

    private fun outcome(body: () -> Byte): Outcome = try {
        Outcome(body(), null)
    } catch (error: IllegalArgumentException) {
        Outcome(null, error.message)
    } catch (error: Throwable) {
        Outcome(null, "${error::class.simpleName}: ${error.message}")
    }

    private data class Outcome(val value: Byte?, val failure: String?)

    private companion object {
        const val SITE = "chr1:1"
        const val SEED = 145
        const val TYPE_INT8 = 1
        const val TYPE_INT16 = 2
        const val TYPE_INT32 = 3
        const val INT8_MISSING = -128
        const val INT8_VECTOR_END = -127
        const val INT16_MISSING = -32768
        const val INT16_VECTOR_END = -32767
        const val INT32_MISSING = Int.MIN_VALUE
        const val INT32_VECTOR_END = Int.MIN_VALUE + 1
    }
}
