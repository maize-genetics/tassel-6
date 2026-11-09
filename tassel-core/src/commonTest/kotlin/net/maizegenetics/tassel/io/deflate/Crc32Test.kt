package net.maizegenetics.tassel.io.deflate

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

class Crc32Test {
    @Test
    fun emptyInputIsZero() {
        val crc = Crc32()
        assertEquals(0, crc.value)
        crc.update(byteArrayOf(1, 2, 3), 1, 1)
        assertEquals(0, crc.value)
    }

    @Test
    fun matchesTheRfc1952CheckValue() {
        val crc = Crc32()
        crc.update("123456789".encodeToByteArray(), 0, 9)
        assertEquals(0xCBF43926.toInt(), crc.value)
    }

    @Test
    fun slicingAgreesWithTheByteLoopAndWithSplitUpdates() {
        val bytes = ByteArray(64) { index -> (index * 17 + 3).toByte() }
        for (size in intArrayOf(0, 1, 7, 8, 9, 15, 16, 17, 64)) {
            val whole = Crc32()
            whole.update(bytes, 0, size)
            val split = Crc32()
            val mid = size / 2
            split.update(bytes, 0, mid)
            split.update(bytes, mid, size)
            assertEquals(whole.value, split.value, "size $size")
        }
    }

    @Test
    fun resetStartsOver() {
        val crc = Crc32()
        val bytes = "123456789".encodeToByteArray()
        crc.update(bytes, 0, bytes.size)
        crc.reset()
        assertEquals(0, crc.value)
        crc.update(bytes, 0, bytes.size)
        assertEquals(0xCBF43926.toInt(), crc.value)
    }

    @Test
    fun rejectsARangeOutsideTheArray() {
        val failure = assertFailsWith<IllegalArgumentException> {
            Crc32().update(byteArrayOf(1), 0, 2)
        }
        assertEquals(true, failure.message?.contains("invalid crc range"), failure.message)
    }
}
