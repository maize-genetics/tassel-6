package net.maizegenetics.tassel.io

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import kotlinx.coroutines.test.runTest
import kotlinx.io.IOException
import net.maizegenetics.tassel.testing.TestData
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotNull
import kotlin.test.assertTrue

/**
 * Positional reads against the htslib tutorial bgzip, an in-memory copy of those bytes, and a sparse file past 2^31.
 */
class RandomAccessSourceTest {
    @Test
    fun emptyInMemorySourceReadsNothing() = runTest {
        ByteArrayRandomAccessSource(ByteArray(0)).use { source ->
            assertEquals(0L, source.size)
            assertEquals(0, source.readAt(0, ByteArray(4), 0, 0))
            assertEquals(0, source.readAt(0, ByteArray(4), 0, 4))
            val dest = ByteArray(4) { SENTINEL }
            assertEquals(0, source.readAt(PAST_INT, dest, 1, 2))
            assertEquals(SENTINEL, dest[0])
            assertEquals(SENTINEL, dest[3])
        }
    }

    @Test
    fun positionPastIntMaxValueDoesNotTruncate() = runTest {
        val bytes = byteArrayOf(9, 8, 7)
        ByteArrayRandomAccessSource(bytes).use { source ->
            val dest = ByteArray(3) { 1 }
            assertEquals(0, source.readAt(3_000_000_000L, dest, 0, 3))
            assertContentEquals(byteArrayOf(1, 1, 1), dest)
            assertEquals(3, source.readAt(0, dest, 0, 3))
            assertContentEquals(bytes, dest)
        }
    }

    @Test
    fun inMemorySlicesMatchTheBytes() = runTest {
        val bytes = TestData.readBytes(TUTORIAL_GZ)
        ByteArrayRandomAccessSource(bytes).use { source ->
            assertEquals(bytes.size.toLong(), source.size)
            assertSlicesMatch(source, bytes)
        }
    }

    @Test
    fun fileSlicesMatchTheBytes() = runTest {
        val bytes = TestData.readBytes(TUTORIAL_GZ)
        openRandomAccessFile("${TestData.root}/$TUTORIAL_GZ").use { source ->
            assertEquals(bytes.size.toLong(), source.size)
            assertSlicesMatch(source, bytes)
        }
    }

    @Test
    fun inMemoryEdges() = runTest {
        val source = ByteArrayRandomAccessSource(byteArrayOf(1, 2, 3, 4))
        assertEquals(1, source.readAt(0, ByteArray(1), 0, 1))
        assertEdges(source)
    }

    @Test
    fun fileEdges() = runTest {
        val source = openRandomAccessFile("${TestData.root}/$TUTORIAL_GZ")
        assertTrue(source.size > 0)
        assertEdges(source)
    }

    @Test
    fun missingFileThrows() {
        val error = assertFailsWith<IOException> {
            openRandomAccessFile("${TestData.root}/htslib/no-such-file")
        }
        assertContains(error.message ?: "", "no-such-file")
    }

    @Test
    fun concurrentReadsMatch() = runTest {
        checkConcurrentReads()
    }

    @Test
    fun readPastTwoToTheThirtyOne() = runTest {
        val boundary = byteArrayOf(1, 2, 3, 4, 5, 6, 7, 8)
        val tail = "tassel-tail".encodeToByteArray()
        val path = createSparseFixture(
            SPARSE_SIZE,
            listOf(PAST_INT to boundary, SPARSE_SIZE - tail.size to tail),
        )
        try {
            openRandomAccessFile(path).use { source ->
                assertEquals(SPARSE_SIZE, source.size)
                val crossed = ByteArray(8)
                assertEquals(8, source.readAt(PAST_INT - 4, crossed, 0, 8))
                assertContentEquals(byteArrayOf(0, 0, 0, 0, 1, 2, 3, 4), crossed)
                val atBoundary = ByteArray(8)
                assertEquals(8, source.readAt(PAST_INT, atBoundary, 0, 8))
                assertContentEquals(boundary, atBoundary)
                assertEquals(0, source.readAt(PAST_INT, ByteArray(1), 0, 0))
                val one = ByteArray(1)
                assertEquals(1, source.readAt(PAST_INT, one, 0, 1))
                assertEquals(1, one[0].toInt())
                val end = ByteArray(tail.size + 4)
                assertEquals(tail.size, source.readAt(SPARSE_SIZE - tail.size, end, 0, tail.size))
                assertContentEquals(tail, end.copyOf(tail.size))
                assertEquals(1, source.readAt(SPARSE_SIZE - 1, one, 0, 1))
                assertEquals(tail.last(), one[0])
                assertEquals(2, source.readAt(SPARSE_SIZE - 2, end, 1, 4))
                assertEquals(0, source.readAt(SPARSE_SIZE, end, 0, 4))
            }
        } finally {
            deleteFixture(path)
        }
    }

    internal suspend fun checkConcurrentReads() {
        val bytes = TestData.readBytes(TUTORIAL_GZ)
        val path = "${TestData.root}/$TUTORIAL_GZ"
        ByteArrayRandomAccessSource(bytes).use { memory -> assertConcurrent(memory, bytes) }
        openRandomAccessFile(path).use { file -> assertConcurrent(file, bytes) }
    }

    private suspend fun assertConcurrent(source: RandomAccessSource, bytes: ByteArray) {
        val shared = minOf(32, bytes.size)
        coroutineScope {
            List(WORKERS) { worker ->
                async {
                    val sharedDest = ByteArray(shared)
                    assertEquals(shared, source.readAt(0, sharedDest, 0, shared))
                    assertContentEquals(bytes.copyOfRange(0, shared), sharedDest)
                    val rng = Random(SEED + worker)
                    repeat(READS_PER_WORKER) {
                        val position = rng.nextLong(bytes.size.toLong())
                        val length = rng.nextInt(1, 65).coerceAtMost(bytes.size - position.toInt())
                        val buffer = ByteArray(length)
                        assertEquals(length, source.readAt(position, buffer, 0, length), "worker $worker at $position")
                        val start = position.toInt()
                        assertContentEquals(bytes.copyOfRange(start, start + length), buffer)
                    }
                }
            }.awaitAll()
        }
    }

    private suspend fun assertSlicesMatch(source: RandomAccessSource, bytes: ByteArray) {
        val explicit = listOf(
            0L to 0,
            0L to 1,
            (bytes.size - 1).toLong() to 1,
            (bytes.size - 1).toLong() to 2,
            bytes.size.toLong() to 0,
            bytes.size.toLong() to 1,
            0L to minOf(64, bytes.size),
        )
        for ((position, length) in explicit) {
            assertSlice(source, bytes, position, length, destinationOffset = 0)
        }
        val rng = Random(SEED)
        repeat(SLICE_COUNT) { i ->
            val position = rng.nextLong(bytes.size.toLong() + 1)
            val room = if (position >= bytes.size) {
                MAX_SLICE
            } else {
                minOf(MAX_SLICE, bytes.size - position.toInt())
            }
            val length = when {
                i % 50 == 0 -> 0
                i % 50 == 1 -> 1
                else -> rng.nextInt(room + 1)
            }
            val offset = if (i % 7 == 0) 3 else 0
            assertSlice(source, bytes, position, length, offset)
        }
    }

    private suspend fun assertSlice(
        source: RandomAccessSource,
        bytes: ByteArray,
        position: Long,
        length: Int,
        destinationOffset: Int,
    ) {
        val dest = ByteArray(destinationOffset + length + 1) { SENTINEL }
        val n = source.readAt(position, dest, destinationOffset, length)
        val available = if (position >= bytes.size) {
            0
        } else {
            minOf(length.toLong(), bytes.size.toLong() - position).toInt()
        }
        assertEquals(available, n, "position=$position length=$length offset=$destinationOffset")
        if (n > 0) {
            val start = position.toInt()
            for (i in 0 until n) {
                assertEquals(bytes[start + i], dest[destinationOffset + i], "byte $i at $position")
            }
        }
        if (destinationOffset > 0) assertEquals(SENTINEL, dest[0])
        assertEquals(SENTINEL, dest[destinationOffset + length])
    }

    private suspend fun assertEdges(source: RandomAccessSource) {
        assertInvalid { source.readAt(-1, ByteArray(1), 0, 0) }
        assertInvalid { source.readAt(-1, ByteArray(1), 0, 1) }
        assertInvalid { source.readAt(0, ByteArray(1), 0, -1) }
        assertInvalid { source.readAt(0, ByteArray(1), -1, 0) }
        assertInvalid { source.readAt(0, ByteArray(1), 1, 1) }
        assertInvalid { source.readAt(0, ByteArray(1), 0, 2) }
        source.close()
        source.close()
        val sizeError = try {
            source.size
            null
        } catch (e: IllegalStateException) {
            e
        }
        assertNotNull(sizeError)
        assertContains(sizeError.message ?: "", "closed")
        val readError = try {
            source.readAt(0, ByteArray(1), 0, 1)
            null
        } catch (e: IllegalStateException) {
            e
        }
        assertNotNull(readError)
        assertContains(readError.message ?: "", "closed")
    }

    private suspend fun assertInvalid(block: suspend () -> Unit) {
        val error = try {
            block()
            null
        } catch (e: IllegalArgumentException) {
            e
        }
        assertNotNull(error)
    }

    private companion object {
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val SEED = 86
        const val SLICE_COUNT = 400
        const val WORKERS = 8
        const val READS_PER_WORKER = 50
        const val MAX_SLICE = 4096
        const val SPARSE_SIZE = 3L * 1024 * 1024 * 1024
        const val PAST_INT = 1L shl 31
        val SENTINEL = 0x5A.toByte()
    }
}
