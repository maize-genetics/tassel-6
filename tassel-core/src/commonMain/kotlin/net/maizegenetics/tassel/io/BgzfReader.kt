package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.io.deflate.BGZF_MAX_BLOCK
import net.maizegenetics.tassel.io.deflate.BGZF_MAX_UNCOMPRESSED
import net.maizegenetics.tassel.io.deflate.Crc32
import net.maizegenetics.tassel.io.deflate.GZIP_FIXED_HEADER
import net.maizegenetics.tassel.io.deflate.GzipMemberHeader
import net.maizegenetics.tassel.io.deflate.InflateDecoder
import net.maizegenetics.tassel.io.deflate.gzipHeaderBytesRequired
import net.maizegenetics.tassel.io.deflate.inflateBgzf
import net.maizegenetics.tassel.io.deflate.parseBgzfPayload
import net.maizegenetics.tassel.io.deflate.parseGzipMemberHeader
import net.maizegenetics.tassel.io.deflate.truncatedBgzf
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report

/**
 * Seekable BGZF reader over a [RandomAccessSource] (SAMv1 section 4.1).
 *
 * A virtual offset is [bgzfVirtualOffset]: the block's compressed file offset in the high bits and the
 * uncompressed byte index in the low 16. [seek] loads that block. [read] and [readLine] then move forward
 * across blocks, and [tell] is the virtual offset of the next unread byte. After the last byte of a block,
 * [tell] is the next block's start, which is the offset htslib 1.25's `bgzf_getline` stores. A [seek] whose
 * low 16 bits equal the block's inflated length (at most 65535) is the same position and continues at the
 * next block. An empty block, including a terminal EOF block and one left in the middle by `cat`, is
 * checked and skipped.
 *
 * [read] fills a caller-owned array the way [Vcf.read] fills the window for [VcfRecordParser]: `0` means
 * there is no more inflated data. The parser stays source-free. A region read seeks here, fills that window,
 * and stops when [tell] reaches the chunk's end. This reader does not parse VCF.
 *
 * One decoder inflates blocks as they are needed. Parallel workers (#89) can replace that one call without
 * changing [seek] or [tell]. Recently inflated blocks are cached, because adjacent index chunks often share
 * a block, and a sequential run is read from the file in one batch rather than one read per block.
 *
 * [progress] is reported every few blocks, which is the suspension point a wasm read needs (ADR-0019).
 * A missing terminal EOF block is recorded on [warnings] with [MISSING_BGZF_EOF] and the bytes already
 * inflated are kept (ADR-0021). CRC32 and ISIZE are checked for every block.
 *
 * [close] does not close [source]. The caller owns the source and must keep it open while reading.
 * Not safe for concurrent calls.
 */
class BgzfReader(
    private val source: RandomAccessSource,
    private val progress: ProgressReporter = ProgressReporter.None,
) : AutoCloseable {
    private val recordedWarnings = ArrayList<String>(1)

    /** Empty until a read passes the end of the file without a terminal EOF block. */
    val warnings: List<String> get() = recordedWarnings

    private val decoder = InflateDecoder()
    private val crc = Crc32()
    private val scratch = ByteArray(BGZF_MAX_UNCOMPRESSED)
    private val ahead = ByteArray(READ_AHEAD)
    private val cache = BlockCache(CACHE_BLOCKS)
    private var lineBuffer = ByteArray(LINE_BUFFER)

    private var aheadFrom = 0L
    private var aheadTo = 0L
    private var blockAddress = 0L
    private var blockOffset = 0
    private var current: InflatedBlock? = null
    private var atEnd = false
    private var sawNonEmpty = false
    private var lastWasEmpty = false
    private var eofNoted = false
    private var blocksPassed = 0
    private var closed = false

    /**
     * Virtual offset of the next unread byte.
     *
     * At the end of a block this is the next block's start, not `uoffset == length`.
     */
    fun tell(): Long {
        ensureOpen()
        return bgzfVirtualOffset(blockAddress, blockOffset)
    }

    /**
     * Moves the next [read] to [virtualOffset].
     *
     * @throws IllegalArgumentException when the compressed offset is not a BGZF block start, or the
     *   uncompressed offset is past that block's length.
     */
    suspend fun seek(virtualOffset: Long) {
        ensureOpen()
        if (virtualOffset < 0) notABlockStart(virtualOffset)
        val placed = place(bgzfCoffset(virtualOffset), bgzfUoffset(virtualOffset), virtualOffset)
        atEnd = false
        when (placed) {
            is Place.End -> {
                blockAddress = placed.address
                blockOffset = 0
                current = null
                finishEnd()
            }
            is Place.At -> {
                current = placed.block
                blockAddress = placed.block.address
                blockOffset = placed.offset
            }
        }
    }

    /**
     * Copies up to [length] inflated bytes into [destination] at [destinationOffset].
     *
     * Returns the number of bytes copied. Zero means the inflated data has ended, including when [length]
     * is zero.
     */
    suspend fun read(destination: ByteArray, destinationOffset: Int, length: Int): Int {
        ensureOpen()
        requireReadArguments(0, destination, destinationOffset, length)
        if (length == 0) return 0
        var copied = 0
        while (copied < length) {
            if (!ensureCurrent()) return copied
            val block = current!!
            val count = minOf(length - copied, block.length - blockOffset)
            block.data.copyInto(destination, destinationOffset + copied, blockOffset, blockOffset + count)
            blockOffset += count
            copied += count
            if (blockOffset == block.length) snapToNext()
        }
        return copied
    }

    /**
     * Reads the next line, without its line break.
     *
     * `\n` ends a line and a preceding `\r` is removed. A line may cross a block boundary. Returns null
     * when there are no more bytes.
     */
    suspend fun readLine(): String? {
        ensureOpen()
        var size = 0
        while (true) {
            val value = nextByte()
            if (value < 0) {
                if (size == 0) return null
                return lineBuffer.decodeToString(0, size)
            }
            if (value == LINE_FEED) {
                val end = if (size > 0 && lineBuffer[size - 1] == CARRIAGE_RETURN) size - 1 else size
                return lineBuffer.decodeToString(0, end)
            }
            if (size == lineBuffer.size) lineBuffer = lineBuffer.copyOf(size * 2)
            lineBuffer[size] = value.toByte()
            size++
        }
    }

    override fun close() {
        closed = true
    }

    private suspend fun ensureCurrent(): Boolean {
        ensureOpen()
        if (atEnd) return false
        val block = current
        if (block != null && blockOffset < block.length) return true
        val virtual = bgzfVirtualOffset(blockAddress, blockOffset)
        return when (val placed = place(blockAddress, blockOffset, virtual)) {
            is Place.End -> {
                blockAddress = placed.address
                blockOffset = 0
                current = null
                finishEnd()
                false
            }
            is Place.At -> {
                current = placed.block
                blockAddress = placed.block.address
                blockOffset = placed.offset
                true
            }
        }
    }

    private suspend fun nextByte(): Int {
        if (!ensureCurrent()) return -1
        val block = current!!
        val value = block.data[blockOffset].toInt() and 0xFF
        blockOffset++
        if (blockOffset == block.length) snapToNext()
        return value
    }

    /** The block just finished. [tell] now names the next block, and the next read loads it. */
    private suspend fun snapToNext() {
        val next = current!!.nextAddress
        sawNonEmpty = true
        lastWasEmpty = false
        current = null
        blockAddress = next
        blockOffset = 0
        blocksPassed++
        if (blocksPassed % PROGRESS_BLOCKS == 0) {
            progress.report(blockAddress.coerceAtMost(source.size), source.size, "Reading BGZF")
        }
    }

    private suspend fun place(address: Long, uoffset: Int, virtualOffset: Long): Place {
        var addr = address
        var off = uoffset
        var virtual = virtualOffset
        while (true) {
            if (addr > source.size || (addr == source.size && off != 0)) notABlockStart(virtual)
            if (addr == source.size) return Place.End(addr)
            val block = blockAt(addr, virtual)
            if (block.length == 0) {
                if (off > 0) pastLength(virtual, 0)
                lastWasEmpty = true
                addr = block.nextAddress
                off = 0
                virtual = bgzfVirtualOffset(addr, 0)
                continue
            }
            sawNonEmpty = true
            lastWasEmpty = false
            if (off > block.length) pastLength(virtual, block.length)
            if (off == block.length) {
                addr = block.nextAddress
                off = 0
                virtual = bgzfVirtualOffset(addr, 0)
                continue
            }
            return Place.At(block, off)
        }
    }

    /**
     * Inflates the block at [address], or returns the cached copy.
     *
     * This is the one-worker decode step. #89 can run it on a worker per block; [seek] and [tell] stay.
     */
    private suspend fun blockAt(address: Long, virtualOffset: Long): InflatedBlock {
        cache[address]?.let { return it }
        val header = try {
            readHeader(address)
        } catch (_: IllegalArgumentException) {
            notABlockStart(virtualOffset)
        }
        if (header.blockSize < 0) notABlockStart(virtualOffset)
        val compressed = header.blockSize + 1
        val index = ensureBlock(address, compressed)
        val payload = parseBgzfPayload(ahead, index, header)
        val written = inflateBgzf(decoder, crc, ahead, payload, scratch)
        val data = if (written == 0) EMPTY else scratch.copyOf(written)
        val block = InflatedBlock(address, data, written, address + compressed)
        cache.put(block)
        return block
    }

    private suspend fun readHeader(address: Long): GzipMemberHeader {
        if (source.size - address < GZIP_FIXED_HEADER) truncatedBgzf("gzip header")
        var needed = GZIP_FIXED_HEADER
        while (true) {
            val index = ensureBlock(address, needed)
            val available = (aheadTo - address).toInt()
            val required = gzipHeaderBytesRequired(ahead, index, available)
            if (required in 0..BGZF_MAX_BLOCK) return parseGzipMemberHeader(ahead, index, available)
            val want = if (required < 0) -required else required
            if (want > BGZF_MAX_BLOCK || address > source.size - want) truncatedBgzf("gzip extra field")
            needed = want
        }
    }

    /**
     * Makes sure [length] bytes at [position] sit in [ahead]. Returns the index of [position].
     *
     * A miss reads one batch of [READ_AHEAD] bytes. A block that runs past the batch but still fits in
     * the buffer is finished with one tail read.
     */
    private suspend fun ensureBlock(position: Long, length: Int): Int {
        val end = position + length
        if (position >= aheadFrom && end <= aheadTo) return (position - aheadFrom).toInt()
        if (position >= aheadFrom && position < aheadTo) {
            val index = (position - aheadFrom).toInt()
            if (index + length <= ahead.size) {
                val missing = (end - aheadTo).toInt()
                if (missing > 0) {
                    val dest = (aheadTo - aheadFrom).toInt()
                    aheadTo += source.readAt(aheadTo, ahead, dest, missing)
                }
                if (aheadTo < end) truncatedBgzf("BGZF block")
                return index
            }
        }
        if (position >= source.size || source.size - position < length) truncatedBgzf("BGZF block")
        val want = minOf(ahead.size.toLong(), source.size - position).toInt()
        val got = source.readAt(position, ahead, 0, want)
        if (got < length) truncatedBgzf("BGZF block")
        aheadFrom = position
        aheadTo = position + got
        return 0
    }

    private fun finishEnd() {
        atEnd = true
        current = null
        blockOffset = 0
        if (eofNoted) return
        eofNoted = true
        if (sawNonEmpty && !lastWasEmpty) recordedWarnings.add(MISSING_BGZF_EOF)
    }

    private fun ensureOpen() {
        if (closed) throw IllegalStateException("BGZF reader is closed")
    }

    private class InflatedBlock(val address: Long, val data: ByteArray, val length: Int, val nextAddress: Long)

    private class BlockCache(private val capacity: Int) {
        private val blocks = HashMap<Long, InflatedBlock>(capacity)
        private val order = ArrayDeque<Long>()

        operator fun get(address: Long): InflatedBlock? {
            val block = blocks[address] ?: return null
            order.remove(address)
            order.addLast(address)
            return block
        }

        fun put(block: InflatedBlock) {
            if (blocks.size == capacity) {
                blocks.remove(order.removeAt(0))
            }
            blocks[block.address] = block
            order.addLast(block.address)
        }
    }

    private sealed class Place {
        class At(val block: InflatedBlock, val offset: Int) : Place()
        class End(val address: Long) : Place()
    }

    private companion object {
        const val READ_AHEAD = 256 * 1024
        const val CACHE_BLOCKS = 8
        const val PROGRESS_BLOCKS = 4
        const val LINE_BUFFER = 1024
        const val LINE_FEED = '\n'.code
        const val CARRIAGE_RETURN = '\r'.code.toByte()
        val EMPTY = ByteArray(0)
    }
}

/** Packs a BGZF virtual offset. [uoffset] is the byte index inside the inflated block. */
internal fun bgzfVirtualOffset(coffset: Long, uoffset: Int): Long {
    require(coffset >= 0) { "coffset is negative: $coffset" }
    require(uoffset in 0..0xFFFF) { "uoffset is outside 0..65535: $uoffset" }
    return (coffset shl 16) or uoffset.toLong()
}

/** Compressed file offset stored in a BGZF virtual offset. */
internal fun bgzfCoffset(virtualOffset: Long): Long = virtualOffset ushr 16

/** Uncompressed byte index stored in a BGZF virtual offset. */
internal fun bgzfUoffset(virtualOffset: Long): Int = (virtualOffset and 0xFFFF).toInt()

private fun notABlockStart(virtualOffset: Long): Nothing =
    throw IllegalArgumentException("BGZF virtual offset $virtualOffset is not a block start")

private fun pastLength(virtualOffset: Long, length: Int): Nothing =
    throw IllegalArgumentException("BGZF virtual offset $virtualOffset is past the block length $length")
