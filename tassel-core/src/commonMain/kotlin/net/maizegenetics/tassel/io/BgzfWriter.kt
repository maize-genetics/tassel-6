package net.maizegenetics.tassel.io

import kotlinx.io.Sink
import net.maizegenetics.tassel.io.deflate.BGZF_BLOCK_SIZE
import net.maizegenetics.tassel.io.deflate.BGZF_EOF_BLOCK
import net.maizegenetics.tassel.io.deflate.Crc32
import net.maizegenetics.tassel.io.deflate.DeflateCompressor
import net.maizegenetics.tassel.io.deflate.writeBgzfMember

/**
 * BGZF writer (SAMv1 section 4.1, ADR-0022).
 *
 * [layout] chooses the block split. [BgzfLayout.Vcf] is htslib 1.25 `bgzip` text mode: a 65280-byte
 * window, the header (`#` and `@` lines) flushed through its last complete line, then each later block
 * ending at the last newline in the window. A line longer than the window is cut, which is the only
 * block that does not end on a newline. A file that fits in one window is a header block plus one body
 * block. [BgzfLayout.Bcf] packs 65280 uncompressed bytes and may split a record. [endHeader] ends the
 * header block.
 *
 * The default [compressor] is one fixed DEFLATE level. A null payload, or one that would not fit in 65536
 * bytes, is stored instead. [close] appends the 28-byte EOF member and leaves [sink] open. [write] does not
 * suspend. [flush] pushes bytes already framed and does not end a block. [tell] is the virtual offset of
 * the next byte still held in the window: the compressed offset of that open block, and how many
 * uncompressed bytes are buffered. A VCF scan can move a held tail into the next block, so a mid-window
 * [tell] is not an index offset. This epic does not write an index.
 *
 * Not safe for concurrent calls. The same input and layout produce the same bytes on every target.
 *
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
 */
class BgzfWriter(
    private val sink: Sink,
    private val layout: BgzfLayout,
    private val compressor: BgzfCompressor = BgzfCompressor.deflate(),
) : AutoCloseable {
    private val window = ByteArray(BGZF_BLOCK_SIZE)
    private val crc = Crc32()
    private var length = 0
    private var blockAddress = 0L
    private var inHeader = layout == BgzfLayout.Vcf
    private var longLine = false
    private var needsScan = false
    private var closed = false

    /**
     * Accepts `source[startIndex, endIndex)`.
     *
     * How the bytes are sliced across calls does not change the blocks.
     */
    fun write(source: ByteArray, startIndex: Int = 0, endIndex: Int = source.size) {
        ensureOpen()
        if (startIndex < 0 || endIndex < startIndex || endIndex > source.size) {
            throw IllegalArgumentException("invalid write range")
        }
        var index = startIndex
        while (index < endIndex) {
            val space = BGZF_BLOCK_SIZE - length
            if (space == 0) throw IllegalStateException("BGZF window was not consumed")
            val take = minOf(space, endIndex - index)
            source.copyInto(window, length, index, index + take)
            length += take
            index += take
            when (layout) {
                BgzfLayout.Bcf -> if (length == BGZF_BLOCK_SIZE) emitOpenBlock()
                BgzfLayout.Vcf -> {
                    needsScan = true
                    if (length == BGZF_BLOCK_SIZE) scanWindow()
                }
            }
        }
    }

    /**
     * Virtual offset of the next uncompressed byte, packed the way [bgzfVirtualOffset] packs a read.
     *
     * The high bits are the compressed offset where the open block will start. The low 16 are the
     * number of uncompressed bytes already accepted into that block.
     */
    fun tell(): Long {
        ensureOpen()
        return bgzfVirtualOffset(blockAddress, length)
    }

    /** Pushes framed bytes to [sink]. Bytes still in the window stay there, and [tell] does not move. */
    fun flush() {
        ensureOpen()
        sink.flush()
    }

    /**
     * Ends the header block.
     *
     * For [BgzfLayout.Bcf], the bytes accepted since the previous block become that block. A header
     * longer than 65280 bytes was already split on that boundary. For [BgzfLayout.Vcf], call this
     * only before any body byte: the buffered header is emitted and later bytes use the body rule.
     * After the header has already ended, this does not cut a body line.
     */
    fun endHeader() {
        ensureOpen()
        when (layout) {
            BgzfLayout.Bcf -> if (length > 0) emitOpenBlock()
            BgzfLayout.Vcf -> {
                if (!inHeader) return
                if (length > 0) emitOpenBlock()
                inHeader = false
                longLine = false
                needsScan = false
            }
        }
    }

    /** Writes the 28-byte EOF member. [sink] stays open. A second call writes nothing. */
    override fun close() {
        if (closed) return
        when (layout) {
            BgzfLayout.Vcf -> {
                if (needsScan && length > 0) scanWindow()
                if (length > 0) emitOpenBlock()
            }
            BgzfLayout.Bcf -> if (length > 0) emitOpenBlock()
        }
        sink.write(BGZF_EOF_BLOCK)
        blockAddress += BGZF_EOF_BLOCK.size
        closed = true
    }

    /**
     * One `bgzip` text-mode pass over the bytes currently in the window.
     *
     * A full window is the read `bgzip` issues. A short window is the last read, or the tail `close`
     * still has to place. Bytes after the chosen newline stay in the window for the next pass.
     */
    private fun scanWindow() {
        val available = length
        var flush = false
        var consumed = available
        if (inHeader && (longLine || (available > 0 && (window[0] == HASH || window[0] == AT)))) {
            var lastStart = 0
            var index = 0
            while (index < available) {
                val newline = window[index] == NEWLINE
                index++
                if (!newline) continue
                lastStart = index
                if (index < available && window[index] != HASH && window[index] != AT) {
                    inHeader = false
                    break
                }
            }
            if (lastStart == 0) {
                longLine = true
            } else {
                consumed = lastStart
                flush = true
                longLine = false
            }
        } else {
            var index = available - 1
            while (index >= 0 && window[index] != NEWLINE) index--
            if (index >= 0) {
                consumed = index + 1
                flush = true
            }
        }
        if (flush || consumed == BGZF_BLOCK_SIZE) {
            emit(consumed)
            val tail = available - consumed
            for (offset in 0 until tail) window[offset] = window[consumed + offset]
            length = tail
        }
        needsScan = false
    }

    private fun emitOpenBlock() {
        emit(length)
        length = 0
    }

    private fun emit(count: Int) {
        if (count == 0) return
        val deflate = compressor.compress(window, 0, count)
        blockAddress += sink.writeBgzfMember(window, 0, count, deflate, crc)
    }

    private fun ensureOpen() {
        if (closed) throw IllegalStateException("BGZF writer is closed")
    }

    private companion object {
        val NEWLINE = '\n'.code.toByte()
        val HASH = '#'.code.toByte()
        val AT = '@'.code.toByte()
    }
}

/** Which htslib block split [BgzfWriter] applies. */
enum class BgzfLayout {
    /** VCF text, split the way `bgzip` 1.25 splits a text file. */
    Vcf,

    /** BCF records. A block holds 65280 uncompressed bytes, and a record may cross it. */
    Bcf,
}

/**
 * Turns one uncompressed BGZF block into raw DEFLATE, or declines and leaves the block stored.
 *
 * [input] is only valid until this call returns. Null, and a payload that would make the member
 * larger than 65536 bytes, both store the block. The split is chosen before this is called.
 */
fun interface BgzfCompressor {
    fun compress(input: ByteArray, from: Int, to: Int): ByteArray?

    companion object {
        /** Every data block is stored. */
        val Stored: BgzfCompressor = BgzfCompressor { _, _, _ -> null }

        /**
         * One fixed DEFLATE level.
         *
         * Each call returns a compressor with its own scratch buffers. [BgzfCompressor.compress] on that
         * instance is not safe to call from more than one worker.
         */
        fun deflate(): BgzfCompressor = DeflateCompressor()
    }
}
