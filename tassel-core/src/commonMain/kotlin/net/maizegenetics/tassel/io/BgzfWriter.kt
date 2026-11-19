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
 * The no-compressor constructor uses one fixed DEFLATE level on [bgzfWorkerCount] workers. One worker,
 * which is the count on wasmJs, compresses on the caller. More workers compress on
 * [kotlinx.coroutines.Dispatchers.Default], and the caller writes the members in block order. At most
 * twice the worker count are copied and not yet framed. A caller-supplied [BgzfCompressor] stays on the
 * caller, because one instance is not safe to share. A null payload, or one that would not fit in 65536
 * bytes, is stored instead. [close] appends the 28-byte EOF member and leaves [sink] open. [write] does not
 * suspend. A block [write] has closed may still be compressing when [write] returns; [tell], [flush], and
 * [close] finish those blocks in order before they observe the sink. [flush] pushes bytes already framed
 * and does not end a block. [tell] is the virtual offset of the next byte still held in the window: the
 * compressed offset of that open block, and how many uncompressed bytes are buffered. A VCF scan can move
 * a held tail into the next block, so a mid-window [tell] is not an index offset. This epic does not write
 * an index.
 *
 * Not safe for concurrent calls. The same input and layout produce the same bytes on every target and for
 * any worker count.
 *
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
 */
class BgzfWriter internal constructor(
    private val sink: Sink,
    private val layout: BgzfLayout,
    private val workers: Int,
    private val maxBlocksInFlight: Int,
    private val newCompressor: () -> BgzfCompressor,
    private val beforeCompress: suspend (Int) -> Unit,
) : AutoCloseable {
    /** Fixed DEFLATE level, on [bgzfWorkerCount] workers. */
    constructor(sink: Sink, layout: BgzfLayout) : this(
        sink,
        layout,
        bgzfWorkerCount(),
        bgzfMaxBlocksInFlight(bgzfWorkerCount()),
        { BgzfCompressor.deflate() },
        { _ -> },
    )

    /**
     * Compresses every block with [compressor], on the caller.
     *
     * One instance is not shared across workers. Pass no compressor for the parallel default.
     */
    constructor(sink: Sink, layout: BgzfLayout, compressor: BgzfCompressor) : this(
        sink,
        layout,
        1,
        1,
        { compressor },
        { _ -> },
    )

    /**
     * Fixed DEFLATE level on [workers] threads.
     *
     * Tests and [openBgzfWriter] use this. The public constructors take no worker count.
     */
    internal constructor(
        sink: Sink,
        layout: BgzfLayout,
        workers: Int,
        beforeCompress: suspend (Int) -> Unit = { _ -> },
    ) : this(
        sink,
        layout,
        workers.coerceAtLeast(1),
        bgzfMaxBlocksInFlight(workers.coerceAtLeast(1)),
        { BgzfCompressor.deflate() },
        beforeCompress,
    )

    private val window = ByteArray(BGZF_BLOCK_SIZE)
    private val crc = Crc32()
    private val inlineCompressor: BgzfCompressor? = if (workers == 1) newCompressor() else null
    private var pipeline: BgzfEncodeWindow? = null
    private var length = 0
    private var blockAddress = 0L
    private var inHeader = layout == BgzfLayout.Vcf
    private var longLine = false
    private var needsScan = false
    private var closed = false
    private var failed = false

    init {
        require(workers >= 1) { "workers must be positive: $workers" }
        require(maxBlocksInFlight >= 1) { "blocks in flight must be positive: $maxBlocksInFlight" }
    }

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
        drainPipeline()
        return bgzfVirtualOffset(blockAddress, length)
    }

    /** Pushes framed bytes to [sink]. Bytes still in the window stay there, and [tell] does not move. */
    fun flush() {
        ensureOpen()
        drainPipeline()
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
        if (closed || failed) return
        try {
            pipeline?.takeFailure()?.let { throw it }
            when (layout) {
                BgzfLayout.Vcf -> {
                    if (needsScan && length > 0) scanWindow()
                    if (length > 0) emitOpenBlock()
                }
                BgzfLayout.Bcf -> if (length > 0) emitOpenBlock()
            }
            drainPipeline()
            pipeline?.takeFailure()?.let { throw it }
            sink.write(BGZF_EOF_BLOCK)
            blockAddress += BGZF_EOF_BLOCK.size
            closed = true
        } catch (e: Throwable) {
            failed = true
            throw e
        } finally {
            pipeline?.shutdown()
        }
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
        val inline = inlineCompressor
        if (inline != null) {
            val deflate = inline.compress(window, 0, count)
            blockAddress += sink.writeBgzfMember(window, 0, count, deflate, crc)
            return
        }
        // Copy before the caller slides a held tail over [window].
        val copy = window.copyOf(count)
        val encode = ensurePipeline()
        if (!encode.trySubmit(copy)) {
            bgzfBlocking {
                encode.submit(copy) { input, deflate ->
                    blockAddress += sink.writeBgzfMember(input, 0, input.size, deflate, crc)
                }
            }
        }
    }

    private fun ensurePipeline(): BgzfEncodeWindow =
        pipeline ?: BgzfEncodeWindow(workers, maxBlocksInFlight, newCompressor, beforeCompress).also { pipeline = it }

    private fun drainPipeline() {
        val encode = pipeline ?: return
        if (encode.blocksInFlight == 0) return
        bgzfBlocking {
            encode.drain { input, deflate ->
                blockAddress += sink.writeBgzfMember(input, 0, input.size, deflate, crc)
            }
        }
    }

    private fun ensureOpen() {
        if (closed || failed) throw IllegalStateException("BGZF writer is closed")
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

/**
 * A [BgzfWriter] that compresses on [workers] threads.
 *
 * [BgzfWriter]'s public constructors take no worker count, and [Bcf.write] uses those. This is the bench
 * seam. Bytes do not depend on [workers].
 */
fun openBgzfWriter(sink: Sink, layout: BgzfLayout, workers: Int): BgzfWriter =
    BgzfWriter(sink, layout, workers.coerceAtLeast(1))
