package net.maizegenetics.tassel.app.export

import kotlinx.io.Sink
import kotlinx.io.writeString

/**
 * Where an exported result goes on this target: a save dialog and a file on desktop, a download in the browser.
 *
 * It deals in a `kotlinx.io.Sink` rather than text because that is what the core writers take between them:
 * `ResultTable`, `PhenotypeTable`, and `DistanceMatrix` write to an `Appendable` (wrap the sink in
 * [SinkAppendable]), and `HapMap.write` writes to the sink directly.
 */
interface ResultFiles {
    /**
     * Saves one file called [fileName], filled in by [write].
     *
     * @return false when the user cancelled, so the caller can say nothing rather than claim a save.
     */
    suspend fun save(fileName: String, write: suspend (Sink) -> Unit): Boolean
}

/**
 * Lets the core's `Appendable` writers stream into a `kotlinx.io.Sink`.
 *
 * Characters are collected in a [StringBuilder] and handed over a chunk at a time, so a 281 x 281 kinship matrix
 * costs one `String` per [chunkSize] characters rather than one per character. `tassel-cli`'s own `SinkAppendable`
 * does the latter and is a queued performance finding on tassel-6#24; this is the shape that finding asks for.
 *
 * Call [flush] when the writer is done. Nothing here is thread safe, and nothing needs to be: one writer fills it.
 *
 * @param chunkSize characters to buffer before writing them through.
 */
class SinkAppendable(private val sink: Sink, private val chunkSize: Int = DEFAULT_CHUNK) : Appendable {
    private val pending = StringBuilder(chunkSize)

    /** Strings handed to the sink so far, which is what makes the batching testable. */
    var writes: Int = 0
        private set

    override fun append(value: Char): Appendable = apply {
        pending.append(value)
        flushIfFull()
    }

    override fun append(value: CharSequence?): Appendable = apply {
        pending.append(value)
        flushIfFull()
    }

    override fun append(value: CharSequence?, startIndex: Int, endIndex: Int): Appendable = apply {
        pending.append(value, startIndex, endIndex)
        flushIfFull()
    }

    /** Writes whatever is buffered. */
    fun flush() {
        if (pending.isEmpty()) return
        sink.writeString(pending.toString())
        writes++
        pending.clear()
    }

    private fun flushIfFull() {
        if (pending.length >= chunkSize) flush()
    }

    companion object {
        /** Characters buffered before a write. Large enough that a matrix row never costs more than one. */
        const val DEFAULT_CHUNK: Int = 8192
    }
}

/** Runs [write] against this sink as an [Appendable] and flushes what is left. */
inline fun Sink.appendText(write: (Appendable) -> Unit) {
    val out = SinkAppendable(this)
    write(out)
    out.flush()
}
