package net.maizegenetics.tassel.io.deflate

/**
 * LSB-first bit packer for one raw DEFLATE block.
 *
 * Huffman codes are written most-significant bit first, which is how RFC 1951 packs them into that stream.
 * The same order is what [InflateDecoder] reads. A partial last byte is padded with zero bits.
 */
internal class DeflateBitSink {
    private var out: ByteArray = EMPTY
    private var index = 0
    private var acc = 0
    private var bitCount = 0
    private var overflow = false

    fun reset(target: ByteArray) {
        out = target
        index = 0
        acc = 0
        bitCount = 0
        overflow = false
    }

    /** Writes the low [length] bits of [value], low bit first. */
    fun writeLsb(value: Int, length: Int) {
        var remaining = value
        var left = length
        while (left > 0) {
            if (overflow) return
            val space = 8 - bitCount
            val take = if (left < space) left else space
            acc = acc or ((remaining and ((1 shl take) - 1)) shl bitCount)
            bitCount += take
            remaining = remaining ushr take
            left -= take
            if (bitCount == 8) flushByte()
        }
    }

    /** Writes the low [length] bits of [value], high bit first. */
    fun writeMsb(value: Int, length: Int) {
        var mask = 1 shl (length - 1)
        repeat(length) {
            writeLsb(if (value and mask != 0) 1 else 0, 1)
            mask = mask ushr 1
        }
    }

    /** Bytes written, or -1 when [out] was too small. */
    fun finish(): Int {
        if (!overflow && bitCount > 0) flushByte()
        return if (overflow) -1 else index
    }

    private fun flushByte() {
        if (index >= out.size) {
            overflow = true
            acc = 0
            bitCount = 0
            return
        }
        out[index++] = acc.toByte()
        acc = 0
        bitCount = 0
    }

    private companion object {
        val EMPTY = ByteArray(0)
    }
}
