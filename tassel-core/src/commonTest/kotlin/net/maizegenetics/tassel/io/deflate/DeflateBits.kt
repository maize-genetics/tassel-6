package net.maizegenetics.tassel.io.deflate

/**
 * Test-only bit packer. Raw DEFLATE fields are written least-significant bit first. Huffman codes are
 * written most-significant bit first, which is how RFC 1951 packs them into that LSB-first stream.
 */
internal class BitWriter {
    private val out = ArrayList<Byte>(64)
    private var acc = 0
    private var nbits = 0

    fun writeLsb(value: Int, length: Int) {
        var remaining = value
        repeat(length) {
            writeBit(remaining and 1)
            remaining = remaining ushr 1
        }
    }

    fun writeMsb(value: Int, length: Int) {
        var mask = 1 shl (length - 1)
        repeat(length) {
            writeBit(if (value and mask != 0) 1 else 0)
            mask = mask ushr 1
        }
    }

    fun align() {
        if (nbits != 0) flushByte()
    }

    fun writeByte(value: Int) {
        check(nbits == 0) { "stored bytes must follow a byte boundary" }
        out.add(value.toByte())
    }

    fun writeShortLe(value: Int) {
        writeByte(value and 0xFF)
        writeByte((value ushr 8) and 0xFF)
    }

    fun toByteArray(): ByteArray {
        if (nbits != 0) flushByte()
        return ByteArray(out.size) { out[it] }
    }

    private fun writeBit(bit: Int) {
        acc = acc or ((bit and 1) shl nbits)
        nbits++
        if (nbits == 8) flushByte()
    }

    private fun flushByte() {
        out.add(acc.toByte())
        acc = 0
        nbits = 0
    }
}

internal fun fixedLiteralCode(symbol: Int): Pair<Int, Int> = when (symbol) {
    in 0..143 -> (symbol + 48) to 8
    in 144..255 -> (symbol - 144 + 400) to 9
    in 256..279 -> (symbol - 256) to 7
    in 280..287 -> (symbol - 280 + 192) to 8
    else -> error("symbol $symbol has no fixed literal code")
}

internal fun storedBlock(data: ByteArray, last: Boolean): ByteArray {
    val writer = BitWriter()
    writeStored(writer, data, last)
    return writer.toByteArray()
}

internal fun writeStored(writer: BitWriter, data: ByteArray, last: Boolean) {
    writer.writeLsb(if (last) 1 else 0, 1)
    writer.writeLsb(0, 2)
    writer.align()
    writer.writeShortLe(data.size)
    writer.writeShortLe(data.size xor 0xFFFF)
    for (byte in data) writer.writeByte(byte.toInt() and 0xFF)
}

internal fun fixedBlock(literals: ByteArray, last: Boolean, match: Match? = null): ByteArray {
    val writer = BitWriter()
    writeFixed(writer, literals, last, match)
    return writer.toByteArray()
}

internal fun writeFixed(writer: BitWriter, literals: ByteArray, last: Boolean, match: Match? = null) {
    writer.writeLsb(if (last) 1 else 0, 1)
    writer.writeLsb(1, 2)
    writeLiterals(writer, literals)
    if (match != null) writeFixedMatch(writer, match.length, match.distance)
    val (endCode, endLength) = fixedLiteralCode(256)
    writer.writeMsb(endCode, endLength)
}

internal data class Match(val length: Int, val distance: Int)

private fun writeLiterals(writer: BitWriter, literals: ByteArray) {
    for (byte in literals) {
        val (code, length) = fixedLiteralCode(byte.toInt() and 0xFF)
        writer.writeMsb(code, length)
    }
}

internal fun writeFixedMatch(writer: BitWriter, length: Int, distance: Int) {
    val lengthIndex = lengthIndex(length)
    val (lengthCode, lengthBits) = fixedLiteralCode(257 + lengthIndex)
    writer.writeMsb(lengthCode, lengthBits)
    val lengthExtra = LENGTH_EXTRA[lengthIndex]
    if (lengthExtra != 0) writer.writeLsb(length - LENGTH_BASE[lengthIndex], lengthExtra)
    val distanceIndex = distanceIndex(distance)
    writer.writeMsb(distanceIndex, 5)
    val distanceExtra = DISTANCE_EXTRA[distanceIndex]
    if (distanceExtra != 0) writer.writeLsb(distance - DISTANCE_BASE[distanceIndex], distanceExtra)
}

private fun lengthIndex(length: Int): Int {
    for (index in LENGTH_BASE.indices) {
        val span = (1 shl LENGTH_EXTRA[index]) - 1
        if (length in LENGTH_BASE[index]..LENGTH_BASE[index] + span) return index
    }
    error("length $length is outside 3..258")
}

private fun distanceIndex(distance: Int): Int {
    for (index in DISTANCE_BASE.indices) {
        val span = (1 shl DISTANCE_EXTRA[index]) - 1
        if (distance in DISTANCE_BASE[index]..DISTANCE_BASE[index] + span) return index
    }
    error("distance $distance is outside 1..32768")
}

private val LENGTH_BASE = intArrayOf(
    3, 4, 5, 6, 7, 8, 9, 10,
    11, 13, 15, 17, 19, 23, 27, 31,
    35, 43, 51, 59, 67, 83, 99, 115,
    131, 163, 195, 227, 258,
)
private val LENGTH_EXTRA = intArrayOf(
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 2, 2, 2, 2,
    3, 3, 3, 3, 4, 4, 4, 4,
    5, 5, 5, 5, 0,
)
private val DISTANCE_BASE = intArrayOf(
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
    257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
)
private val DISTANCE_EXTRA = intArrayOf(
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
    7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13,
)
