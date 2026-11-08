package net.maizegenetics.tassel.io.deflate

/**
 * One-shot raw DEFLATE decoder (RFC 1951).
 *
 * [inflate] decodes a complete stream from `input[from, to)` into `output[at, at + capacity)` and returns
 * the number of bytes written. The same instance is reused across blocks: the fixed tables are built once,
 * and each dynamic block overwrites the scratch tables. [inflate] does not allocate.
 *
 * Not safe for concurrent calls on one instance. A failed call can leave the scratch tables half-built;
 * the next [inflate] rebuilds whatever it needs.
 *
 * Gzip and BGZF framing stay outside this type. Streaming, with a 32 KiB window, is a later variant.
 */
internal class InflateDecoder {
    private val fixedLiteral = HuffmanTable(
        HuffmanTable.LITERAL_ENOUGH,
        HuffmanTable.LITERAL_TABLE_BITS,
        "invalid length symbol",
    )
    private val fixedDistance = HuffmanTable(
        HuffmanTable.DISTANCE_ENOUGH,
        HuffmanTable.DISTANCE_TABLE_BITS,
        "invalid distance symbol",
    )
    private val literal = HuffmanTable(
        HuffmanTable.LITERAL_ENOUGH,
        HuffmanTable.LITERAL_TABLE_BITS,
        "invalid length symbol",
    )
    private val distance = HuffmanTable(
        HuffmanTable.DISTANCE_ENOUGH,
        HuffmanTable.DISTANCE_TABLE_BITS,
        "invalid distance symbol",
    )
    private val precode = HuffmanTable(
        HuffmanTable.PRECODE_ENOUGH,
        HuffmanTable.PRECODE_TABLE_BITS,
        "invalid length symbol",
    )
    private val precodeLengths = IntArray(PRECODE_SYMBOLS)
    private val codeLengths = IntArray(MAX_LITERAL_CODES + MAX_DISTANCE_CODES)

    private var input: ByteArray = EMPTY
    private var inputPos: Int = 0
    private var inputEnd: Int = 0
    private var output: ByteArray = EMPTY
    private var outputPos: Int = 0
    private var outputEnd: Int = 0
    private var origin: Int = 0
    private var bits: Long = 0
    private var bitCount: Int = 0

    init {
        val literalLengths = IntArray(MAX_LITERAL_CODES)
        for (symbol in 0..143) literalLengths[symbol] = 8
        for (symbol in 144..255) literalLengths[symbol] = 9
        for (symbol in 256..279) literalLengths[symbol] = 7
        for (symbol in 280 until MAX_LITERAL_CODES) literalLengths[symbol] = 8
        fixedLiteral.build(literalLengths, 0, MAX_LITERAL_CODES, allowEmpty = false)
        fixedDistance.build(IntArray(MAX_DISTANCE_CODES) { 5 }, 0, MAX_DISTANCE_CODES, allowEmpty = false)
    }

    /**
     * Decodes one raw DEFLATE stream.
     *
     * Padding bits in the last consumed byte are ignored. Any byte after that is an error.
     *
     * @throws IllegalArgumentException on a bad range or a stream that is not valid raw DEFLATE.
     */
    fun inflate(input: ByteArray, from: Int, to: Int, output: ByteArray, at: Int, capacity: Int): Int {
        if (from < 0 ||
            to < from ||
            to > input.size ||
            at < 0 ||
            capacity < 0 ||
            at > output.size ||
            capacity > output.size - at
        ) {
            fail("invalid inflate range")
        }
        this.input = input
        inputPos = from
        inputEnd = to
        this.output = output
        outputPos = at
        outputEnd = at + capacity
        origin = at
        bits = 0
        bitCount = 0
        try {
            while (true) {
                val last = take(1) == 1
                when (val type = take(2)) {
                    0 -> decodeStored()
                    1 -> decodeHuffman(fixedLiteral, fixedDistance)
                    2 -> decodeDynamic()
                    else -> fail("invalid block type")
                }
                if (last) break
            }
            if (inputPos - bitCount / 8 != inputEnd) fail("trailing input")
            return outputPos - origin
        } finally {
            this.input = EMPTY
            this.output = EMPTY
        }
    }

    private fun decodeStored() {
        alignToByte()
        if (inputEnd - inputPos < 4) fail("truncated input")
        val length = readU16(inputPos)
        val complement = readU16(inputPos + 2)
        inputPos += 4
        if ((length xor complement) != 0xFFFF) {
            fail("stored block length does not match its complement")
        }
        if (inputEnd - inputPos < length) fail("truncated input")
        if (outputEnd - outputPos < length) fail("output overflow")
        input.copyInto(output, outputPos, inputPos, inputPos + length)
        inputPos += length
        outputPos += length
    }

    private fun decodeDynamic() {
        val literalCount = take(5) + 257
        val distanceCount = take(5) + 1
        val precodeCount = take(4) + 4
        precodeLengths.fill(0)
        for (index in 0 until precodeCount) {
            precodeLengths[CODE_LENGTH_ORDER[index]] = take(3)
        }
        precode.build(precodeLengths, 0, PRECODE_SYMBOLS, allowEmpty = false)
        val total = literalCount + distanceCount
        var index = 0
        while (index < total) {
            when (val symbol = symbol(precode)) {
                in 0..15 -> codeLengths[index++] = symbol
                16 -> {
                    if (index == 0) fail("invalid length symbol")
                    val repeat = take(2) + 3
                    if (index + repeat > total) fail("invalid length symbol")
                    val previous = codeLengths[index - 1]
                    repeat(repeat) { codeLengths[index++] = previous }
                }
                17 -> {
                    val repeat = take(3) + 3
                    if (index + repeat > total) fail("invalid length symbol")
                    repeat(repeat) { codeLengths[index++] = 0 }
                }
                18 -> {
                    val repeat = take(7) + 11
                    if (index + repeat > total) fail("invalid length symbol")
                    repeat(repeat) { codeLengths[index++] = 0 }
                }
                else -> fail("invalid length symbol")
            }
        }
        literal.build(codeLengths, 0, literalCount, allowEmpty = false)
        distance.build(codeLengths, literalCount, distanceCount, allowEmpty = true)
        decodeHuffman(literal, distance)
    }

    private fun decodeHuffman(literalTable: HuffmanTable, distanceTable: HuffmanTable) {
        while (inputEnd - inputPos >= FAST_INPUT_BYTES && outputEnd - outputPos >= MAX_MATCH) {
            when (val symbol = symbol(literalTable)) {
                in 0..255 -> output[outputPos++] = symbol.toByte()
                END_OF_BLOCK -> return
                in 257..285 -> writeMatch(symbol, distanceTable, checked = false)
                else -> fail("invalid length symbol")
            }
        }
        while (true) {
            when (val symbol = symbol(literalTable)) {
                in 0..255 -> {
                    if (outputPos == outputEnd) fail("output overflow")
                    output[outputPos++] = symbol.toByte()
                }
                END_OF_BLOCK -> return
                in 257..285 -> writeMatch(symbol, distanceTable, checked = true)
                else -> fail("invalid length symbol")
            }
        }
    }

    private fun writeMatch(symbol: Int, distanceTable: HuffmanTable, checked: Boolean) {
        val length = LENGTH_BASE[symbol - 257] + extra(LENGTH_EXTRA[symbol - 257])
        if (distanceTable.empty) fail("invalid distance symbol")
        val distanceSymbol = symbol(distanceTable)
        if (distanceSymbol > 29) fail("invalid distance symbol")
        val distance = DISTANCE_BASE[distanceSymbol] + extra(DISTANCE_EXTRA[distanceSymbol])
        if (distance > outputPos - origin) fail("distance past the start of the output")
        if (checked && length > outputEnd - outputPos) fail("output overflow")
        copyMatch(distance, length)
    }

    private fun copyMatch(distance: Int, length: Int) {
        val dest = outputPos
        val src = dest - distance
        var index = 0
        if (distance >= 8) {
            while (index + 8 <= length) {
                val from = src + index
                output[dest + index] = output[from]
                output[dest + index + 1] = output[from + 1]
                output[dest + index + 2] = output[from + 2]
                output[dest + index + 3] = output[from + 3]
                output[dest + index + 4] = output[from + 4]
                output[dest + index + 5] = output[from + 5]
                output[dest + index + 6] = output[from + 6]
                output[dest + index + 7] = output[from + 7]
                index += 8
            }
        }
        while (index < length) {
            output[dest + index] = output[src + index]
            index++
        }
        outputPos = dest + length
    }

    /**
     * Next symbol. Uses the flat table when 15 bits are already buffered or two more input bytes can
     * supply them, and walks the canonical codes one bit at a time at the end of the buffer.
     */
    private fun symbol(table: HuffmanTable): Int {
        if (bitCount >= 15 || inputEnd - inputPos >= 2) {
            ensure(15)
            val entry = table.lookup((bits and 0x7FFF).toInt())
            if (entry and HuffmanTable.FLAG_INVALID != 0) fail(table.invalidMessage)
            val length = (entry ushr HuffmanTable.LEN_SHIFT) and 0xF
            if (length == 0 || length > bitCount) fail("truncated input")
            bits = bits ushr length
            bitCount -= length
            return entry and HuffmanTable.PAYLOAD_MASK
        }
        return symbolSlow(table)
    }

    private fun symbolSlow(table: HuffmanTable): Int {
        var code = 0
        var first = 0
        var index = 0
        for (length in 1..HuffmanTable.MAX_CODE_LENGTH) {
            code = code or nextBit()
            val count = table.counts[length]
            if (code - count < first) return table.symbols[index + (code - first)]
            index += count
            first += count
            first = first shl 1
            code = code shl 1
        }
        fail(table.invalidMessage)
    }

    private fun extra(count: Int): Int = if (count == 0) 0 else take(count)

    private fun take(count: Int): Int {
        ensure(count)
        val value = (bits and ((1L shl count) - 1)).toInt()
        bits = bits ushr count
        bitCount -= count
        return value
    }

    private fun nextBit(): Int {
        if (bitCount == 0) {
            if (inputPos >= inputEnd) fail("truncated input")
            bits = input[inputPos].toLong() and 0xFF
            inputPos++
            bitCount = 8
        }
        val bit = (bits and 1L).toInt()
        bits = bits ushr 1
        bitCount--
        return bit
    }

    private fun ensure(count: Int) {
        while (bitCount < count) {
            if (inputPos >= inputEnd) fail("truncated input")
            bits = bits or ((input[inputPos].toLong() and 0xFF) shl bitCount)
            inputPos++
            bitCount += 8
        }
    }

    private fun alignToByte() {
        val partial = bitCount and 7
        if (partial != 0) {
            bits = bits ushr partial
            bitCount -= partial
        }
        inputPos -= bitCount / 8
        bits = 0
        bitCount = 0
    }

    private fun readU16(at: Int): Int = (input[at].toInt() and 0xFF) or ((input[at + 1].toInt() and 0xFF) shl 8)

    private fun fail(message: String): Nothing = throw IllegalArgumentException(message)

    private companion object {
        val EMPTY = ByteArray(0)

        const val END_OF_BLOCK = 256
        const val MAX_MATCH = 258
        const val FAST_INPUT_BYTES = 8
        const val MAX_LITERAL_CODES = 288
        const val MAX_DISTANCE_CODES = 32
        const val PRECODE_SYMBOLS = 19

        val CODE_LENGTH_ORDER = intArrayOf(
            16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
        )

        val LENGTH_BASE = intArrayOf(
            3, 4, 5, 6, 7, 8, 9, 10,
            11, 13, 15, 17, 19, 23, 27, 31,
            35, 43, 51, 59, 67, 83, 99, 115,
            131, 163, 195, 227, 258,
        )
        val LENGTH_EXTRA = intArrayOf(
            0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 2, 2, 2, 2,
            3, 3, 3, 3, 4, 4, 4, 4,
            5, 5, 5, 5, 0,
        )
        val DISTANCE_BASE = intArrayOf(
            1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
            257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
        )
        val DISTANCE_EXTRA = intArrayOf(
            0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
            7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13,
        )
    }
}
