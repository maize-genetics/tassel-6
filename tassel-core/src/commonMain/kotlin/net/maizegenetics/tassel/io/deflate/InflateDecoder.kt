package net.maizegenetics.tassel.io.deflate

/**
 * Raw DEFLATE decoder (RFC 1951).
 *
 * [inflate] decodes a complete stream from `input[from, to)` into `output[at, at + capacity)` and returns
 * the number of bytes written. That is the BGZF path: the whole member is in hand, and the output buffer is
 * the history. It does not allocate.
 *
 * [inflateSome] is the plain-gzip path. It stops when the input or the output runs out and continues on the
 * next call. History lives in a 32 KiB window, allocated on the first [resetStream] or [inflateSome], so a
 * member of any size stays bounded by that window and the caller's buffers. Bytes reported as consumed are
 * not presented again; leftover bits stay in this decoder. [InflateStep.returned] is only set when a stream
 * finishes and a previous call had already consumed bytes that belong to the caller (the gzip trailer).
 *
 * The same instance is reused across blocks: the fixed tables are built once, and each dynamic block
 * overwrites the scratch tables. Not safe for concurrent calls on one instance. A failed call can leave the
 * scratch tables half-built; the next [inflate] or [resetStream] rebuilds whatever it needs.
 *
 * Gzip and BGZF framing stay outside this type.
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

    private var streamFrom: Int = 0
    private var endOfInput: Boolean = false
    private var window: ByteArray = EMPTY
    private var windowPos: Int = 0
    private var windowFilled: Int = 0
    private var phase: Phase = Phase.BLOCK
    private var blockStage: Int = 0
    private var blockLast: Boolean = false
    private var streamFixed: Boolean = false
    private var storedStage: Int = 0
    private var storedLeft: Int = 0
    private var dynStage: Int = 0
    private var dynIndex: Int = 0
    private var literalCount: Int = 0
    private var distanceCount: Int = 0
    private var precodeCount: Int = 0
    private var codeTotal: Int = 0
    private var pendingRepeat: Int = -1
    private var matchStage: Int = 0
    private var lengthSymbol: Int = 0
    private var distanceSymbol: Int = 0
    private var matchLength: Int = 0
    private var matchDistance: Int = 0
    private var matchLeft: Int = 0

    /** Bytes a finished stream had already taken that belong back to the caller, oldest first. */
    internal val returnedBytes = ByteArray(RETURNED_BYTES)
    internal var returnedCount: Int = 0
        private set

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

    /**
     * Starts a new raw DEFLATE stream on this instance. Call it before each gzip member.
     * The 32 KiB window is allocated here if this decoder has not streamed before.
     */
    fun resetStream() {
        if (window.size != WINDOW) window = ByteArray(WINDOW)
        windowPos = 0
        windowFilled = 0
        bits = 0
        bitCount = 0
        phase = Phase.BLOCK
        blockStage = 0
        blockLast = false
        streamFixed = false
        storedStage = 0
        storedLeft = 0
        dynStage = 0
        dynIndex = 0
        pendingRepeat = -1
        matchStage = 0
        matchLeft = 0
        returnedCount = 0
    }

    /**
     * Decodes the next slice of the stream started by [resetStream].
     *
     * [endOfInput] is true only when [input] from [from] to [to] is the last input that will ever arrive.
     * A short buffer with [endOfInput] false returns; the same condition with [endOfInput] true is a truncated
     * stream. A full output buffer returns with the match or stored copy still open.
     *
     * @throws IllegalArgumentException on a bad range or a stream that is not valid raw DEFLATE.
     */
    fun inflateSome(
        input: ByteArray,
        from: Int,
        to: Int,
        output: ByteArray,
        at: Int,
        capacity: Int,
        endOfInput: Boolean,
    ): InflateStep {
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
        if (window.size != WINDOW) window = ByteArray(WINDOW)
        this.input = input
        inputPos = from
        inputEnd = to
        this.output = output
        outputPos = at
        outputEnd = at + capacity
        streamFrom = from
        this.endOfInput = endOfInput
        returnedCount = 0
        try {
            if (phase != Phase.DONE) runStream()
            return InflateStep(inputPos - from, outputPos - at, phase == Phase.DONE, returnedCount)
        } finally {
            this.input = EMPTY
            this.output = EMPTY
        }
    }

    private fun runStream() {
        while (phase != Phase.DONE) {
            val more = when (phase) {
                Phase.BLOCK -> decodeBlock()
                Phase.STORED_LEN -> decodeStoredHeader()
                Phase.STORED_COPY -> copyStored()
                Phase.DYNAMIC -> decodeDynamicStream()
                Phase.HUFFMAN -> decodeHuffmanStream()
                Phase.MATCH -> decodeMatch()
                Phase.DONE -> false
            }
            if (!more) return
        }
    }

    private fun decodeBlock(): Boolean {
        when (blockStage) {
            0 -> {
                blockLast = (pull(1) ?: return false) == 1
                blockStage = 1
                return true
            }
            else -> {
                when (val type = pull(2) ?: return false) {
                    0 -> {
                        byteAlignStream()
                        storedStage = 0
                        phase = Phase.STORED_LEN
                    }
                    1 -> {
                        streamFixed = true
                        phase = Phase.HUFFMAN
                    }
                    2 -> {
                        streamFixed = false
                        dynStage = 0
                        pendingRepeat = -1
                        phase = Phase.DYNAMIC
                    }
                    else -> fail("invalid block type")
                }
                blockStage = 0
                return true
            }
        }
    }

    private fun decodeStoredHeader(): Boolean {
        when (storedStage) {
            0 -> {
                storedLeft = pull(16) ?: return false
                storedStage = 1
                return true
            }
            else -> {
                val complement = pull(16) ?: return false
                if ((storedLeft xor complement) != 0xFFFF) {
                    fail("stored block length does not match its complement")
                }
                storedStage = 0
                if (storedLeft == 0) completeBlock() else phase = Phase.STORED_COPY
                return true
            }
        }
    }

    private fun copyStored(): Boolean {
        val history = window
        while (storedLeft > 0) {
            if (outputPos == outputEnd) return false
            if (bitCount >= 8) {
                val value = (bits and 0xFF).toInt()
                bits = bits ushr 8
                bitCount -= 8
                emit(value)
                storedLeft--
                continue
            }
            if (bitCount != 0) fail("truncated input")
            if (inputPos >= inputEnd) {
                if (endOfInput) fail("truncated input")
                return false
            }
            val chunk = minOf(storedLeft, outputEnd - outputPos, inputEnd - inputPos)
            var index = 0
            while (index < chunk) {
                val value = input[inputPos + index].toInt() and 0xFF
                history[windowPos] = value.toByte()
                windowPos = (windowPos + 1) and WINDOW_MASK
                if (windowFilled < WINDOW) windowFilled++
                output[outputPos + index] = value.toByte()
                index++
            }
            inputPos += chunk
            outputPos += chunk
            storedLeft -= chunk
        }
        completeBlock()
        return true
    }

    private fun decodeDynamicStream(): Boolean {
        when (dynStage) {
            0 -> {
                literalCount = (pull(5) ?: return false) + 257
                dynStage = 1
                return true
            }
            1 -> {
                distanceCount = (pull(5) ?: return false) + 1
                dynStage = 2
                return true
            }
            2 -> {
                precodeCount = (pull(4) ?: return false) + 4
                precodeLengths.fill(0)
                dynIndex = 0
                dynStage = 3
                return true
            }
            3 -> {
                while (dynIndex < precodeCount) {
                    precodeLengths[CODE_LENGTH_ORDER[dynIndex]] = pull(3) ?: return false
                    dynIndex++
                }
                precode.build(precodeLengths, 0, PRECODE_SYMBOLS, allowEmpty = false)
                codeTotal = literalCount + distanceCount
                dynIndex = 0
                pendingRepeat = -1
                dynStage = 4
                return true
            }
            else -> {
                while (dynIndex < codeTotal) {
                    if (pendingRepeat >= 0) {
                        val extra = pull(repeatBits(pendingRepeat)) ?: return false
                        applyRepeat(pendingRepeat, extra)
                        pendingRepeat = -1
                        continue
                    }
                    when (val symbol = symbolStream(precode) ?: return false) {
                        in 0..15 -> codeLengths[dynIndex++] = symbol
                        16, 17, 18 -> {
                            if (symbol == 16 && dynIndex == 0) fail("invalid length symbol")
                            val extra = pull(repeatBits(symbol))
                            if (extra == null) {
                                pendingRepeat = symbol
                                return false
                            }
                            applyRepeat(symbol, extra)
                        }
                        else -> fail("invalid length symbol")
                    }
                }
                literal.build(codeLengths, 0, literalCount, allowEmpty = false)
                distance.build(codeLengths, literalCount, distanceCount, allowEmpty = true)
                streamFixed = false
                phase = Phase.HUFFMAN
                dynStage = 0
                return true
            }
        }
    }

    private fun applyRepeat(symbol: Int, extra: Int) {
        val value = if (symbol == 16) codeLengths[dynIndex - 1] else 0
        val count = extra + if (symbol == 18) 11 else 3
        if (dynIndex + count > codeTotal) fail("invalid length symbol")
        repeat(count) { codeLengths[dynIndex++] = value }
    }

    private fun repeatBits(symbol: Int): Int = when (symbol) {
        16 -> 2
        17 -> 3
        else -> 7
    }

    private fun decodeHuffmanStream(): Boolean {
        val literals = if (streamFixed) fixedLiteral else literal
        while (phase == Phase.HUFFMAN) {
            if (outputPos == outputEnd) return false
            when (val symbol = symbolStream(literals) ?: return false) {
                in 0..255 -> emit(symbol)
                END_OF_BLOCK -> {
                    completeBlock()
                    return true
                }
                in 257..285 -> {
                    lengthSymbol = symbol
                    matchLength = LENGTH_BASE[symbol - 257]
                    matchStage = 0
                    phase = Phase.MATCH
                    return true
                }
                else -> fail("invalid length symbol")
            }
        }
        return true
    }

    private fun decodeMatch(): Boolean {
        val distances = if (streamFixed) fixedDistance else distance
        when (matchStage) {
            0 -> {
                val extraCount = LENGTH_EXTRA[lengthSymbol - 257]
                if (extraCount != 0) matchLength += pull(extraCount) ?: return false
                matchStage = 1
                return true
            }
            1 -> {
                if (distances.empty) fail("invalid distance symbol")
                val symbol = symbolStream(distances) ?: return false
                if (symbol > 29) fail("invalid distance symbol")
                distanceSymbol = symbol
                matchDistance = DISTANCE_BASE[symbol]
                matchStage = 2
                return true
            }
            2 -> {
                val extraCount = DISTANCE_EXTRA[distanceSymbol]
                if (extraCount != 0) matchDistance += pull(extraCount) ?: return false
                if (matchDistance > windowFilled) fail("distance past the start of the output")
                matchLeft = matchLength
                matchStage = 3
                return true
            }
            else -> {
                val history = window
                while (matchLeft > 0) {
                    if (outputPos == outputEnd) return false
                    val src = (windowPos - matchDistance) and WINDOW_MASK
                    emit(history[src].toInt() and 0xFF)
                    matchLeft--
                }
                phase = Phase.HUFFMAN
                return true
            }
        }
    }

    private fun completeBlock() {
        if (blockLast) finishStream() else phase = Phase.BLOCK
    }

    /**
     * Drops the 0–7 padding bits in the last DEFLATE byte and gives back any whole bytes that were prefetched
     * past the end of the stream. Those bytes are the gzip trailer. Rewinding only reaches bytes pulled during
     * this call; older ones are copied to [returnedBytes] because the caller has already been told they were
     * consumed.
     */
    private fun finishStream() {
        val partial = bitCount and 7
        if (partial != 0) {
            bits = bits ushr partial
            bitCount -= partial
        }
        val extra = bitCount / 8
        val pulled = inputPos - streamFrom
        val rewind = minOf(extra, pulled)
        inputPos -= rewind
        val putback = extra - rewind
        if (putback > returnedBytes.size) fail("truncated input")
        var hold = bits
        var index = 0
        while (index < putback) {
            returnedBytes[index] = (hold and 0xFF).toByte()
            hold = hold ushr 8
            index++
        }
        returnedCount = putback
        bits = 0
        bitCount = 0
        phase = Phase.DONE
    }

    private fun emit(value: Int) {
        window[windowPos] = value.toByte()
        windowPos = (windowPos + 1) and WINDOW_MASK
        if (windowFilled < WINDOW) windowFilled++
        output[outputPos++] = value.toByte()
    }

    /** Drops 0–7 bits so the next pull is byte-aligned. Whole bytes already in the accumulator stay there. */
    private fun byteAlignStream() {
        val partial = bitCount and 7
        if (partial != 0) {
            bits = bits ushr partial
            bitCount -= partial
        }
    }

    private fun symbolStream(table: HuffmanTable): Int? {
        if (bitCount < 15) {
            while (bitCount < 15 && inputPos < inputEnd) {
                bits = bits or ((input[inputPos].toLong() and 0xFF) shl bitCount)
                inputPos++
                bitCount += 8
            }
        }
        if (bitCount >= 15) {
            val entry = table.lookup((bits and 0x7FFF).toInt())
            if (entry and HuffmanTable.FLAG_INVALID != 0) fail(table.invalidMessage)
            val length = (entry ushr HuffmanTable.LEN_SHIFT) and 0xF
            if (length == 0 || length > bitCount) fail("truncated input")
            bits = bits ushr length
            bitCount -= length
            return entry and HuffmanTable.PAYLOAD_MASK
        }
        if (bitCount == 0 && inputPos >= inputEnd) {
            if (endOfInput) fail("truncated input")
            return null
        }
        return symbolSlowStream(table)
    }

    /** Peeks bits and commits them only when a symbol is complete, so a short input can be resumed. */
    private fun symbolSlowStream(table: HuffmanTable): Int? {
        var code = 0
        var first = 0
        var index = 0
        var taken = 0
        for (length in 1..HuffmanTable.MAX_CODE_LENGTH) {
            if (taken == bitCount) {
                if (inputPos >= inputEnd) {
                    if (endOfInput) fail("truncated input")
                    return null
                }
                bits = bits or ((input[inputPos].toLong() and 0xFF) shl bitCount)
                inputPos++
                bitCount += 8
            }
            code = code or ((bits ushr taken) and 1L).toInt()
            taken++
            val count = table.counts[length]
            if (code - count < first) {
                val symbol = table.symbols[index + (code - first)]
                bits = bits ushr taken
                bitCount -= taken
                return symbol
            }
            index += count
            first += count
            first = first shl 1
            code = code shl 1
        }
        fail(table.invalidMessage)
    }

    /** Next [count] bits, or null when more input is required. [endOfInput] turns that shortage into an error. */
    private fun pull(count: Int): Int? {
        if (!ensureBits(count)) return null
        val value = (bits and ((1L shl count) - 1)).toInt()
        bits = bits ushr count
        bitCount -= count
        return value
    }

    private fun ensureBits(count: Int): Boolean {
        while (bitCount < count) {
            if (inputPos >= inputEnd) {
                if (endOfInput) fail("truncated input")
                return false
            }
            bits = bits or ((input[inputPos].toLong() and 0xFF) shl bitCount)
            inputPos++
            bitCount += 8
        }
        return true
    }

    private enum class Phase {
        BLOCK,
        STORED_LEN,
        STORED_COPY,
        DYNAMIC,
        HUFFMAN,
        MATCH,
        DONE,
    }

    private companion object {
        val EMPTY = ByteArray(0)

        const val WINDOW = 32768
        const val WINDOW_MASK = WINDOW - 1
        const val RETURNED_BYTES = 8
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

/**
 * One resume of [InflateDecoder.inflateSome].
 *
 * [consumed] bytes of this call's input are kept by the decoder. [written] bytes were appended to the output.
 * [finished] means the DEFLATE stream has ended and the next input byte is outside it. [returned] is the number of
 * bytes in [InflateDecoder.returnedBytes] that a previous call consumed and that the caller must read before the
 * unconsumed tail of this call.
 */
internal class InflateStep(val consumed: Int, val written: Int, val finished: Boolean, val returned: Int)
