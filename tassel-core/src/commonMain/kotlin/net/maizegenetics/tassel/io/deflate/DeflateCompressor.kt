package net.maizegenetics.tassel.io.deflate

import net.maizegenetics.tassel.io.BgzfCompressor

/**
 * One fixed DEFLATE level for a single BGZF slice (RFC 1951).
 *
 * The search is zlib's level-6 lazy table: 15-bit hash chains, chain cap 128, and a chain stops once a match
 * of 128 bytes is found. A match shorter than 16 bytes is kept only when the next position is not a longer one.
 * When a held match is already 8 bytes, the next chain is a quarter of the cap. Matches are 3 to 258 bytes and
 * look back at most 32,768 bytes. Each accepted match inserts every position inside it. One call emits one final
 * block, fixed or dynamic Huffman, whichever uses fewer bits. A payload that would not fit in a BGZF member
 * returns null so the writer stores the slice.
 *
 * Scratch buffers live on the instance. [compress] is not safe to call concurrently on the same instance.
 */
internal class DeflateCompressor : BgzfCompressor {
    private val head = IntArray(HASH_SIZE)
    private val prev = IntArray(BGZF_BLOCK_SIZE)
    private val tokenLen = IntArray(BGZF_BLOCK_SIZE)
    private val tokenDist = IntArray(BGZF_BLOCK_SIZE)
    private val litFreq = IntArray(DeflateCodes.LITERAL_SYMBOLS)
    private val distFreq = IntArray(DeflateCodes.DISTANCE_SYMBOLS)
    private val litLengths = IntArray(DeflateCodes.LITERAL_SYMBOLS)
    private val distLengths = IntArray(DeflateCodes.DISTANCE_SYMBOLS)
    private val litCodes = IntArray(DeflateCodes.LITERAL_SYMBOLS)
    private val distCodes = IntArray(DeflateCodes.DISTANCE_SYMBOLS)
    private val lengthSequence = IntArray(SEQUENCE_SYMBOLS)
    private val rleSymbol = IntArray(SEQUENCE_SYMBOLS)
    private val rleExtra = IntArray(SEQUENCE_SYMBOLS)
    private val clFreq = IntArray(DeflateCodes.CODE_LENGTH_SYMBOLS)
    private val clLengths = IntArray(DeflateCodes.CODE_LENGTH_SYMBOLS)
    private val clCodes = IntArray(DeflateCodes.CODE_LENGTH_SYMBOLS)
    private val payload = ByteArray(PAYLOAD_LIMIT)
    private val bits = DeflateBitSink()
    private val huffman = CodeLengthBuilder()

    private var tokenCount = 0
    private var litCount = 0
    private var distCount = 0
    private var codeLengthCount = 0
    private var rleCount = 0

    override fun compress(input: ByteArray, from: Int, to: Int): ByteArray? {
        require(from >= 0 && to >= from && to <= input.size) { "invalid deflate range" }
        val count = to - from
        require(count <= BGZF_BLOCK_SIZE) { "deflate block is longer than $BGZF_BLOCK_SIZE bytes" }
        tokenize(input, from, count)
        countFrequencies()
        prepareDynamic()
        val fixedBits = encode(null, fixed = true)
        val dynamicBits = encode(null, fixed = false)
        val useFixed = fixedBits <= dynamicBits
        val bitCount = if (useFixed) fixedBits else dynamicBits
        val bytes = (bitCount + 7) ushr 3
        if (bytes > PAYLOAD_LIMIT) return null
        bits.reset(payload)
        val writtenBits = encode(bits, fixed = useFixed)
        val written = bits.finish()
        check(writtenBits == bitCount && written == bytes) {
            "deflate wrote $written bytes ($writtenBits bits), counted $bytes bytes ($bitCount bits)"
        }
        return payload.copyOf(written)
    }

    private fun tokenize(input: ByteArray, from: Int, count: Int) {
        tokenCount = 0
        if (count < MIN_MATCH) {
            for (offset in 0 until count) emitLiteral(input[from + offset].toInt() and 0xFF)
            return
        }
        head.fill(0)
        var hash = input[from].toInt() and 0xFF
        hash = ((hash shl HASH_SHIFT) xor (input[from + 1].toInt() and 0xFF)) and HASH_MASK
        var pos = 0
        var matchLen = 0
        var matchDist = 0
        var pending = false
        while (pos < count) {
            val available = count - pos
            if (available < MIN_MATCH) break
            hash = ((hash shl HASH_SHIFT) xor (input[from + pos + 2].toInt() and 0xFF)) and HASH_MASK
            val chainHead = head[hash]
            head[hash] = pos + 1
            prev[pos] = chainHead
            val prevLen = matchLen
            val prevDist = matchDist
            if (prevLen >= MAX_LAZY) {
                matchLen = 0
                matchDist = 0
            } else {
                val chainLimit = if (prevLen >= GOOD_LENGTH) MAX_CHAIN ushr 2 else MAX_CHAIN
                val maxLen = minOf(MAX_MATCH, available)
                var bestLen = 0
                var bestDist = 0
                var chain = chainHead
                var steps = 0
                while (chain != 0 && steps < chainLimit && bestLen < maxLen) {
                    val matchPos = chain - 1
                    val distance = pos - matchPos
                    if (distance <= 0 || distance > WINDOW) break
                    if (input[from + matchPos + bestLen] == input[from + pos + bestLen]) {
                        var length = 0
                        while (length < maxLen && input[from + pos + length] == input[from + matchPos + length]) {
                            length++
                        }
                        if (length > bestLen) {
                            bestLen = length
                            bestDist = distance
                            if (bestLen >= NICE_LENGTH) break
                        }
                    }
                    chain = prev[matchPos]
                    steps++
                }
                matchLen = bestLen
                matchDist = bestDist
            }
            if (prevLen >= MIN_MATCH && matchLen <= prevLen) {
                emitMatch(prevLen, prevDist)
                val matchEnd = pos - 1 + prevLen
                var at = pos + 1
                while (at < matchEnd && count - at >= MIN_MATCH) {
                    val next = input[from + at + 2].toInt() and 0xFF
                    hash = ((hash shl HASH_SHIFT) xor next) and HASH_MASK
                    val previous = head[hash]
                    head[hash] = at + 1
                    prev[at] = previous
                    at++
                }
                pos = matchEnd
                matchLen = 0
                matchDist = 0
                pending = false
            } else if (pending) {
                emitLiteral(input[from + pos - 1].toInt() and 0xFF)
                pos++
            } else {
                pending = true
                pos++
            }
        }
        if (pending && matchLen >= MIN_MATCH) {
            emitMatch(matchLen, matchDist)
            pos = pos - 1 + matchLen
        } else if (pending && pos > 0) {
            emitLiteral(input[from + pos - 1].toInt() and 0xFF)
        }
        while (pos < count) {
            emitLiteral(input[from + pos].toInt() and 0xFF)
            pos++
        }
    }

    private fun emitLiteral(symbol: Int) {
        tokenLen[tokenCount] = symbol
        tokenDist[tokenCount] = 0
        tokenCount++
    }

    private fun emitMatch(length: Int, distance: Int) {
        tokenLen[tokenCount] = length
        tokenDist[tokenCount] = distance
        tokenCount++
    }

    private fun countFrequencies() {
        litFreq.fill(0)
        distFreq.fill(0)
        litFreq[DeflateCodes.END_OF_BLOCK] = 1
        for (index in 0 until tokenCount) {
            if (tokenDist[index] == 0) {
                litFreq[tokenLen[index]]++
            } else {
                val lengthIndex = DeflateCodes.lengthCode[tokenLen[index]]
                val distanceIndex = DeflateCodes.distanceCode[tokenDist[index]]
                litFreq[257 + lengthIndex]++
                distFreq[distanceIndex]++
            }
        }
        var usedDistance = false
        for (freq in distFreq) {
            if (freq != 0) usedDistance = true
        }
        if (!usedDistance) distFreq[0] = 1
    }

    private fun prepareDynamic() {
        huffman.lengths(litFreq, DeflateCodes.LITERAL_SYMBOLS, DeflateCodes.MAX_CODE_LENGTH, litLengths)
        huffman.lengths(distFreq, DeflateCodes.DISTANCE_SYMBOLS, DeflateCodes.MAX_CODE_LENGTH, distLengths)
        huffman.codes(litLengths, DeflateCodes.LITERAL_SYMBOLS, litCodes)
        huffman.codes(distLengths, DeflateCodes.DISTANCE_SYMBOLS, distCodes)
        check(litLengths[DeflateCodes.END_OF_BLOCK] != 0) { "end of block has no code" }

        var literals = DeflateCodes.LITERAL_SYMBOLS
        while (literals > 257 && litLengths[literals - 1] == 0) literals--
        var distances = DeflateCodes.DISTANCE_SYMBOLS
        while (distances > 1 && distLengths[distances - 1] == 0) distances--
        litCount = literals
        distCount = distances

        for (index in 0 until literals) lengthSequence[index] = litLengths[index]
        for (index in 0 until distances) lengthSequence[literals + index] = distLengths[index]
        rleCount = runLengthEncode(literals + distances)

        clFreq.fill(0)
        for (index in 0 until rleCount) clFreq[rleSymbol[index]]++
        // The code-length alphabet is written as 3-bit lengths, so those codes are at most 7 bits.
        huffman.lengths(clFreq, DeflateCodes.CODE_LENGTH_SYMBOLS, MAX_CODE_LENGTH_BITS, clLengths)
        huffman.codes(clLengths, DeflateCodes.CODE_LENGTH_SYMBOLS, clCodes)

        var codes = DeflateCodes.CODE_LENGTH_SYMBOLS
        val order = DeflateCodes.codeLengthOrder
        while (codes > 4 && clLengths[order[codes - 1]] == 0) codes--
        codeLengthCount = codes
    }

    private fun runLengthEncode(count: Int): Int {
        var commands = 0
        var index = 0
        while (index < count) {
            val length = lengthSequence[index]
            var run = 1
            while (index + run < count && lengthSequence[index + run] == length) run++
            if (length == 0) {
                var left = run
                while (left >= 11) {
                    val take = minOf(left, 138)
                    rleSymbol[commands] = 18
                    rleExtra[commands] = take - 11
                    commands++
                    left -= take
                }
                if (left >= 3) {
                    rleSymbol[commands] = 17
                    rleExtra[commands] = left - 3
                    commands++
                } else {
                    repeat(left) {
                        rleSymbol[commands] = 0
                        rleExtra[commands] = 0
                        commands++
                    }
                }
            } else {
                rleSymbol[commands] = length
                rleExtra[commands] = 0
                commands++
                var left = run - 1
                while (left >= 3) {
                    val take = minOf(left, 6)
                    rleSymbol[commands] = 16
                    rleExtra[commands] = take - 3
                    commands++
                    left -= take
                }
                while (left > 0) {
                    rleSymbol[commands] = length
                    rleExtra[commands] = 0
                    commands++
                    left--
                }
            }
            index += run
        }
        return commands
    }

    /**
     * Bit length of one final block. When [sink] is null, only the length is counted, so the fixed and dynamic
     * choices can be compared before either is written.
     */
    private fun encode(sink: DeflateBitSink?, fixed: Boolean): Int {
        var total = 0
        fun lsb(value: Int, length: Int) {
            total += length
            sink?.writeLsb(value, length)
        }
        fun msb(code: Int, length: Int) {
            check(length > 0) { "missing Huffman code" }
            total += length
            sink?.writeMsb(code, length)
        }
        lsb(1, 1)
        if (fixed) {
            lsb(1, 2)
        } else {
            lsb(2, 2)
            lsb(litCount - 257, 5)
            lsb(distCount - 1, 5)
            lsb(codeLengthCount - 4, 4)
            val order = DeflateCodes.codeLengthOrder
            for (index in 0 until codeLengthCount) lsb(clLengths[order[index]], 3)
            for (index in 0 until rleCount) {
                val symbol = rleSymbol[index]
                msb(clCodes[symbol], clLengths[symbol])
                val extra = repeatExtra(symbol)
                if (extra != 0) lsb(rleExtra[index], extra)
            }
        }
        for (index in 0 until tokenCount) {
            val distance = tokenDist[index]
            if (distance == 0) {
                val symbol = tokenLen[index]
                if (fixed) {
                    msb(DeflateCodes.fixedCode(symbol), DeflateCodes.fixedLength(symbol))
                } else {
                    msb(litCodes[symbol], litLengths[symbol])
                }
            } else {
                writeMatch(tokenLen[index], distance, fixed, ::lsb, ::msb)
            }
        }
        if (fixed) {
            msb(DeflateCodes.fixedCode(DeflateCodes.END_OF_BLOCK), 7)
        } else {
            msb(litCodes[DeflateCodes.END_OF_BLOCK], litLengths[DeflateCodes.END_OF_BLOCK])
        }
        return total
    }

    private fun writeMatch(
        length: Int,
        distance: Int,
        fixed: Boolean,
        lsb: (Int, Int) -> Unit,
        msb: (Int, Int) -> Unit,
    ) {
        val lengthIndex = DeflateCodes.lengthCode[length]
        val distanceIndex = DeflateCodes.distanceCode[distance]
        val lengthSymbol = 257 + lengthIndex
        if (fixed) {
            msb(DeflateCodes.fixedCode(lengthSymbol), DeflateCodes.fixedLength(lengthSymbol))
        } else {
            msb(litCodes[lengthSymbol], litLengths[lengthSymbol])
        }
        val lengthExtra = DeflateCodes.lengthExtra[lengthIndex]
        if (lengthExtra != 0) lsb(length - DeflateCodes.lengthBase[lengthIndex], lengthExtra)
        if (fixed) {
            msb(distanceIndex, 5)
        } else {
            msb(distCodes[distanceIndex], distLengths[distanceIndex])
        }
        val distanceExtra = DeflateCodes.distanceExtra[distanceIndex]
        if (distanceExtra != 0) lsb(distance - DeflateCodes.distanceBase[distanceIndex], distanceExtra)
    }

    private fun repeatExtra(symbol: Int): Int = when (symbol) {
        16 -> 2
        17 -> 3
        18 -> 7
        else -> 0
    }

    private companion object {
        const val MIN_MATCH = 3
        const val MAX_MATCH = DeflateCodes.MAX_MATCH
        const val WINDOW = DeflateCodes.MAX_DISTANCE
        const val HASH_BITS = 15
        const val HASH_SHIFT = 5
        const val HASH_SIZE = 1 shl HASH_BITS
        const val HASH_MASK = HASH_SIZE - 1
        const val MAX_CHAIN = 128
        const val NICE_LENGTH = 128
        const val MAX_LAZY = 16
        const val GOOD_LENGTH = 8
        const val MAX_CODE_LENGTH_BITS = 7
        const val SEQUENCE_SYMBOLS = DeflateCodes.LITERAL_SYMBOLS + DeflateCodes.DISTANCE_SYMBOLS
        const val PAYLOAD_LIMIT = BGZF_MAX_BLOCK - BGZF_HEADER_LENGTH - GZIP_TRAILER
    }
}
