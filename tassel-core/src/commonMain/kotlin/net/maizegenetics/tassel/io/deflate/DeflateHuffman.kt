package net.maizegenetics.tassel.io.deflate

/**
 * RFC 1951 length and distance alphabets, and length-limited Huffman codes.
 *
 * Code lengths are at most 15. A single used symbol gets length 1, which is the one incomplete set RFC 1951
 * allows. Every other set is a complete prefix code. Ties in the heap break by symbol index.
 */
internal object DeflateCodes {
    const val LITERAL_SYMBOLS = 286
    const val DISTANCE_SYMBOLS = 30
    const val CODE_LENGTH_SYMBOLS = 19
    const val END_OF_BLOCK = 256
    const val MAX_CODE_LENGTH = 15
    const val MAX_MATCH = 258
    const val MAX_DISTANCE = 32768

    val lengthBase = intArrayOf(
        3, 4, 5, 6, 7, 8, 9, 10,
        11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115,
        131, 163, 195, 227, 258,
    )
    val lengthExtra = intArrayOf(
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4,
        5, 5, 5, 5, 0,
    )
    val distanceBase = intArrayOf(
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
    )
    val distanceExtra = intArrayOf(
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13,
    )

    /** Permutation of the 19 code-length symbols. RFC 1951 section 3.2.7. */
    val codeLengthOrder = intArrayOf(
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
    )

    val lengthCode = IntArray(MAX_MATCH + 1)
    val distanceCode = IntArray(MAX_DISTANCE + 1)

    init {
        for (index in lengthBase.indices) {
            val span = 1 shl lengthExtra[index]
            for (length in lengthBase[index] until lengthBase[index] + span) {
                lengthCode[length] = index
            }
        }
        for (index in distanceBase.indices) {
            val span = 1 shl distanceExtra[index]
            val last = distanceBase[index] + span
            for (distance in distanceBase[index] until last) {
                if (distance <= MAX_DISTANCE) distanceCode[distance] = index
            }
        }
    }

    fun fixedCode(symbol: Int): Int = when (symbol) {
        in 0..143 -> symbol + 48
        in 144..255 -> symbol - 144 + 400
        in 256..279 -> symbol - 256
        in 280..287 -> symbol - 280 + 192
        else -> error("symbol $symbol has no fixed literal code")
    }

    fun fixedLength(symbol: Int): Int = when (symbol) {
        in 0..143 -> 8
        in 144..255 -> 9
        in 256..279 -> 7
        in 280..287 -> 8
        else -> error("symbol $symbol has no fixed literal code")
    }
}

/**
 * Builds canonical Huffman code lengths and the codes those lengths assign.
 *
 * [lengths] overwrites `into[0, symbolCount)`. Symbols with a zero frequency get length 0. When a length would
 * exceed [maxLength], counts are redistributed until every leaf fits: one code at a shorter length is moved
 * down, which frees a prefix for two codes at the next length, and one of those replaces a code that was past
 * the limit. The longest codes then go to the least frequent symbols.
 */
internal class CodeLengthBuilder {
    private val freq = IntArray(MAX_NODES)
    private val minSymbol = IntArray(MAX_NODES)
    private val parent = IntArray(MAX_NODES)
    private val heap = IntArray(DeflateCodes.LITERAL_SYMBOLS)
    private val leaves = IntArray(DeflateCodes.LITERAL_SYMBOLS)
    private val counts = IntArray(MAX_TRACKED_LENGTH + 1)
    private val nextCode = IntArray(DeflateCodes.MAX_CODE_LENGTH + 1)

    fun lengths(frequencies: IntArray, symbolCount: Int, maxLength: Int, into: IntArray) {
        require(maxLength in 1..MAX_TRACKED_LENGTH) { "code length limit $maxLength is out of range" }
        into.fill(0, 0, symbolCount)
        var leafCount = 0
        for (symbol in 0 until symbolCount) {
            if (frequencies[symbol] != 0) leaves[leafCount++] = symbol
        }
        if (leafCount == 0) return
        if (leafCount == 1) {
            into[leaves[0]] = 1
            return
        }

        var nextNode = symbolCount
        var heapSize = 0
        for (index in 0 until leafCount) {
            val symbol = leaves[index]
            freq[symbol] = frequencies[symbol]
            minSymbol[symbol] = symbol
            parent[symbol] = -1
            heap[heapSize] = symbol
            siftUp(heapSize)
            heapSize++
        }
        while (heapSize > 1) {
            val first = heap[0]
            heapSize = removeRoot(heapSize)
            val second = heap[0]
            heapSize = removeRoot(heapSize)
            val node = nextNode++
            freq[node] = freq[first] + freq[second]
            minSymbol[node] = minOf(minSymbol[first], minSymbol[second])
            parent[first] = node
            parent[second] = node
            parent[node] = -1
            heap[heapSize] = node
            siftUp(heapSize)
            heapSize++
        }

        counts.fill(0, 0, maxLength + 1)
        var overflow = 0
        for (index in 0 until leafCount) {
            val symbol = leaves[index]
            var depth = 0
            var node = symbol
            while (parent[node] >= 0) {
                node = parent[node]
                depth++
            }
            if (depth > maxLength) {
                overflow++
                counts[maxLength]++
            } else {
                counts[depth]++
                into[symbol] = depth
            }
        }
        if (overflow == 0) {
            checkComplete(into, symbolCount, maxLength)
            return
        }

        do {
            var bits = maxLength - 1
            while (bits > 0 && counts[bits] == 0) bits--
            check(bits > 0) { "code lengths cannot be limited to $maxLength" }
            counts[bits]--
            counts[bits + 1] += 2
            counts[maxLength]--
            overflow -= 2
        } while (overflow > 0)

        sortLeaves(frequencies, leafCount)
        var cursor = 0
        for (bits in maxLength downTo 1) {
            var remaining = counts[bits]
            while (remaining > 0) {
                into[leaves[cursor]] = bits
                cursor++
                remaining--
            }
        }
        check(cursor == leafCount) { "code length redistribution lost a symbol" }
        checkComplete(into, symbolCount, maxLength)
    }

    /** Canonical codes for [lengths], most-significant bit first, as RFC 1951 section 3.2.2 numbers them. */
    fun codes(lengths: IntArray, symbolCount: Int, into: IntArray) {
        counts.fill(0, 0, DeflateCodes.MAX_CODE_LENGTH + 1)
        for (symbol in 0 until symbolCount) {
            val length = lengths[symbol]
            if (length != 0) counts[length]++
        }
        var code = 0
        for (length in 1..DeflateCodes.MAX_CODE_LENGTH) {
            code = (code + counts[length - 1]) shl 1
            nextCode[length] = code
        }
        for (symbol in 0 until symbolCount) {
            val length = lengths[symbol]
            if (length != 0) {
                into[symbol] = nextCode[length]
                nextCode[length]++
            }
        }
    }

    private fun sortLeaves(frequencies: IntArray, leafCount: Int) {
        for (index in 1 until leafCount) {
            val symbol = leaves[index]
            var cursor = index
            while (cursor > 0 && comesAfter(frequencies, leaves[cursor - 1], symbol)) {
                leaves[cursor] = leaves[cursor - 1]
                cursor--
            }
            leaves[cursor] = symbol
        }
    }

    private fun comesAfter(frequencies: IntArray, left: Int, right: Int): Boolean {
        val leftFreq = frequencies[left]
        val rightFreq = frequencies[right]
        if (leftFreq != rightFreq) return leftFreq > rightFreq
        return left > right
    }

    private fun removeRoot(heapSize: Int): Int {
        val last = heapSize - 1
        heap[0] = heap[last]
        if (last > 0) siftDown(last)
        return last
    }

    private fun siftUp(start: Int) {
        var index = start
        val node = heap[index]
        while (index > 0) {
            val parentIndex = (index - 1) ushr 1
            if (!less(node, heap[parentIndex])) break
            heap[index] = heap[parentIndex]
            index = parentIndex
        }
        heap[index] = node
    }

    private fun siftDown(size: Int) {
        var index = 0
        val node = heap[0]
        while (true) {
            val left = index * 2 + 1
            if (left >= size) break
            val right = left + 1
            val child = if (right < size && less(heap[right], heap[left])) right else left
            if (!less(heap[child], node)) break
            heap[index] = heap[child]
            index = child
        }
        heap[index] = node
    }

    private fun less(left: Int, right: Int): Boolean {
        val leftFreq = freq[left]
        val rightFreq = freq[right]
        if (leftFreq != rightFreq) return leftFreq < rightFreq
        return minSymbol[left] < minSymbol[right]
    }

    private fun checkComplete(lengths: IntArray, symbolCount: Int, maxLength: Int) {
        if (maxLength > DeflateCodes.MAX_CODE_LENGTH) return
        counts.fill(0, 0, DeflateCodes.MAX_CODE_LENGTH + 1)
        var used = 0
        for (symbol in 0 until symbolCount) {
            val length = lengths[symbol]
            if (length != 0) {
                counts[length]++
                used++
            }
        }
        if (used <= 1) return
        var left = 1
        for (length in 1..DeflateCodes.MAX_CODE_LENGTH) {
            left = (left shl 1) - counts[length]
            check(left >= 0) { "over-subscribed code lengths" }
        }
        check(left == 0) { "incomplete code lengths" }
    }

    private companion object {
        const val MAX_TRACKED_LENGTH = 64
        const val MAX_NODES = DeflateCodes.LITERAL_SYMBOLS * 2
    }
}
