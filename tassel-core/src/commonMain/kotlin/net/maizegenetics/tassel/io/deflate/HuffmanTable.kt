package net.maizegenetics.tassel.io.deflate

/**
 * Canonical Huffman decode table for one DEFLATE alphabet.
 *
 * The primary table is indexed by the next [tableBits] bits of the LSB-first stream. Codes longer than that
 * spill into one subtable, allocated from the same [entries] array. [build] overwrites the table in place.
 *
 * A length set is accepted when it is complete. It is also accepted when exactly one symbol has a non-zero
 * length: RFC 1951 allows that for distances, and real compressors emit it for a literal alphabet that uses
 * one symbol. Every other incomplete or over-subscribed set is an error. An empty set is accepted only when
 * [build] is called with `allowEmpty` (a distance alphabet with no matches).
 */
internal class HuffmanTable(capacity: Int, val tableBits: Int, val invalidMessage: String) {
    val entries = IntArray(capacity)
    val primaryMask: Int = (1 shl tableBits) - 1

    /** Code-length counts, kept for the bit-by-bit tail decoder. Index 0 stays zero. */
    val counts = IntArray(MAX_CODE_LENGTH + 1)

    /** Symbols in canonical order, grouped by increasing code length. */
    val symbols = IntArray(MAX_SYMBOLS)

    var empty: Boolean = true
        private set

    private val nextCode = IntArray(MAX_CODE_LENGTH + 1)
    private val prefixMax = IntArray(1 shl tableBits)
    private val reversed = IntArray(MAX_SYMBOLS)
    private val offsets = IntArray(MAX_CODE_LENGTH + 1)

    /**
     * Builds a table from `lengths[offset, offset + symbolCount)`.
     *
     * @throws IllegalArgumentException if the lengths are over-subscribed or incomplete, or if `symbolCount`
     *   exceeds [MAX_SYMBOLS].
     */
    fun build(lengths: IntArray, offset: Int, symbolCount: Int, allowEmpty: Boolean) {
        if (symbolCount > MAX_SYMBOLS) fail("invalid length symbol")
        entries.fill(INVALID)
        counts.fill(0)
        var used = 0
        var maxLength = 0
        for (symbol in 0 until symbolCount) {
            val length = lengths[offset + symbol]
            if (length !in 0..MAX_CODE_LENGTH) fail("invalid length symbol")
            if (length != 0) {
                counts[length]++
                used++
                if (length > maxLength) maxLength = length
            }
        }
        if (used == 0) {
            empty = true
            if (!allowEmpty) fail("incomplete code lengths")
            return
        }
        empty = false
        var left = 1
        for (length in 1..MAX_CODE_LENGTH) {
            left = left shl 1
            left -= counts[length]
            if (left < 0) fail("over-subscribed code lengths")
        }
        if (left > 0 && used != 1) fail("incomplete code lengths")

        var code = 0
        for (length in 1..MAX_CODE_LENGTH) {
            code = (code + counts[length - 1]) shl 1
            nextCode[length] = code
        }
        for (symbol in 0 until symbolCount) {
            val length = lengths[offset + symbol]
            if (length != 0) {
                reversed[symbol] = bitReverse(nextCode[length], length)
                nextCode[length]++
            }
        }

        val primarySize = 1 shl tableBits
        var nextFree = primarySize
        if (maxLength > tableBits) {
            prefixMax.fill(0)
            for (symbol in 0 until symbolCount) {
                val length = lengths[offset + symbol]
                if (length > tableBits) {
                    val prefix = reversed[symbol] and primaryMask
                    if (length > prefixMax[prefix]) prefixMax[prefix] = length
                }
            }
            for (prefix in 0 until primarySize) {
                val length = prefixMax[prefix]
                if (length > tableBits) {
                    val extra = length - tableBits
                    val size = 1 shl extra
                    if (nextFree + size > entries.size) fail("huffman table overflow")
                    entries[prefix] = nextFree or (extra shl EXTRA_SHIFT) or FLAG_SUBTABLE
                    entries.fill(INVALID, nextFree, nextFree + size)
                    nextFree += size
                }
            }
        }

        for (symbol in 0 until symbolCount) {
            val length = lengths[offset + symbol]
            if (length == 0) continue
            val rev = reversed[symbol]
            val packed = symbol or (length shl LEN_SHIFT)
            if (length <= tableBits) {
                val step = 1 shl length
                var index = rev
                while (index < primarySize) {
                    entries[index] = packed
                    index += step
                }
            } else {
                val pointer = entries[rev and primaryMask]
                val extra = (pointer ushr EXTRA_SHIFT) and 0xF
                val base = pointer and PAYLOAD_MASK
                val usedExtra = length - tableBits
                val step = 1 shl usedExtra
                var index = base + ((rev ushr tableBits) and ((1 shl extra) - 1))
                val end = base + (1 shl extra)
                while (index < end) {
                    entries[index] = packed
                    index += step
                }
            }
        }

        var sum = 0
        for (length in 1..MAX_CODE_LENGTH) {
            offsets[length] = sum
            sum += counts[length]
        }
        for (symbol in 0 until symbolCount) {
            val length = lengths[offset + symbol]
            if (length != 0) symbols[offsets[length]++] = symbol
        }
    }

    /** Entry for the low 15 bits [peek], following one subtable when the primary slot says to. */
    fun lookup(peek: Int): Int {
        var entry = entries[peek and primaryMask]
        if (entry and FLAG_SUBTABLE != 0) {
            val extra = (entry ushr EXTRA_SHIFT) and 0xF
            val extraValue = (peek ushr tableBits) and ((1 shl extra) - 1)
            entry = entries[(entry and PAYLOAD_MASK) + extraValue]
        }
        return entry
    }

    private fun bitReverse(code: Int, length: Int): Int {
        var value = code
        var reversed = 0
        repeat(length) {
            reversed = (reversed shl 1) or (value and 1)
            value = value ushr 1
        }
        return reversed
    }

    private fun fail(message: String): Nothing = throw IllegalArgumentException(message)

    companion object {
        const val MAX_CODE_LENGTH = 15
        const val MAX_SYMBOLS = 288

        const val PAYLOAD_MASK = 0xFFF
        const val LEN_SHIFT = 12
        const val EXTRA_SHIFT = 16
        const val FLAG_SUBTABLE = 1 shl 20
        const val FLAG_INVALID = 1 shl 21
        const val INVALID = FLAG_INVALID

        /** Literal/length primary bits, and the worst-case table size (`enough 288 11 15`). */
        const val LITERAL_TABLE_BITS = 11
        const val LITERAL_ENOUGH = 2342

        /** Distance primary bits (`enough 32 8 15`). */
        const val DISTANCE_TABLE_BITS = 8
        const val DISTANCE_ENOUGH = 402

        /** Code-length codes are at most 7 bits, so the primary table never needs a subtable. */
        const val PRECODE_TABLE_BITS = 7
        const val PRECODE_ENOUGH = 128
    }
}
