package net.maizegenetics.tassel.association

/*
 * Primitive columns that grow a segment at a time. A segment is allocated when the previous one fills and nothing is
 * ever copied, so growing a column never needs twice its memory, and a row costs exactly its cells.
 */

internal const val COLUMN_SEGMENT_SHIFT = 13
private const val SEGMENT_SIZE = 1 shl COLUMN_SEGMENT_SHIFT
private const val SEGMENT_MASK = SEGMENT_SIZE - 1

internal class DoubleColumn {
    private val segments = ArrayList<DoubleArray>()
    var size = 0
        private set

    fun add(value: Double) {
        val offset = size and SEGMENT_MASK
        if (offset == 0) segments.add(DoubleArray(SEGMENT_SIZE))
        segments[segments.size - 1][offset] = value
        size++
    }

    operator fun get(index: Int): Double = segments[index ushr COLUMN_SEGMENT_SHIFT][index and SEGMENT_MASK]

    fun seal() {
        val used = size and SEGMENT_MASK
        if (used != 0) segments[segments.size - 1] = segments[segments.size - 1].copyOf(used)
    }

    val bytes: Long get() = segments.sumOf { it.size.toLong() * Double.SIZE_BYTES }
}

internal class IntColumn {
    private val segments = ArrayList<IntArray>()
    var size = 0
        private set

    fun add(value: Int) {
        val offset = size and SEGMENT_MASK
        if (offset == 0) segments.add(IntArray(SEGMENT_SIZE))
        segments[segments.size - 1][offset] = value
        size++
    }

    operator fun get(index: Int): Int = segments[index ushr COLUMN_SEGMENT_SHIFT][index and SEGMENT_MASK]

    fun seal() {
        val used = size and SEGMENT_MASK
        if (used != 0) segments[segments.size - 1] = segments[segments.size - 1].copyOf(used)
    }

    val bytes: Long get() = segments.sumOf { it.size.toLong() * Int.SIZE_BYTES }
}

internal class ByteColumn {
    private val segments = ArrayList<ByteArray>()
    var size = 0
        private set

    fun add(value: Int) {
        val offset = size and SEGMENT_MASK
        if (offset == 0) segments.add(ByteArray(SEGMENT_SIZE))
        segments[segments.size - 1][offset] = value.toByte()
        size++
    }

    operator fun get(index: Int): Int = segments[index ushr COLUMN_SEGMENT_SHIFT][index and SEGMENT_MASK].toInt()

    fun seal() {
        val used = size and SEGMENT_MASK
        if (used != 0) segments[segments.size - 1] = segments[segments.size - 1].copyOf(used)
    }

    val bytes: Long get() = segments.sumOf { it.size.toLong() }
}

/** Each distinct string once, numbered in order of first use; the empty string is always 0. */
internal class StringDictionary {
    private val codes = HashMap<String, Int>()
    private val values = ArrayList<String>()

    init {
        code("")
    }

    fun code(value: String): Int = codes.getOrPut(value) {
        values.add(value)
        values.size - 1
    }

    operator fun get(code: Int): String = values[code]

    val size: Int get() = values.size

    /** A JVM-style estimate of the strings' memory: 48 bytes of object and map overhead plus two per character. */
    val estimatedBytes: Long get() = values.sumOf { 48L + 2L * it.length }
}
