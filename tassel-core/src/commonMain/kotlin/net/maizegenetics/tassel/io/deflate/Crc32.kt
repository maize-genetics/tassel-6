package net.maizegenetics.tassel.io.deflate

/**
 * CRC-32 used by gzip and BGZF (RFC 1952 section 8).
 *
 * Reflected polynomial `0xEDB88320`, initial value 0, final xor `0xFFFFFFFF`. The tables are built once
 * and [update] slices the input eight bytes at a time. A caller owns one instance and reuses it across
 * members; [update] does not allocate.
 */
internal class Crc32 {
    private var crc: Int = EMPTY

    /** Current checksum, after the final xor. */
    val value: Int
        get() = crc xor EMPTY

    fun reset() {
        crc = EMPTY
    }

    /** Folds `data[from, to)` into the checksum. */
    fun update(data: ByteArray, from: Int, to: Int) {
        if (from < 0 || to < from || to > data.size) {
            throw IllegalArgumentException("invalid crc range")
        }
        var c = crc
        var index = from
        val tables = TABLES
        while (index + 8 <= to) {
            val b0 = (data[index].toInt() xor c) and 0xFF
            val b1 = ((data[index + 1].toInt() and 0xFF) xor (c ushr 8)) and 0xFF
            val b2 = ((data[index + 2].toInt() and 0xFF) xor (c ushr 16)) and 0xFF
            val b3 = ((data[index + 3].toInt() and 0xFF) xor (c ushr 24)) and 0xFF
            val b4 = data[index + 4].toInt() and 0xFF
            val b5 = data[index + 5].toInt() and 0xFF
            val b6 = data[index + 6].toInt() and 0xFF
            val b7 = data[index + 7].toInt() and 0xFF
            c = tables[7][b0] xor tables[6][b1] xor tables[5][b2] xor tables[4][b3] xor
                tables[3][b4] xor tables[2][b5] xor tables[1][b6] xor tables[0][b7]
            index += 8
        }
        while (index < to) {
            c = tables[0][(c xor data[index].toInt()) and 0xFF] xor (c ushr 8)
            index++
        }
        crc = c
    }

    private companion object {
        const val EMPTY = -1 // 0xFFFFFFFF, the reflected initializer and the final xor

        val TABLES: Array<IntArray> = buildTables()

        private fun buildTables(): Array<IntArray> {
            val polynomial = 0xEDB88320.toInt()
            val tables = Array(8) { IntArray(256) }
            for (i in 0..255) {
                var entry = i
                repeat(8) {
                    entry = if (entry and 1 != 0) (entry ushr 1) xor polynomial else entry ushr 1
                }
                tables[0][i] = entry
            }
            for (i in 0..255) {
                var entry = tables[0][i]
                for (slice in 1..7) {
                    entry = tables[0][entry and 0xFF] xor (entry ushr 8)
                    tables[slice][i] = entry
                }
            }
            return tables
        }
    }
}
