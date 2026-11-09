package net.maizegenetics.tassel.io.deflate

/**
 * One BGZF member, enough for a test to inflate it. Framing itself belongs to the gzip source; this only
 * finds the raw DEFLATE slice and the trailer.
 */
internal class BgzfBlock(val file: ByteArray, val from: Int, val to: Int, val uncompressedSize: Int, val crc: Int)

internal fun bgzfBlocks(file: ByteArray): List<BgzfBlock> {
    val blocks = ArrayList<BgzfBlock>()
    var start = 0
    while (start < file.size) {
        if (start + 18 > file.size) error("truncated gzip header at $start")
        if ((file[start].toInt() and 0xFF) != 0x1F || (file[start + 1].toInt() and 0xFF) != 0x8B) {
            error("gzip magic at $start")
        }
        if ((file[start + 2].toInt() and 0xFF) != 8) error("gzip method at $start")
        val flags = file[start + 3].toInt() and 0xFF
        var cursor = start + 10
        var blockSize = -1
        if (flags and 4 != 0) {
            val extraLength = u16(file, cursor)
            cursor += 2
            val extraEnd = cursor + extraLength
            var sub = cursor
            while (sub + 4 <= extraEnd) {
                val id1 = file[sub].toInt() and 0xFF
                val id2 = file[sub + 1].toInt() and 0xFF
                val subLength = u16(file, sub + 2)
                if (id1 == 'B'.code && id2 == 'C'.code && subLength == 2) {
                    blockSize = u16(file, sub + 4)
                }
                sub += 4 + subLength
            }
            cursor = extraEnd
        }
        if (flags and 8 != 0) cursor = skipCString(file, cursor)
        if (flags and 16 != 0) cursor = skipCString(file, cursor)
        if (flags and 2 != 0) cursor += 2
        if (blockSize < 0) error("BGZF block at $start has no BC subfield")
        val end = start + blockSize + 1
        if (end > file.size || end - 8 < cursor) error("BGZF block at $start overruns the file")
        val crc = u32(file, end - 8)
        val uncompressed = u32(file, end - 4)
        if (uncompressed < 0) error("BGZF ISIZE does not fit in a block at $start")
        blocks.add(BgzfBlock(file, cursor, end - 8, uncompressed, crc))
        start = end
    }
    return blocks
}

private fun skipCString(file: ByteArray, from: Int): Int {
    var cursor = from
    while (file[cursor].toInt() != 0) cursor++
    return cursor + 1
}

private fun u16(file: ByteArray, at: Int): Int = (file[at].toInt() and 0xFF) or ((file[at + 1].toInt() and 0xFF) shl 8)

private fun u32(file: ByteArray, at: Int): Int = u16(file, at) or (u16(file, at + 2) shl 16)
