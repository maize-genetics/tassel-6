package net.maizegenetics.tassel.io.index

import kotlinx.io.Source
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.decompressing

/**
 * One BGZF virtual-offset interval.
 *
 * [start] and [end] are unsigned 64-bit virtual offsets held in a [Long] (SAMv1 section 4.1). A value at or
 * above 2^63 is ordered with [unsignedCompare], not with signed [Long] comparison. htslib's offsets for the
 * files this reader opens are below 2^63; a bin number above 2^31 is held the same way, in a [Long], because
 * a bin is an unsigned 32-bit value.
 */
data class BgzfChunk(val start: Long, val end: Long)

/**
 * One reference sequence in index order.
 *
 * [mapped] and [unmapped] come from the pseudo-bin, whose number is `binLimit + 1` with
 * `binLimit = ((1 << (depth + 1) * 3) - 1) / 7` (37450 for TBI). That bin is not data: its first chunk is
 * the reference's first and last file offsets, and its second chunk is these counts.
 */
data class IndexedContig(val name: String, val mapped: Long, val unmapped: Long)

/**
 * A CSI or TBI index.
 *
 * Coordinates on [chunks] are 0-based and half-open, as in htslib's `hts_itr_query`. A caller that has a
 * 1-based inclusive range converts it before asking. [maximumCoordinate] is `2^(minShift + 3 * depth)`:
 * 2^29 for TBI, and whatever the file's `min_shift` and `depth` give for CSI. Neither is assumed.
 *
 * @see net.maizegenetics.tassel.io.index.readGenomicIndex
 */
class GenomicIndex internal constructor(
    val minShift: Int,
    val depth: Int,
    val maximumCoordinate: Long,
    val contigs: List<IndexedContig>,
    /** Records that belong to no contig, the optional trailing `n_no_coor`. Zero when the field is absent. */
    val unplaced: Long,
    private val references: List<Reference>,
    private val names: Map<String, Int>,
) {
    /** Index of [name], or -1 when the index has no contig by that exact name. */
    fun contigIndex(name: String): Int = names[name] ?: -1

    /**
     * BGZF chunks of [contig] that can hold records overlapping `[beg], [end)`.
     *
     * The chunks are sorted, contained chunks are dropped, overlaps are trimmed, and chunks that meet
     * inside one compressed block are merged. That is the list `hts_itr_query` walks.
     *
     * An empty list is a contig with no bins, or a region that starts at or past [maximumCoordinate].
     *
     * @throws IllegalArgumentException when [contig] is out of range, [beg] is negative, or [end] is before [beg].
     */
    fun chunks(contig: Int, beg: Long, end: Long): List<BgzfChunk> {
        if (contig !in contigs.indices) {
            throw IllegalArgumentException("contig index $contig is outside 0..${contigs.lastIndex}")
        }
        if (beg < 0) throw IllegalArgumentException("region start is negative: $beg")
        if (end < beg) throw IllegalArgumentException("region end $end is before start $beg")
        val reference = references[contig]
        if (reference.bins.isEmpty() || beg >= maximumCoordinate) return emptyList()
        if (beg == end) return emptyList()

        val minOff = minimumOffset(reference, beg)
        val maxOff = maximumOffset(reference, end)
        val queryEnd = if (end > maximumCoordinate) maximumCoordinate else end
        val collected = ArrayList<BgzfChunk>()
        // Walk bins the index actually stores. reg2bins would list every bin in the coordinate
        // space, and a depth-8 CSI queried over the whole contig is 2^24 leaves.
        for ((bin, stored) in reference.bins) {
            if (!binOverlapsQuery(bin, beg, queryEnd)) continue
            for (chunk in stored.chunks) {
                if (unsignedLess(minOff, chunk.end) && unsignedLess(chunk.start, maxOff)) {
                    val start = if (unsignedLess(chunk.start, minOff)) minOff else chunk.start
                    val stop = if (unsignedLess(maxOff, chunk.end)) maxOff else chunk.end
                    collected.add(BgzfChunk(start, stop))
                }
            }
        }
        return mergeChunks(collected)
    }

    /**
     * The lowest virtual offset that can still hold a record at [beg].
     *
     * Starts at the lowest-level bin containing [beg] and walks to earlier siblings, then to parents,
     * until a bin exists. TBI then raises that offset to the linear index.
     */
    private fun minimumOffset(reference: Reference, beg: Long): Long {
        var bin = binFirst(depth) + (beg ushr minShift)
        var found: StoredBin? = null
        while (true) {
            val hit = reference.bins[bin]
            if (hit != null) {
                found = hit
                break
            }
            if (bin == 0L) break
            val firstSibling = (binParent(bin) shl 3) + 1
            bin = if (bin > firstSibling) bin - 1 else binParent(bin)
        }
        if (bin == 0L) found = reference.bins[0L]
        var minOff = found?.loff ?: 0L
        val linear = reference.linear
        val window = beg ushr minShift
        if (linear.isNotEmpty() && window <= Int.MAX_VALUE && window < linear.size.toLong()) {
            val at = window.toInt()
            if (unsignedLess(minOff, linear[at])) minOff = linear[at]
            if (reference.unmappedPresent) {
                var earlier = at - 1
                while (earlier >= 0) {
                    if (unsignedLess(linear[earlier], minOff)) {
                        minOff = linear[earlier]
                        break
                    }
                    earlier--
                }
                val loff = found?.loff ?: 0L
                if (found != null && (unsignedLess(minOff, loff) || earlier < 0)) minOff = loff
            }
        }
        return minOff
    }

    /** True when [bin]'s interval overlaps `[beg], [endExclusive)`, the same set [reg2bins] names. */
    private fun binOverlapsQuery(bin: Long, beg: Long, endExclusive: Long): Boolean {
        val level = binLevel(bin)
        if (level > depth) return false
        val shift = minShift + 3 * (depth - level)
        if (shift !in 0..62) return false
        val offset = bin - binFirst(level)
        if (offset < 0 || offset > (Long.MAX_VALUE ushr shift)) return false
        val start = offset shl shift
        val stop = start + (1L shl shift)
        return start < endExclusive && stop > beg
    }

    /** A virtual offset from the first existing bin to the right of [end], or every remaining byte. */
    private fun maximumOffset(reference: Reference, end: Long): Long {
        if (end > maximumCoordinate) return UNSIGNED_MAX
        var bin = binFirst(depth) + ((end - 1) ushr minShift) + 1
        if (bin >= binLimit(depth)) bin = 0
        while (true) {
            while (bin % 8L == 1L) bin = binParent(bin)
            if (bin == 0L) return UNSIGNED_MAX
            val hit = reference.bins[bin]
            if (hit != null && hit.chunks.isNotEmpty()) return hit.chunks[0].start
            bin++
        }
    }
}

internal class StoredBin(val loff: Long, val chunks: List<BgzfChunk>)

internal class Reference(
    val bins: Map<Long, StoredBin>,
    val linear: LongArray,
    val unmappedPresent: Boolean,
    val mapped: Long,
    val unmapped: Long,
)

/**
 * Reads a CSI or TBI index from [source].
 *
 * The index is BGZF-compressed and small, so it is inflated whole. [source] is not closed.
 *
 * @throws IllegalArgumentException when the magic is wrong, the index is truncated, the tabix
 *   configuration names an unknown format, or a CSI index has no contig names.
 */
fun readGenomicIndex(source: Source): GenomicIndex {
    val bytes = source.decompressing().source.readByteArray()
    if (bytes.size < 4) throw IllegalArgumentException("truncated genomic index")
    return when {
        isMagic(bytes, 'C', 'S', 'I') -> parseCsi(bytes)
        isMagic(bytes, 'T', 'B', 'I') -> parseTbi(bytes)
        else -> throw IllegalArgumentException("genomic index magic is not CSI or TBI")
    }
}

private fun parseCsi(bytes: ByteArray): GenomicIndex {
    val cursor = Cursor(bytes)
    cursor.skip(4)
    val minShift = cursor.i32()
    val depth = cursor.i32()
    val auxLength = cursor.i32()
    if (auxLength < 0) throw IllegalArgumentException("truncated genomic index")
    val aux = cursor.take(auxLength)
    val references = cursor.i32()
    val names = contigNames(aux)
    if (names.size != references) {
        throw IllegalArgumentException("genomic index contig names do not match its references")
    }
    return assemble(minShift, depth, names, readReferences(cursor, references, csi = true, depth), readUnplaced(cursor))
}

private fun parseTbi(bytes: ByteArray): GenomicIndex {
    val cursor = Cursor(bytes)
    cursor.skip(4)
    val references = cursor.i32()
    val header = cursor.take(TABIX_HEADER)
    val nameBytes = Cursor(header).let { names ->
        names.skip(24)
        names.u32()
    }
    if (nameBytes > Int.MAX_VALUE) throw IllegalArgumentException("truncated genomic index")
    val names = contigNames(header + cursor.take(nameBytes.toInt()))
    if (names.size != references) {
        throw IllegalArgumentException("genomic index contig names do not match its references")
    }
    val depth = TBI_DEPTH
    val bins = readReferences(cursor, references, csi = false, depth)
    return assemble(TBI_MIN_SHIFT, depth, names, bins, readUnplaced(cursor))
}

private fun assemble(
    minShift: Int,
    depth: Int,
    names: List<String>,
    references: List<Reference>,
    unplaced: Long,
): GenomicIndex {
    if (minShift !in 0..28 || depth !in 1..10 || minShift + 3 * depth !in 1..62) {
        throw IllegalArgumentException("genomic index geometry min_shift=$minShift depth=$depth is not usable")
    }
    val contigs = names.mapIndexed { index, name ->
        val mapped = references[index].mapped
        val unmapped = references[index].unmapped
        IndexedContig(name, mapped, unmapped)
    }
    val lookup = HashMap<String, Int>(names.size)
    names.forEachIndexed { index, name ->
        if (lookup.put(name, index) != null) {
            throw IllegalArgumentException("genomic index repeats contig $name")
        }
    }
    return GenomicIndex(
        minShift,
        depth,
        1L shl (minShift + 3 * depth),
        contigs,
        unplaced,
        references,
        lookup,
    )
}

private fun readReferences(cursor: Cursor, count: Int, csi: Boolean, depth: Int): List<Reference> {
    if (count < 0) throw IllegalArgumentException("truncated genomic index")
    val pseudo = binLimit(depth) + 1
    return List(count) {
        val binCount = cursor.i32()
        if (binCount < 0) throw IllegalArgumentException("truncated genomic index")
        val bins = HashMap<Long, StoredBin>(binCount)
        var mapped = 0L
        var unmapped = 0L
        var unmappedPresent = true
        repeat(binCount) {
            val bin = cursor.u32()
            val loff = if (csi) cursor.u64() else 0L
            val chunks = readChunks(cursor)
            if (bin == pseudo) {
                if (chunks.size > 1) {
                    mapped = chunks[1].start
                    unmapped = chunks[1].end
                    unmappedPresent = unmapped != 0L
                }
            } else if (bins.put(bin, StoredBin(loff, chunks)) != null) {
                throw IllegalArgumentException("genomic index repeats bin $bin")
            }
        }
        val linear = if (csi) {
            LongArray(0)
        } else {
            val windows = cursor.i32()
            if (windows < 0) throw IllegalArgumentException("truncated genomic index")
            LongArray(windows) { cursor.u64() }.also { fillLinearGaps(it) }
        }
        val resolved = if (csi) bins else linearOffsets(bins, linear, depth)
        Reference(resolved, linear, unmappedPresent, mapped, unmapped)
    }
}

/** TBI stores no per-bin loffset. htslib fills it from the linear index after the gaps are closed. */
private fun linearOffsets(bins: Map<Long, StoredBin>, linear: LongArray, depth: Int): Map<Long, StoredBin> {
    val limit = binLimit(depth)
    return bins.mapValues { (bin, stored) ->
        if (bin >= limit) {
            stored
        } else {
            val bottom = binBottom(bin, depth)
            StoredBin(if (bottom in linear.indices) linear[bottom] else 0L, stored.chunks)
        }
    }
}

private fun binBottom(bin: Long, depth: Int): Int {
    val level = binLevel(bin)
    if (level > depth) return Int.MAX_VALUE
    val shifted = (bin - binFirst(level)) shl ((depth - level) * 3)
    if (shifted < 0L || shifted > Int.MAX_VALUE) return Int.MAX_VALUE
    return shifted.toInt()
}

private fun binLevel(bin: Long): Int {
    var level = 0
    var value = bin
    while (value != 0L && level < 32) {
        level++
        value = binParent(value)
    }
    return level
}

/**
 * Tabix `format`, `col_seq`, `col_beg`, `col_end`, `meta`, `skip`, and `l_nm`, then the contig names.
 *
 * A CSI written for BCF has no tabix configuration. That is rejected: this reader names contigs from `aux`.
 */
private fun contigNames(aux: ByteArray): List<String> {
    if (aux.size < TABIX_HEADER) throw IllegalArgumentException("genomic index has no contig names")
    val cursor = Cursor(aux)
    val format = cursor.i32()
    cursor.skip(20)
    val nameBytes = cursor.u32()
    if ((format and 0xFFFF) !in 0..3) throw IllegalArgumentException("unknown tabix format $format")
    if (nameBytes <= 0L || nameBytes > (aux.size - TABIX_HEADER).toLong()) {
        throw IllegalArgumentException("genomic index has no contig names")
    }
    val names = ArrayList<String>()
    var at = TABIX_HEADER
    val end = TABIX_HEADER + nameBytes.toInt()
    while (at < end) {
        var stop = at
        while (stop < end && aux[stop] != 0.toByte()) stop++
        if (stop == end) throw IllegalArgumentException("truncated genomic index")
        names.add(aux.decodeToString(at, stop))
        at = stop + 1
    }
    return names
}

private fun readChunks(cursor: Cursor): List<BgzfChunk> {
    val count = cursor.i32()
    if (count < 0) throw IllegalArgumentException("truncated genomic index")
    return List(count) { BgzfChunk(cursor.u64(), cursor.u64()) }
}

/** Leading empty windows stay empty. A later hole copies the next offset, as htslib does for older tabix. */
private fun fillLinearGaps(linear: LongArray) {
    var first = 0
    while (first < linear.size && linear[first] == 0L) first++
    var index = linear.size - 1
    while (index > first) {
        if (linear[index - 1] == 0L) linear[index - 1] = linear[index]
        index--
    }
}

private fun readUnplaced(cursor: Cursor): Long = if (cursor.remaining() >= 8) cursor.u64() else 0L

private fun mergeChunks(chunks: List<BgzfChunk>): List<BgzfChunk> {
    if (chunks.isEmpty()) return emptyList()
    val sorted = chunks.sortedWith { left, right ->
        val byStart = unsignedCompare(left.start, right.start)
        if (byStart != 0) byStart else unsignedCompare(left.end, right.end)
    }
    val kept = ArrayList<BgzfChunk>()
    for (chunk in sorted) {
        if (kept.isEmpty() || unsignedLess(kept.last().end, chunk.end)) kept.add(chunk)
    }
    for (index in 1 until kept.size) {
        val previous = kept[index - 1]
        val current = kept[index]
        if (!unsignedLess(previous.end, current.start)) {
            kept[index - 1] = BgzfChunk(previous.start, current.start)
        }
    }
    val merged = ArrayList<BgzfChunk>()
    for (chunk in kept) {
        val last = merged.lastOrNull()
        if (last != null && (last.end ushr 16) == (chunk.start ushr 16)) {
            merged[merged.lastIndex] = BgzfChunk(last.start, chunk.end)
        } else {
            merged.add(chunk)
        }
    }
    return merged
}

private fun binLimit(depth: Int): Long = ((1L shl ((depth + 1) * 3)) - 1) / 7

/**
 * Smallest bin containing `[beg], [end)`.
 *
 * This is the generalized `reg2bin` from the CSI specification: `min_shift` and `depth` come from the file.
 * [end] is exclusive. The TBI form is `min_shift` 14 and `depth` 5.
 */
internal fun reg2bin(beg: Long, endExclusive: Long, minShift: Int, depth: Int): Long {
    var end = endExclusive - 1
    var shift = minShift
    var base = ((1L shl ((depth shl 1) + depth)) - 1) / 7
    var level = depth
    while (level > 0) {
        if ((beg ushr shift) == (end ushr shift)) return base + (beg ushr shift)
        level--
        shift += 3
        base -= 1L shl ((level shl 1) + level)
    }
    return 0
}

/**
 * Every bin overlapping `[beg], [end)`, from the root down to the leaves.
 *
 * [end] is clipped to the geometry's limit. An empty range returns no bins. Bins the index does not
 * contain are still listed; a query skips those.
 */
internal fun reg2bins(beg: Long, endExclusive: Long, minShift: Int, depth: Int): List<Long> {
    val limit = 1L shl (minShift + 3 * depth)
    var end = endExclusive
    if (end >= limit) end = limit
    if (beg >= end) return emptyList()
    val last = end - 1
    var shift = minShift + 3 * depth
    var base = 0L
    val bins = ArrayList<Long>()
    for (level in 0..depth) {
        val from = base + (beg ushr shift)
        val to = base + (last ushr shift)
        var bin = from
        while (bin <= to) {
            bins.add(bin)
            bin++
        }
        base += 1L shl ((level shl 1) + level)
        shift -= 3
    }
    return bins
}

internal fun binFirst(level: Int): Long = ((1L shl ((level shl 1) + level)) - 1) / 7

internal fun binParent(bin: Long): Long = (bin - 1) ushr 3

/** Signed [Long] ordering with the high bit flipped, which is unsigned 64-bit order. */
internal fun unsignedCompare(left: Long, right: Long): Int =
    (left xor Long.MIN_VALUE).compareTo(right xor Long.MIN_VALUE)

internal fun unsignedLess(left: Long, right: Long): Boolean = unsignedCompare(left, right) < 0

private const val TABIX_HEADER = 28
private const val TBI_MIN_SHIFT = 14
private const val TBI_DEPTH = 5
private const val UNSIGNED_MAX = -1L

private fun isMagic(bytes: ByteArray, first: Char, second: Char, third: Char): Boolean =
    bytes[0] == first.code.toByte() &&
        bytes[1] == second.code.toByte() &&
        bytes[2] == third.code.toByte() &&
        bytes[3] == 1.toByte()

private class Cursor(private val bytes: ByteArray) {
    private var at = 0

    fun remaining(): Int = bytes.size - at

    fun skip(count: Int) {
        need(count)
        at += count
    }

    fun take(count: Int): ByteArray {
        need(count)
        val slice = bytes.copyOfRange(at, at + count)
        at += count
        return slice
    }

    fun i32(): Int {
        need(4)
        val value = (bytes[at].toInt() and 0xFF) or
            ((bytes[at + 1].toInt() and 0xFF) shl 8) or
            ((bytes[at + 2].toInt() and 0xFF) shl 16) or
            ((bytes[at + 3].toInt() and 0xFF) shl 24)
        at += 4
        return value
    }

    fun u32(): Long {
        need(4)
        val value = (bytes[at].toLong() and 0xFF) or
            ((bytes[at + 1].toLong() and 0xFF) shl 8) or
            ((bytes[at + 2].toLong() and 0xFF) shl 16) or
            ((bytes[at + 3].toLong() and 0xFF) shl 24)
        at += 4
        return value
    }

    fun u64(): Long {
        need(8)
        var value = 0L
        for (index in 0 until 8) {
            value = value or ((bytes[at + index].toLong() and 0xFF) shl (8 * index))
        }
        at += 8
        return value
    }

    private fun need(count: Int) {
        if (count < 0 || at > bytes.size - count) throw IllegalArgumentException("truncated genomic index")
    }
}
