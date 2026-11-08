package net.maizegenetics.tassel.io

import kotlinx.io.Source
import kotlinx.io.indexOf
import kotlinx.io.readByteArray

/**
 * One `##contig` entry, in file order. [length] is absent when the line has no integer `length`.
 *
 * A `##contig` line with no id is dropped, which is what [Vcf.read] already did with every meta line.
 */
internal data class VcfContig(val id: String, val length: Long?)

/**
 * The `#CHROM` layout and the contig list from one VCF header.
 *
 * Column indexes match TASSEL 5's `HeaderPositions`: by name, not by order. [fixed] is the number of columns
 * before the first taxon, and the taxa start after whichever of `INFO` and `FORMAT` comes last.
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
internal class VcfHeader(
    val contigs: List<VcfContig>,
    val taxaNames: List<String>,
    val chromosomeColumn: Int,
    val positionColumn: Int,
    val snpIdColumn: Int,
    val referenceColumn: Int,
    val alternateColumn: Int,
    val formatColumn: Int,
    val fixed: Int,
)

/**
 * Reads meta lines and the `#CHROM` line, and leaves [source] at the first record.
 *
 * `##contig=<ID=...,length=...>` keeps its id and optional length, in order. Other keys are ignored, and a quoted
 * value may contain commas. `##SAMPLE` and every other meta line are read and dropped. A `##contig` line that has
 * no id is dropped too, so a file the old reader accepted still loads.
 */
internal fun readVcfHeader(source: Source): VcfHeader {
    val contigs = ArrayList<VcfContig>()
    var line = source.readVcfLine() ?: throw IllegalArgumentException("VCF input is empty")
    while (line.startsWith(META_PREFIX)) {
        contigOrNull(line)?.let { contigs.add(it) }
        line = source.readVcfLine() ?: throw IllegalArgumentException("VCF input has no #CHROM header line")
    }
    return headerFrom(line.decodeToString(), contigs)
}

/**
 * [readVcfHeader] for a reader that yields lines instead of a [Source].
 *
 * BGZF uses this so header parsing suspends with the block reader. Line breaks match [BgzfReader.readLine].
 */
internal suspend fun readVcfHeaderLines(nextLine: suspend () -> String?): VcfHeader {
    val contigs = ArrayList<VcfContig>()
    var line = nextLine() ?: throw IllegalArgumentException("VCF input is empty")
    while (line.startsWith("##")) {
        contigOrNull(line.encodeToByteArray())?.let { contigs.add(it) }
        line = nextLine() ?: throw IllegalArgumentException("VCF input has no #CHROM header line")
    }
    return headerFrom(line, contigs)
}

/**
 * One line, without its `\n` or the `\r` of a `\r\n`. A final line with no break is returned whole, including a
 * bare `\r`, matching [kotlinx.io.readLine]. Returns null when [this] is already exhausted.
 */
private fun Source.readVcfLine(): ByteArray? {
    if (exhausted()) return null
    val newline = indexOf(LINE_FEED)
    if (newline == -1L) return readByteArray()
    val raw = readByteArray(newline.toInt())
    skip(1)
    return if (raw.isNotEmpty() && raw[raw.lastIndex] == CARRIAGE_RETURN) raw.copyOf(raw.size - 1) else raw
}

private fun headerFrom(text: String, contigs: List<VcfContig>): VcfHeader {
    val names = text.split('\t')
    // TASSEL 5 accepts #CHR as well as #CHROM.
    val chromosomeColumn = maxOf(names.indexOf("#CHROM"), names.indexOf("#CHR"))
    val positionColumn = names.indexOf("POS")
    val snpIdColumn = names.indexOf("ID")
    val referenceColumn = names.indexOf("REF")
    val alternateColumn = names.indexOf("ALT")
    val formatColumn = names.indexOf("FORMAT")
    val fixed = maxOf(names.indexOf("INFO"), formatColumn) + 1
    require(chromosomeColumn >= 0 && positionColumn > 0 && referenceColumn > 0 && alternateColumn > 0) {
        "VCF header must name #CHROM, POS, REF, and ALT: ${vcfExcerpt(text)}"
    }
    require(fixed > maxOf(chromosomeColumn, positionColumn, referenceColumn, alternateColumn)) {
        "VCF header must name INFO or FORMAT after the other columns: ${vcfExcerpt(text)}"
    }
    require(names.size > fixed) { "VCF header must name at least one taxon: ${vcfExcerpt(text)}" }
    return VcfHeader(
        contigs = contigs,
        taxaNames = names.subList(fixed, names.size).toList(),
        chromosomeColumn = chromosomeColumn,
        positionColumn = positionColumn,
        snpIdColumn = snpIdColumn,
        referenceColumn = referenceColumn,
        alternateColumn = alternateColumn,
        formatColumn = formatColumn,
        fixed = fixed,
    )
}

/**
 * The contig on a `##contig=<...>` line, or null when this is not one or it has no id.
 * Other keys, including a non-integer `length`, do not drop the entry.
 */
private fun contigOrNull(line: ByteArray): VcfContig? {
    if (!line.startsWith(CONTIG_PREFIX) || line.isEmpty() || line[line.lastIndex] != ANGLE_CLOSE) return null
    val stop = line.size - 1
    var index = CONTIG_PREFIX.size
    var id: String? = null
    var length: Long? = null
    var sawLength = false
    while (index < stop) {
        val equals = line.indexOfByte(index, stop, EQUALS)
        if (equals < 0) break
        val key = line.decodeToString(index, equals)
        index = equals + 1
        val value: String
        if (index < stop && line[index] == QUOTE) {
            val close = line.indexOfByte(index + 1, stop, QUOTE)
            if (close < 0) return null
            value = line.decodeToString(index + 1, close)
            index = close + 1
        } else {
            val comma = line.indexOfByte(index, stop, COMMA)
            val valueEnd = if (comma < 0) stop else comma
            value = line.decodeToString(index, valueEnd)
            index = valueEnd
        }
        if (index < stop && line[index] == COMMA) index++
        when (key) {
            "ID" -> if (id == null && value.isNotEmpty()) id = value
            "length" -> if (!sawLength) {
                sawLength = true
                length = value.toLongOrNull()
            }
        }
    }
    val contigId = id ?: return null
    return VcfContig(contigId, length)
}

private fun ByteArray.startsWith(prefix: ByteArray): Boolean {
    if (size < prefix.size) return false
    for (i in prefix.indices) {
        if (this[i] != prefix[i]) return false
    }
    return true
}

private fun ByteArray.indexOfByte(from: Int, until: Int, target: Byte): Int {
    var index = from
    while (index < until) {
        if (this[index] == target) return index
        index++
    }
    return -1
}

private val META_PREFIX = "##".encodeToByteArray()
private val CONTIG_PREFIX = "##contig=<".encodeToByteArray()

private val LINE_FEED = '\n'.code.toByte()
private val CARRIAGE_RETURN = '\r'.code.toByte()
private val EQUALS = '='.code.toByte()
private val COMMA = ','.code.toByte()
private val QUOTE = '"'.code.toByte()
private val ANGLE_CLOSE = '>'.code.toByte()
