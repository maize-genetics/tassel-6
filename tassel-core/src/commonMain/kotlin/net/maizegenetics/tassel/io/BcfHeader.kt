package net.maizegenetics.tassel.io

import kotlinx.io.Source
import kotlinx.io.readByteArray

/**
 * A BCF 2.2 header: the VCF column layout plus the dictionaries a record indexes by `IDX`.
 *
 * [strings] is the shared FILTER, INFO, and FORMAT dictionary. [contigs] is a separate dictionary. A null slot is
 * an `IDX` gap. [gtId] is the string-dictionary index of `FORMAT/GT`, or null when the header has no GT field.
 * [vcf] is the column layout from [readVcfHeaderLines]. A header with no samples has an empty taxa list.
 */
internal class BcfHeader(val vcf: VcfHeader, val strings: List<String?>, val contigs: List<VcfContig?>, val gtId: Int?)

/**
 * Reads the BCF preamble and header from an inflated byte stream and stops after the trailing NUL.
 *
 * The source is raw BCF, with the magic at byte 0. The caller inflates BGZF. BCF 2.1 and any other magic are
 * refused with htslib 1.25's messages. Bytes after the header, including the records, are left unread.
 */
internal fun readBcfHeader(source: Source): BcfHeader {
    val magic = source.readExact(BCF_MAGIC_LENGTH)
    if (!magic.startsWith(BCF_LETTERS)) throw IllegalArgumentException(INVALID_BCF_MAGIC)
    if (magic[3] != BCF_MAJOR || magic[4] != BCF_MINOR) throw IllegalArgumentException(UNSUPPORTED_BCF_VERSION)
    val textLength = source.readExact(UINT32_BYTES).uint32Le(0)
    if (textLength == 0L || textLength > Int.MAX_VALUE.toLong() || !source.request(textLength)) {
        throw IllegalArgumentException(TRUNCATED_BCF_HEADER)
    }
    val textBytes = source.readByteArray(textLength.toInt())
    if (textBytes[textBytes.lastIndex] != 0.toByte()) {
        throw IllegalArgumentException("BCF header text is missing its trailing NUL")
    }
    val lines = headerLines(textBytes.decodeToString(0, textBytes.size - 1))
    val vcf = readVcfHeaderLines(lines.iterator(), requireTaxon = false)
    val fields = dictionaryLines(lines, STRING_KINDS)
    val strings = stringDictionary(fields)
    val gt = fields.firstOrNull { it.kind == FORMAT && it.id == GT }
    return BcfHeader(
        vcf = vcf,
        strings = strings.ids,
        contigs = contigDictionary(lines, vcf.contigs),
        gtId = if (gt == null) null else strings.indexOf(GT).takeIf { it >= 0 },
    )
}

/**
 * FILTER, INFO, and FORMAT share one dictionary, keyed the way htslib 1.25's `bcf_hdr_set_idx` keys it.
 *
 * `FILTER/PASS` is index 0 when the header has that line, which is htslib inserting PASS before the other
 * lines. An explicit `IDX` occupies that slot and leaves a gap. A line with no `IDX` takes the current length,
 * so it does not fill a gap. The same id keeps the index it already has.
 */
private fun stringDictionary(entries: List<DictLine>): DictionarySlots {
    val slots = DictionarySlots()
    if (entries.any { it.kind == FILTER && it.id == PASS }) slots.place(PASS, 0)
    for (entry in entries) slots.place(entry)
    return slots
}

private fun contigDictionary(lines: List<String>, fileOrder: List<VcfContig>): List<VcfContig?> {
    val lengthOf = HashMap<String, Long?>()
    for (contig in fileOrder) {
        if (contig.id !in lengthOf) lengthOf[contig.id] = contig.length
    }
    val slots = DictionarySlots()
    for (entry in dictionaryLines(lines, CONTIG_KINDS)) slots.place(entry)
    return slots.ids.map { id -> if (id == null) null else VcfContig(id, lengthOf[id]) }
}

private fun dictionaryLines(lines: List<String>, kinds: Set<String>): List<DictLine> {
    val entries = ArrayList<DictLine>()
    for (line in lines) {
        val kind = dictionaryKind(line)
        if (kind == null || kind !in kinds) continue
        val keys = dictionaryKeys(line) ?: continue
        entries.add(DictLine(kind, keys.id, keys.idx))
    }
    return entries
}

private class DictLine(val kind: String, val id: String, val idx: Int?)

private class DictKeys(val id: String, val idx: Int?)

/**
 * `ID` and `IDX` from one structured header line. A quoted value may contain commas. A line with no id, or a
 * quote that never closes, is skipped, which is what [readVcfHeaderLines] does for a `##contig` line.
 */
private fun dictionaryKeys(line: String): DictKeys? {
    if (!line.endsWith('>')) return null
    val open = line.indexOf('<')
    if (open < 0) return null
    var index = open + 1
    val stop = line.length - 1
    var id: String? = null
    var idx: Int? = null
    var sawIdx = false
    while (index < stop) {
        val equals = line.indexOf('=', index)
        if (equals < 0 || equals >= stop) break
        val key = line.substring(index, equals)
        index = equals + 1
        val value: String
        if (index < stop && line[index] == '"') {
            val close = line.indexOf('"', index + 1)
            if (close < 0 || close > stop) return null
            value = line.substring(index + 1, close)
            index = close + 1
        } else {
            var end = index
            while (end < stop && line[end] != ',') end++
            value = line.substring(index, end)
            index = end
        }
        if (index < stop && line[index] == ',') index++
        when (key) {
            "ID" -> if (id == null && value.isNotEmpty()) id = value
            "IDX" -> if (!sawIdx) {
                sawIdx = true
                val parsed = value.toIntOrNull()
                if (parsed == null || parsed < 0) {
                    throw IllegalArgumentException("BCF dictionary IDX is not a valid index: [$value]")
                }
                idx = parsed
            }
        }
    }
    val found = id ?: return null
    return DictKeys(found, if (sawIdx) idx else null)
}

private fun dictionaryKind(line: String): String? = when {
    line.startsWith(FILTER_PREFIX) -> FILTER
    line.startsWith(INFO_PREFIX) -> INFO
    line.startsWith(FORMAT_PREFIX) -> FORMAT
    line.startsWith(CONTIG_PREFIX) -> CONTIG
    else -> null
}

private fun headerLines(text: String): List<String> {
    val lines = ArrayList<String>()
    var start = 0
    var index = 0
    while (index < text.length) {
        if (text[index] == '\n') {
            lines.add(text.substring(start, index).trimEnd('\r'))
            start = index + 1
        }
        index++
    }
    if (start < text.length) lines.add(text.substring(start).trimEnd('\r'))
    return lines
}

private class DictionarySlots {
    val ids = ArrayList<String?>()
    private val indexOfId = HashMap<String, Int>()

    fun place(id: String, index: Int) {
        if (indexOfId.containsKey(id)) return
        while (ids.size <= index) ids.add(null)
        if (ids[index] != null) {
            throw IllegalArgumentException(
                "Conflicting IDX=$index lines in the header dictionary, the new tag is $id",
            )
        }
        ids[index] = id
        indexOfId[id] = index
    }

    fun place(entry: DictLine) {
        val index = entry.idx
        if (index == null) place(entry.id, ids.size) else place(entry.id, index)
    }

    fun indexOf(id: String): Int = indexOfId[id] ?: -1
}

private fun Source.readExact(count: Int): ByteArray {
    val out = ByteArray(count)
    var at = 0
    while (at < count) {
        val n = readAtMostTo(out, at, count)
        if (n <= 0) throw IllegalArgumentException(TRUNCATED_BCF_HEADER)
        at += n
    }
    return out
}

internal fun ByteArray.uint32Le(at: Int): Long = (this[at].toLong() and 0xff) or
    ((this[at + 1].toLong() and 0xff) shl 8) or
    ((this[at + 2].toLong() and 0xff) shl 16) or
    ((this[at + 3].toLong() and 0xff) shl 24)

private fun ByteArray.startsWith(prefix: ByteArray): Boolean {
    if (size < prefix.size) return false
    for (index in prefix.indices) {
        if (this[index] != prefix[index]) return false
    }
    return true
}

internal const val BCF_MAGIC_LENGTH = 5
internal const val UINT32_BYTES = 4
private const val BCF_MAJOR: Byte = 2
private const val BCF_MINOR: Byte = 2
private val BCF_LETTERS = "BCF".encodeToByteArray()

private const val FILTER = "FILTER"
private const val INFO = "INFO"
private const val FORMAT = "FORMAT"
private const val CONTIG = "CONTIG"
private const val PASS = "PASS"
private const val GT = "GT"
private val STRING_KINDS = setOf(FILTER, INFO, FORMAT)
private val CONTIG_KINDS = setOf(CONTIG)
private const val FILTER_PREFIX = "##FILTER=<"
private const val INFO_PREFIX = "##INFO=<"
private const val FORMAT_PREFIX = "##FORMAT=<"
private const val CONTIG_PREFIX = "##contig=<"

internal const val INVALID_BCF_MAGIC = "Invalid BCF2 magic string"
internal const val UNSUPPORTED_BCF_VERSION = "Invalid BCF2 magic string: only BCFv2.2 is supported"
internal const val TRUNCATED_BCF_HEADER = "BCF header is truncated"
