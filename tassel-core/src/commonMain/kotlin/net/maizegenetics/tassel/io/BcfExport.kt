package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.position.Position

/**
 * BCF 2.2 bytes for the records [Vcf.write] prints.
 *
 * The header is [VCF_EXPORT_HEADER] plus the three things htslib adds when it encodes that VCF: the PASS
 * filter line, `IDX` on each FORMAT and INFO line, and one `##contig` line per chromosome. A repeated
 * dictionary id keeps the index it already has, which is why both `DP` lines are `IDX=3`. Contig lines
 * have an id and an `IDX` and no length.
 *
 * A record is the layout [readBcfRecord] decodes. QUAL is missing, FILTER is PASS, INFO is empty, and GT
 * is the only FORMAT field. GT is int8: an unphased allele index `a` is stored as `(a + 1) shl 1`, and a
 * missing call is `0`, which is the value htslib writes for `.`.
 *
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
 */
internal class BcfExportHeader(val preamble: ByteArray, val contigs: Map<String, Int>, val gtId: Int)

/** Header preamble for [table], including the magic, the text length, and the trailing NUL. */
internal fun bcfExportHeader(table: GenotypeTable): BcfExportHeader {
    val contigs = LinkedHashMap<String, Int>()
    for (site in 0 until table.numberOfSites) {
        val name = table.positions[site].chromosome.name
        if (name !in contigs) contigs[name] = contigs.size
    }

    val indexOf = HashMap<String, Int>()
    indexOf[PASS] = 0
    var next = 1
    var gtId = -1
    fun assign(id: String): Int {
        val existing = indexOf[id]
        if (existing != null) return existing
        val assigned = next
        next += 1
        indexOf[id] = assigned
        if (id == GT) gtId = assigned
        return assigned
    }

    val lines = ArrayList<String>(VCF_EXPORT_HEADER.size + 1 + contigs.size)
    for (line in VCF_EXPORT_HEADER) {
        if (line.startsWith(FILEFORMAT_PREFIX)) {
            lines.add(line)
            lines.add(PASS_LINE)
            continue
        }
        val id = dictionaryId(line)
        if (id == null) {
            lines.add(line)
        } else {
            lines.add(line.dropLast(1) + ",IDX=${assign(id)}>")
        }
    }
    check(gtId >= 0) { "VCF export header has no GT field" }
    for ((name, index) in contigs) {
        lines.add("##contig=<ID=$name,IDX=$index>")
    }

    val text = buildString {
        for (line in lines) {
            append(line)
            append('\n')
        }
        append(Vcf.FIXED_COLUMNS.joinToString("\t"))
        for (taxon in table.taxa) {
            append('\t')
            append(taxon.name.trim())
        }
        append('\n')
    }
    return BcfExportHeader(preamble(text), contigs, gtId)
}

/**
 * One BCF record, rewritten in place for each site.
 *
 * [bytes] may be longer than [size]. The caller writes `bytes[0, size)`.
 */
internal class BcfRecordEncoder(
    private val sampleCount: Int,
    private val gtId: Int,
    private val contigs: Map<String, Int>,
) {
    private var data = ByteArray(64 + sampleCount * 2)
    var size: Int = 0
        private set

    val bytes: ByteArray get() = data

    init {
        require(sampleCount in 0..0xFFFFFF) { "BCF export sample count does not fit in 24 bits: $sampleCount" }
        require(gtId in 0..127) { "BCF export GT index does not fit in an int8 dictionary key: $gtId" }
    }

    fun encode(position: Position, record: VcfSiteRecord) {
        require(record.calls.size == sampleCount * 2) {
            "BCF export has ${record.calls.size} calls, expected ${sampleCount * 2}"
        }
        require(record.alleleCount in 1..0xFFFF) {
            "BCF export allele count does not fit in 16 bits: ${record.alleleCount}"
        }
        val contig = contigs[position.chromosome.name]
            ?: throw IllegalArgumentException("BCF export has no contig for ${position.chromosome.name}")

        size = 0
        val prefix = size
        int(0)
        int(0)
        val sharedAt = size
        int(contig)
        int(position.position - 1)
        val reference = record.reference.encodeToByteArray()
        int(reference.size)
        int(BCF_FLOAT_MISSING)
        int(record.alleleCount shl 16)
        int((1 shl 24) or sampleCount)
        writeId(position.name)
        writeChars(reference)
        for (alternate in record.alternates) writeChars(alternate.encodeToByteArray())
        require(1 + record.alternates.size == record.alleleCount) {
            "BCF export wrote ${1 + record.alternates.size} alleles, expected ${record.alleleCount}"
        }
        byte((1 shl 4) or TYPE_INT8)
        byte(0)
        val sharedLength = size - sharedAt

        val indivAt = size
        byte((1 shl 4) or TYPE_INT8)
        byte(gtId)
        byte((2 shl 4) or TYPE_INT8)
        for (call in record.calls) byte(gtValue(call))
        val indivLength = size - indivAt

        putInt(prefix, sharedLength)
        putInt(prefix + UINT32_BYTES, indivLength)
    }

    private fun writeId(name: String) {
        if (name.isEmpty() || name == ".") writeChars(ByteArray(0)) else writeChars(name.encodeToByteArray())
    }

    /** A character vector. A count of 15 or more is the following typed integer, which is how [BcfCursor] reads it. */
    private fun writeChars(value: ByteArray) {
        if (value.size < OVERFLOW_COUNT) {
            byte((value.size shl 4) or TYPE_CHAR)
        } else {
            byte((OVERFLOW_COUNT shl 4) or TYPE_CHAR)
            writeTypedSize(value.size)
        }
        writeBytes(value)
    }

    private fun writeTypedSize(count: Int) {
        when {
            count <= 127 -> {
                byte((1 shl 4) or TYPE_INT8)
                byte(count)
            }
            count <= 32767 -> {
                byte((1 shl 4) or TYPE_INT16)
                byte(count)
                byte(count ushr 8)
            }
            else -> {
                byte((1 shl 4) or TYPE_INT32)
                int(count)
            }
        }
    }

    private fun gtValue(index: Byte): Int {
        if (index == VcfSiteRecord.MISSING) return 0
        val allele = index.toInt()
        require(allele >= 0) { "BCF export allele index is negative: $allele" }
        val stored = (allele + 1) shl 1
        require(stored <= 127) { "BCF export allele index does not fit in an int8 genotype: $allele" }
        return stored
    }

    private fun byte(value: Int) {
        ensure(1)
        data[size] = value.toByte()
        size += 1
    }

    private fun int(value: Int) {
        ensure(UINT32_BYTES)
        putInt(size, value)
        size += UINT32_BYTES
    }

    private fun writeBytes(value: ByteArray) {
        ensure(value.size)
        value.copyInto(data, size)
        size += value.size
    }

    private fun putInt(at: Int, value: Int) {
        data[at] = value.toByte()
        data[at + 1] = (value ushr 8).toByte()
        data[at + 2] = (value ushr 16).toByte()
        data[at + 3] = (value ushr 24).toByte()
    }

    private fun ensure(extra: Int) {
        val needed = size + extra
        if (needed <= data.size) return
        var capacity = data.size
        while (capacity < needed) capacity *= 2
        data = data.copyOf(capacity)
    }
}

private fun preamble(text: String): ByteArray {
    val encoded = text.encodeToByteArray()
    val length = encoded.size + 1
    val out = ByteArray(BCF_MAGIC_LENGTH + UINT32_BYTES + length)
    out[0] = 'B'.code.toByte()
    out[1] = 'C'.code.toByte()
    out[2] = 'F'.code.toByte()
    out[3] = 2
    out[4] = 2
    out[BCF_MAGIC_LENGTH] = length.toByte()
    out[BCF_MAGIC_LENGTH + 1] = (length ushr 8).toByte()
    out[BCF_MAGIC_LENGTH + 2] = (length ushr 16).toByte()
    out[BCF_MAGIC_LENGTH + 3] = (length ushr 24).toByte()
    encoded.copyInto(out, BCF_MAGIC_LENGTH + UINT32_BYTES)
    return out
}

/** `ID` from a FORMAT or INFO line. Any other export line, including `##Tassel`, is not a dictionary entry. */
private fun dictionaryId(line: String): String? {
    if (!line.startsWith(FORMAT_PREFIX) && !line.startsWith(INFO_PREFIX)) return null
    require(line.endsWith('>')) { "VCF export dictionary line does not end with '>': $line" }
    val idAt = line.indexOf("ID=")
    require(idAt > 0) { "VCF export dictionary line has no ID: $line" }
    val start = idAt + 3
    var end = start
    while (end < line.length && line[end] != ',' && line[end] != '>') end++
    require(end > start) { "VCF export dictionary line has an empty ID: $line" }
    return line.substring(start, end)
}

private const val FILEFORMAT_PREFIX = "##fileformat="
private const val FORMAT_PREFIX = "##FORMAT=<"
private const val INFO_PREFIX = "##INFO=<"
private const val PASS_LINE = "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>"
private const val PASS = "PASS"
private const val GT = "GT"

/** Typed-value type codes from VCF 4.3 section 6. [BcfCursor] reads these same codes. */
private const val TYPE_INT8 = 1
private const val TYPE_INT16 = 2
private const val TYPE_INT32 = 3
private const val TYPE_CHAR = 7

/** A count nibble of 15 means the next typed integer is the real count. */
private const val OVERFLOW_COUNT = 15
