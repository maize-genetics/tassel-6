package net.maizegenetics.tassel.io

/**
 * Reference span of one VCF record, in htslib's 0-based half-open coordinates.
 *
 * [pos] is the 1-based `POS` column. [beg] is `POS - 1` and [end] is exclusive. This is `tbx_parse1` for
 * `TBX_VCF` in htslib 1.25, the same rule as `get_rlen` in `vcf.c` (ADR-0021).
 */
internal data class VcfSpan(val pos: Long, val beg: Long, val end: Long)

/** True when this span overlaps [region]. An open-ended region runs through the end of the contig. */
internal fun VcfSpan.overlaps(region: VcfRegion): Boolean {
    val queryBeg = (region.start ?: 1L) - 1L
    val queryEnd = region.end ?: Long.MAX_VALUE
    return beg < queryEnd && end > queryBeg
}

/**
 * The reference span htslib uses when it decides that a record overlaps a region.
 *
 * The span starts at `POS` and covers the longest of the `REF` allele, a valid `INFO/END`, the
 * absolute `INFO/SVLEN` of a symbolic `<DEL>`, `<DUP>`, `<CNV>`, or `<INV>` allele (and a subtype
 * such as `<DEL:ME>`), and `FORMAT/LEN` when an allele is `<*>` or `<NON_REF>`.
 */
internal fun vcfReferenceSpan(line: String): VcfSpan {
    val fields = TabFields(line)
    val posText = fields.text(1)
        ?: throw IllegalArgumentException("VCF record is missing POS: ${vcfExcerpt(line)}")
    val pos = posText.toLongOrNull()
        ?: throw IllegalArgumentException("VCF position must be an integer: [${posText.trim()}]")
    var beg = pos
    var end = beg
    beg -= 1
    if (beg < 0) beg = 0
    if (end < 1) end = 1
    val ref = fields.text(3) ?: throw IllegalArgumentException("VCF record is missing REF: ${vcfExcerpt(line)}")
    val reflen = ref.length.toLong()
    if (ref.isNotEmpty()) end = beg + reflen

    val alt = fields.text(4).orEmpty()
    val alleles = if (alt.isEmpty()) emptyList() else alt.split(',')
    val useSvlen = BooleanArray(alleles.size)
    var readFormatLength = false
    for (index in alleles.indices) {
        val allele = alleles[index]
        if (svlenOnRef(allele)) {
            useSvlen[index] = true
        } else if (allele == "<*>" || allele == "<NON_REF>") {
            readFormatLength = true
        }
    }

    var svlen = 0L
    val info = fields.text(7)
    if (info != null && info != ".") {
        val endText = infoValue(info, "END")
        if (endText != null && endText != ".") {
            val infoEnd = endText.toLongOrNull()
            if (infoEnd != null && infoEnd > beg) end = infoEnd
        }
        val svText = infoValue(info, "SVLEN")
        if (svText != null) {
            val parts = svText.split(',')
            val limit = minOf(alleles.size, parts.size)
            for (index in 0 until limit) {
                val length = if (useSvlen[index]) absoluteLong(parts[index]) else 1L
                if (length > svlen) svlen = length
            }
        }
    }

    var formatLength = 0L
    if (readFormatLength) {
        val format = fields.text(8)
        if (format != null) {
            val lenAt = format.split(':').indexOf("LEN")
            if (lenAt >= 0) {
                for (sample in fields.samples()) {
                    val value = colonField(sample, lenAt) ?: break
                    val length = value.toLongOrNull() ?: 0L
                    if (length > formatLength) formatLength = length
                }
            }
        }
    }

    val fromLengths = beg + maxOf(reflen, svlen, formatLength)
    if (fromLengths > end) end = fromLengths
    return VcfSpan(pos, beg, end)
}

/** htslib `svlen_on_ref_for_vcf_alt`: symbolic CNV, DEL, DUP, or INV, including a subtype before `>`. */
private fun svlenOnRef(alt: String): Boolean {
    if (alt.length < 5 || alt[0] != '<') return false
    if (alt[4] != '>' && alt[4] != ':') return false
    if (!alt.startsWith("<CNV") && !alt.startsWith("<DEL") && !alt.startsWith("<DUP") && !alt.startsWith("<INV")) {
        return false
    }
    return alt.last() == '>'
}

private fun absoluteLong(text: String): Long {
    val value = text.toLongOrNull() ?: return 0L
    if (value >= 0) return value
    if (value == Long.MIN_VALUE) return Long.MAX_VALUE
    return -value
}

/** The value of `key=` at the start of [info] or after a semicolon, without the trailing fields. */
private fun infoValue(info: String, key: String): String? {
    val token = "$key="
    var from = 0
    while (from <= info.length - token.length) {
        val at = info.indexOf(token, from)
        if (at < 0) return null
        if (at == 0 || info[at - 1] == ';') {
            val start = at + token.length
            val semi = info.indexOf(';', start)
            return if (semi < 0) info.substring(start) else info.substring(start, semi)
        }
        from = at + token.length
    }
    return null
}

/** The [index]th `:` field of a sample, or null when the sample has fewer fields. */
private fun colonField(sample: String, index: Int): String? {
    var start = 0
    var seen = 0
    var cursor = 0
    while (cursor <= sample.length) {
        if (cursor == sample.length || sample[cursor] == ':') {
            if (seen == index) return sample.substring(start, cursor)
            seen++
            start = cursor + 1
        }
        if (cursor == sample.length) break
        cursor++
    }
    return null
}

/** Tab boundaries of one VCF line. Fields are sliced only when a caller asks for them. */
private class TabFields(private val line: String) {
    private val starts = ArrayList<Int>()
    private val ends = ArrayList<Int>()

    init {
        var start = 0
        var cursor = 0
        while (cursor <= line.length) {
            if (cursor == line.length || line[cursor] == '\t') {
                starts.add(start)
                ends.add(cursor)
                start = cursor + 1
            }
            if (cursor == line.length) break
            cursor++
        }
    }

    fun text(index: Int): String? = if (index in starts.indices) line.substring(starts[index], ends[index]) else null

    fun samples(): List<String> {
        if (starts.size <= FORMAT) return emptyList()
        val values = ArrayList<String>(starts.size - FORMAT)
        for (index in FORMAT until starts.size) values.add(line.substring(starts[index], ends[index]))
        return values
    }

    private companion object {
        const val FORMAT = 9
    }
}
