package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.position.Chromosome

/**
 * One VCF region, 1-based and inclusive, as `bcftools view -r` reads it (ADR-0021).
 *
 * [start] and [end] are both null for a whole contig (`chr3B`). [end] is null when the query runs
 * through the end of the contig (`chr3B:1000-`). A single position (`chr3B:600000000`) has [start] and
 * [end] equal. Coordinates are [Long] until a kept site becomes a [net.maizegenetics.tassel.position.Position].
 */
data class VcfRegion(
    /** Contig name as the query wrote it, before it is matched to the index. */
    val contig: String,
    val start: Long?,
    val end: Long?,
)

/**
 * Parses one region.
 *
 * `chr3B`, `chr3B:1000-2000`, and `chr3B:1000-` match `bcftools view -r`. `chr3B:600000000` is that
 * single position, not the range through the end of the contig. Commas with a digit on both sides are
 * thousands separators, so `chr3B:599,999,999-600,000,001` is one region. Any other comma is an error.
 * `{name}:1000-2000` is the form for a contig name that contains `:`.
 *
 * This is htslib 1.25 `hts_parse_region` with `HTS_PARSE_ONE_COORD` and thousands separators, and it
 * refuses a list of regions instead of ignoring the comma.
 *
 * @throws IllegalArgumentException when [text] is malformed, a coordinate does not fit in a [Long],
 *   or the range is inverted.
 */
fun parseVcfRegion(text: String): VcfRegion {
    if (text.isEmpty()) malformed(text)
    rejectBadCommas(text)
    val (name, range) = splitName(text)
    if (name.isEmpty()) malformed(text)
    if (range == null) return VcfRegion(name, start = null, end = null)
    val (start, end) = parseRange(text, range)
    return VcfRegion(name, start, end)
}

/**
 * The index contig named by [query].
 *
 * The comparison is exact and case-sensitive first. If that misses, both sides are compared with
 * [Chromosome.of]. One normalized hit is used, and the name returned is the index's own spelling.
 * Zero hits lists [names]. Two or more hits names those index contigs.
 */
internal fun resolveContigName(query: String, names: List<String>): String {
    val exact = names.indexOf(query)
    if (exact >= 0) return names[exact]
    val normalized = normalizedName(query)
    val hits = if (normalized == null) emptyList() else names.filter { normalizedName(it) == normalized }
    if (hits.size == 1) return hits[0]
    if (hits.size > 1) {
        throw IllegalArgumentException(
            "VCF contig $query is ambiguous; it matches ${hits.joinToString(", ")}",
        )
    }
    throw IllegalArgumentException("VCF contig $query is not in the index: ${names.joinToString(", ")}")
}

private fun splitName(text: String): Pair<String, String?> {
    if (text[0] == '{') {
        val close = text.indexOf('}')
        if (close < 1) malformed(text)
        val name = text.substring(1, close)
        if (name.isEmpty()) malformed(text)
        val rest = text.substring(close + 1)
        if (rest.isEmpty()) return name to null
        if (rest.length < 2 || rest[0] != ':') malformed(text)
        return name to rest.substring(1)
    }
    val colon = text.lastIndexOf(':')
    if (colon < 0) return text to null
    if (colon == 0 || colon == text.lastIndex) malformed(text)
    return text.substring(0, colon) to text.substring(colon + 1)
}

private fun parseRange(text: String, range: String): Pair<Long, Long?> {
    val dash = range.indexOf('-')
    if (dash < 0) {
        val position = coordinate(text, range)
        return position to position
    }
    if (dash == 0) malformed(text)
    val start = coordinate(text, range.substring(0, dash))
    val endText = range.substring(dash + 1)
    if (endText.isEmpty()) return start to null
    val end = coordinate(text, endText)
    if (start > end) {
        throw IllegalArgumentException("VCF region $text is inverted: $start is after $end")
    }
    return start to end
}

private fun coordinate(text: String, raw: String): Long {
    if (raw.isEmpty()) malformed(text)
    val digits = StringBuilder(raw.length)
    for (index in raw.indices) {
        val char = raw[index]
        if (char == ',') continue
        if (char !in '0'..'9') malformed(text)
        digits.append(char)
    }
    if (digits.isEmpty()) malformed(text)
    val value = digits.toString().toLongOrNull()
        ?: throw IllegalArgumentException("VCF region $text has a number that does not fit in a Long")
    if (value < 1) malformed(text)
    return value
}

private fun rejectBadCommas(text: String) {
    for (index in text.indices) {
        if (text[index] != ',') continue
        val before = index > 0 && text[index - 1].isDigit()
        val after = index + 1 < text.length && text[index + 1].isDigit()
        if (!before || !after) {
            throw IllegalArgumentException("VCF region $text has a comma that is not a thousands separator")
        }
    }
}

private fun normalizedName(name: String): String? {
    if (name.isEmpty()) return null
    return try {
        Chromosome.of(name).name
    } catch (_: IllegalArgumentException) {
        null
    }
}

private fun malformed(text: String): Nothing = throw IllegalArgumentException("VCF region $text is malformed")
