package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position

/** What [VcfRecordParser.tryParse] found in the caller's window. */
internal sealed class VcfParse {
    /** One site. [consumed] runs through the line break, and the caller's genotype array has been filled. */
    class Site(val consumed: Int, val position: Position) : VcfParse()

    /** A `#` line in the body. [consumed] runs through the line break. */
    class Skip(val consumed: Int) : VcfParse()

    /** [tryParse] saw no line break and the caller has not reached EOF. Nothing was consumed. */
    object Incomplete : VcfParse()

    /** No bytes left at EOF. */
    object End : VcfParse()
}

/**
 * Reads one VCF record from a window of bytes the caller owns.
 *
 * This parser does not touch a `Source`. A window with no line break returns [VcfParse.Incomplete], and the
 * caller keeps that tail, appends the next fill, and calls again. A window that begins on any record is a valid
 * start, so a later BGZF reader can seek to a record and stop when its own `tell()` passes a chunk end without
 * the parser holding bytes the caller cannot see.
 *
 * [calls] is filled in place and must have one byte per taxon. The tab indexes are reused across records.
 * `CHROM`, `POS`, `ID`, `REF`, and `ALT` are the only fields decoded to text. `FORMAT` and each `GT` stay bytes.
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
internal class VcfRecordParser(private val header: VcfHeader) {
    private val chromosomes = HashMap<String, Chromosome>()
    private val tabs = IntArray(header.fixed)

    /**
     * Parses the line that starts at [from] in [window], looking no further than [until].
     *
     * [eof] makes a trailing line with no break a record, matching `readLine`. `\n` and `\r\n` end a line. A bare
     * `\r` does not.
     */
    fun tryParse(window: ByteArray, from: Int, until: Int, eof: Boolean, calls: ByteArray): VcfParse {
        if (from >= until) return if (eof) VcfParse.End else VcfParse.Incomplete
        var newline = from
        while (newline < until && window[newline] != LINE_FEED) newline++
        if (newline >= until && !eof) return VcfParse.Incomplete
        val hasBreak = newline < until
        val contentEnd = if (hasBreak && newline > from && window[newline - 1] == CARRIAGE_RETURN) {
            newline - 1
        } else if (hasBreak) {
            newline
        } else {
            until
        }
        val consumed = if (hasBreak) newline + 1 - from else until - from
        if (contentEnd > from && window[from] == HASH) return VcfParse.Skip(consumed)
        return VcfParse.Site(consumed, readSite(window, from, contentEnd, calls))
    }

    private fun readSite(window: ByteArray, from: Int, to: Int, calls: ByteArray): Position {
        val length = to - from
        var found = 0
        var index = 0
        while (found < tabs.size && index < length) {
            if (window[from + index] == TAB) tabs[found++] = index
            index++
        }
        require(found == tabs.size) {
            "VCF site line must have ${header.fixed} tab-separated columns before the genotypes: " +
                vcfExcerpt(window, from, to)
        }

        val chromosome = field(window, from, header.chromosomeColumn)
        val snpId = if (header.snpIdColumn <= 0) {
            null
        } else {
            field(window, from, header.snpIdColumn).takeUnless { it == "." }
        }
        val text = field(window, from, header.positionColumn)
        val site = snpId ?: "$chromosome:$text"
        val asLong = text.toLongOrNull()
            ?: throw IllegalArgumentException("VCF position must be an integer: [${text.trim()}] in $site")
        // The index and BGZF layers keep Long. Position is the GenotypeTable boundary (ADR-0021).
        if (asLong !in Int.MIN_VALUE.toLong()..Int.MAX_VALUE.toLong()) {
            throw IllegalArgumentException(
                "VCF position $asLong is past TASSEL's position limit of ${Int.MAX_VALUE} in $site",
            )
        }
        val position = asLong.toInt()
        val variants = splitVariants(
            variantsOf(field(window, from, header.referenceColumn), field(window, from, header.alternateColumn)),
            site,
        )
        val parsed = Position(
            chromosome = chromosomes.getOrPut(chromosome) { Chromosome.of(chromosome) },
            position = position,
            snpId = snpId,
            // TASSEL 5 also stores the first two alleles as the position's reference and alternate alleles, and
            // the INFO column as position annotations. Position has neither (#10), so a VCF-backed site summary
            // prints N for Ref and Alt (#43); see docs/porting-ledger.md. knownVariants keeps the variant text,
            // as HapMap.read does, which is also all TASSEL 5's knownVariants() stores.
            knownVariants = variants,
        )
        val alleles = allelesOf(variants, site)
        if (header.formatColumn >= 0) requireGt(window, from, site)
        readCalls(window, from + tabs[header.fixed - 1] + 1, to, alleles, site, calls)
        return parsed
    }

    private fun field(window: ByteArray, from: Int, column: Int): String {
        val start = if (column == 0) from else from + tabs[column - 1] + 1
        return window.decodeToString(start, from + tabs[column])
    }

    private fun requireGt(window: ByteArray, from: Int, site: String) {
        val start = if (header.formatColumn == 0) from else from + tabs[header.formatColumn - 1] + 1
        val end = from + tabs[header.formatColumn]
        if (startsWithGt(window, start, end)) return
        val format = window.decodeToString(start, end)
        val message = if (containsGt(window, start, end)) {
            "VCF site $site does not have GT as its first FORMAT field: [$format]"
        } else {
            "VCF site $site has no GT field in FORMAT: [$format]"
        }
        throw IllegalArgumentException(message)
    }

    /** The genotypes of one site, from the tab-separated taxon columns in `window[start until to]`. */
    private fun readCalls(window: ByteArray, start: Int, to: Int, alleles: ByteArray, site: String, calls: ByteArray) {
        val numberOfTaxa = calls.size
        var taxon = 0
        var at = start
        while (true) {
            var end = at
            while (end < to && window[end] != TAB) end++
            require(taxon < numberOfTaxa) { "VCF site $site has more than $numberOfTaxa genotypes" }
            calls[taxon++] = call(window, at, end, alleles, site)
            if (end >= to) break
            at = end + 1
        }
        require(taxon == numberOfTaxa) { "VCF site $site has $taxon genotypes, expected $numberOfTaxa" }
    }

    private companion object {
        val LINE_FEED = '\n'.code.toByte()
        val CARRIAGE_RETURN = '\r'.code.toByte()
        val TAB = '\t'.code.toByte()
        val HASH = '#'.code.toByte()

        fun startsWithGt(window: ByteArray, from: Int, until: Int): Boolean =
            until - from >= 2 && window[from] == 'G'.code.toByte() && window[from + 1] == 'T'.code.toByte()

        fun containsGt(window: ByteArray, from: Int, until: Int): Boolean {
            var index = from
            while (index < until - 1) {
                if (window[index] == 'G'.code.toByte() && window[index + 1] == 'T'.code.toByte()) return true
                index++
            }
            return false
        }
    }
}

internal const val VCF_EXCERPT_LENGTH = 25

internal fun vcfExcerpt(text: String): String =
    text.take(VCF_EXCERPT_LENGTH) + if (text.length > VCF_EXCERPT_LENGTH) "..." else ""

/** First [VCF_EXCERPT_LENGTH] bytes, which is the first 25 characters for the ASCII lines the reader reports. */
internal fun vcfExcerpt(bytes: ByteArray, from: Int, until: Int): String {
    val length = until - from
    val take = minOf(VCF_EXCERPT_LENGTH, length)
    return bytes.decodeToString(from, from + take) + if (length > VCF_EXCERPT_LENGTH) "..." else ""
}
