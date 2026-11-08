package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
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
        val position = text.toIntOrNull()
            ?: throw IllegalArgumentException("VCF position must be an integer: [${text.trim()}] in $site")
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

    /**
     * One taxon's call, from the `GT` field of `window[from until to]`.
     *
     * The separator between the two allele indexes is ignored, so `0|1` and `0/1` give the same genotype. TASSEL 5
     * also reads only the characters at offsets 0 and 2, which is why an allele index above 9 is not supported.
     */
    private fun call(window: ByteArray, from: Int, to: Int, alleles: ByteArray, site: String): Byte {
        require(from < to) { "VCF site $site has an empty genotype column" }
        // TASSEL 5 treats any column starting with a dot as missing, which covers ./. and .|. and a bare dot.
        if (window[from] == DOT) return Genotype.UNKNOWN
        var end = from
        while (end < to && window[end] != COLON) end++
        val first: Int
        val second: Int
        if (end - from == 1) {
            first = digit(window[from])
            second = first
        } else {
            require(end - from >= 3) {
                // TASSEL 5 indexes character 2 unguarded here and dies with a string-index error instead.
                "VCF site $site has a genotype that is not diploid: [${window.decodeToString(from, end)}]"
            }
            first = digit(window[from])
            second = digit(window[from + 2])
        }
        require(first <= alleles.size - 1 && second <= alleles.size - 1) {
            "VCF site $site has a genotype naming an allele it does not list: " +
                "[${window.decodeToString(from, end)}] with ${alleles.size} variants"
        }
        // A half-missing call such as 0/. gives a negative index, which TASSEL 5 checks only after the bound above.
        return if (first < 0 || second < 0) Genotype.UNKNOWN else Genotype.of(alleles[first], alleles[second])
    }

    /**
     * The variants a site's calls index into, from its `REF` and `ALT` columns as TASSEL 5 builds them. Symbolic
     * alleles become one character (`<INS>` to `+`, `<DEL>` to `-`, `*` to `N`), and a missing `ALT` leaves the
     * reference as the only variant.
     */
    private fun variantsOf(reference: String, alternate: String): String {
        if (alternate == ".") return reference
        // TASSEL 5's substitutions, in its order: <INS> and any I become +, <DEL> and any D become -, * becomes N.
        return (reference + "/" + alternate).replace(',', '/')
            .replace("<INS>", "+").replace('I', '+')
            .replace("<DEL>", "-").replace('D', '-')
            .replace("*", "N")
    }

    /** Splits the variants on `/` as Java's `String.split` does, which drops trailing empty parts. */
    private fun splitVariants(variants: String, site: String): List<String> {
        val parts = variants.split('/')
        var size = parts.size
        while (size > 0 && parts[size - 1].isEmpty()) size--
        require(size > 0) { "VCF site $site has no REF or ALT allele: [$variants]" }
        return if (size == parts.size) parts else parts.subList(0, size)
    }

    /**
     * One allele code per variant. A site is an indel when the reference is longer than one character or some
     * alternate is longer than the reference; TASSEL 5 then drops the shared first base and calls what is left of a
     * one-character variant a gap, so `AT`/`A` becomes `T` and `-`.
     */
    private fun allelesOf(variants: List<String>, site: String): ByteArray {
        val indel = variants[0].length > 1 || variants.drop(1).any { it.length > variants[0].length }
        return ByteArray(variants.size) { index ->
            val variant = variants[index]
            val text = when {
                !indel -> variant
                variant.length > 1 -> variant.substring(1)
                else -> "-"
            }
            require(text.isNotEmpty()) { "VCF site $site has an empty allele in [${variants.joinToString("/")}]" }
            Nucleotide.allele(text[0])
        }
    }

    private companion object {
        val LINE_FEED = '\n'.code.toByte()
        val CARRIAGE_RETURN = '\r'.code.toByte()
        val TAB = '\t'.code.toByte()
        val HASH = '#'.code.toByte()
        val DOT = '.'.code.toByte()
        val COLON = ':'.code.toByte()

        fun digit(byte: Byte): Int = (byte.toInt() and 0xFF) - '0'.code

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
