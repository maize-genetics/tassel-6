package net.maizegenetics.tassel.io

import kotlinx.io.Source
import kotlinx.io.readLine
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * VCF genotype files, genotype (`GT`) field only.
 *
 * [read] takes the variants of a site from its `REF` and `ALT` columns and the calls from the first `FORMAT` field,
 * which TASSEL 5 also requires to be `GT`. Both phased (`0|1`) and unphased (`0/1`) diploid calls are read, and the
 * two are indistinguishable afterwards, because the genotype byte has no phase bit. Every other `FORMAT` field,
 * including `AD`, `DP`, `GQ`, and `PL`, is skipped, so read depth does not survive (tassel-6#14 scope).
 *
 * Compression is not handled here: hand [read] an already-decompressed [Source].
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
object Vcf {
    /** The column names TASSEL 5 looks for in the `#CHROM` header line, before the first taxon. */
    val FIXED_COLUMNS: List<String> = listOf(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    )

    /**
     * Reads a VCF file.
     *
     * Lines starting with `##` are skipped, including the `##SAMPLE=` taxa annotations TASSEL 5 turns into `Taxon`
     * annotations: `Taxon` carries only a name, so they are read and dropped (the maintainer's decision on
     * tassel-6#14, recorded in `docs/porting-ledger.md`). The next line is the `#CHROM` header, which names the taxa,
     * and every line after it that does not start with `#` is one site.
     *
     * Sites that are not in [PositionList] order are sorted with their genotypes, as [HapMap.read] does. TASSEL 5
     * instead refuses an unsorted VCF unless the caller asks for `buildAndSortInMemory()`.
     *
     * @throws IllegalArgumentException if the file is empty, the header is missing a column the reader needs, a
     *   position is not an integer, `GT` is not the first `FORMAT` field, a call names an allele the site does not
     *   list, or a line has the wrong number of genotypes.
     */
    fun read(source: Source, progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
        var line = source.readLine() ?: throw IllegalArgumentException("VCF input is empty")
        while (line.startsWith("##")) {
            line = source.readLine() ?: throw IllegalArgumentException("VCF input has no #CHROM header line")
        }
        val columns = Columns(line)
        val taxa = TaxaList.of(columns.taxaNames)

        val positions = ArrayList<Position>()
        val genotypes = ArrayList<ByteArray>()
        val chromosomes = HashMap<String, Chromosome>()
        var siteLine = source.readLine()
        while (siteLine != null) {
            // TASSEL 5 skips every # line in the body, not only the ones ahead of the header.
            if (!siteLine.startsWith("#")) {
                readSite(siteLine, taxa.size, columns, chromosomes, positions, genotypes)
                if (positions.size % PROGRESS_SITES == 0) progress.report(positions.size, 0, "Reading VCF")
            }
            siteLine = source.readLine()
        }

        require(positions.size.toLong() * taxa.size <= Int.MAX_VALUE) {
            "VCF file is too large to hold in memory: ${positions.size} sites x ${taxa.size} taxa"
        }
        val order = sortedOrder(positions)
        val combined = ByteArray(positions.size * taxa.size)
        order.forEachIndexed { site, from -> genotypes[from].copyInto(combined, site * taxa.size) }
        val positionList = PositionList(order.map { positions[it] })
        progress.report(positions.size, positions.size, "Reading VCF")
        return InMemoryGenotypeTable.wrap(taxa, positionList, combined)
    }

    /**
     * Where each column the reader needs sits in the `#CHROM` header line, as TASSEL 5's `HeaderPositions` finds
     * them: by name, not by order. The taxa start after whichever of `INFO` and `FORMAT` comes last.
     */
    private class Columns(header: String) {
        private val chromosomeColumn: Int
        private val positionColumn: Int
        private val snpIdColumn: Int
        private val referenceColumn: Int
        private val alternateColumn: Int
        private val formatColumn: Int

        /** Number of columns before the first taxon. */
        val fixed: Int

        /** Taxon names, in file order. */
        val taxaNames: List<String>

        init {
            val names = header.split('\t')
            // TASSEL 5 accepts #CHR as well as #CHROM.
            chromosomeColumn = maxOf(names.indexOf("#CHROM"), names.indexOf("#CHR"))
            positionColumn = names.indexOf("POS")
            snpIdColumn = names.indexOf("ID")
            referenceColumn = names.indexOf("REF")
            alternateColumn = names.indexOf("ALT")
            formatColumn = names.indexOf("FORMAT")
            fixed = maxOf(names.indexOf("INFO"), formatColumn) + 1
            require(chromosomeColumn >= 0 && positionColumn > 0 && referenceColumn > 0 && alternateColumn > 0) {
                "VCF header must name #CHROM, POS, REF, and ALT: ${excerpt(header)}"
            }
            require(fixed > maxOf(chromosomeColumn, positionColumn, referenceColumn, alternateColumn)) {
                "VCF header must name INFO or FORMAT after the other columns: ${excerpt(header)}"
            }
            require(names.size > fixed) { "VCF header must name at least one taxon: ${excerpt(header)}" }
            taxaNames = names.subList(fixed, names.size)
        }

        fun chromosome(line: String, tabs: IntArray): String = column(line, tabs, chromosomeColumn)

        fun position(line: String, tabs: IntArray): String = column(line, tabs, positionColumn)

        fun reference(line: String, tabs: IntArray): String = column(line, tabs, referenceColumn)

        fun alternate(line: String, tabs: IntArray): String = column(line, tabs, alternateColumn)

        /** The `ID` column, or null when the file has none. TASSEL 5 also treats `.` as no name. */
        fun snpId(line: String, tabs: IntArray): String? =
            if (snpIdColumn <= 0) null else column(line, tabs, snpIdColumn).takeUnless { it == "." }

        /** The `FORMAT` column, or null when the file has none. */
        fun format(line: String, tabs: IntArray): String? =
            if (formatColumn < 0) null else column(line, tabs, formatColumn)

        private fun column(line: String, tabs: IntArray, index: Int): String =
            line.substring(if (index == 0) 0 else tabs[index - 1] + 1, tabs[index])
    }

    private fun readSite(
        line: String,
        numberOfTaxa: Int,
        columns: Columns,
        chromosomes: HashMap<String, Chromosome>,
        positions: MutableList<Position>,
        genotypes: MutableList<ByteArray>,
    ) {
        val tabs = IntArray(columns.fixed)
        var found = 0
        var index = 0
        while (found < tabs.size && index < line.length) {
            if (line[index] == '\t') tabs[found++] = index
            index++
        }
        require(found == tabs.size) {
            "VCF site line must have ${columns.fixed} tab-separated columns before the genotypes: ${excerpt(line)}"
        }

        val chromosome = columns.chromosome(line, tabs)
        val snpId = columns.snpId(line, tabs)
        val text = columns.position(line, tabs)
        val site = snpId ?: "$chromosome:$text"
        val position = text.toIntOrNull()
            ?: throw IllegalArgumentException("VCF position must be an integer: [${text.trim()}] in $site")
        val variants = splitVariants(variantsOf(columns.reference(line, tabs), columns.alternate(line, tabs)), site)
        positions.add(
            Position(
                chromosome = chromosomes.getOrPut(chromosome) { Chromosome.of(chromosome) },
                position = position,
                snpId = snpId,
                // TASSEL 5 also stores the first two alleles as the position's reference and alternate alleles, and
                // the INFO column as position annotations. Position has neither (#10), so a VCF-backed site summary
                // prints N for Ref and Alt (#43); see docs/porting-ledger.md. knownVariants keeps the variant text,
                // as HapMap.read does, which is also all TASSEL 5's knownVariants() stores.
                knownVariants = variants,
            ),
        )
        val alleles = allelesOf(variants, site)

        columns.format(line, tabs)?.let { format ->
            require(format.startsWith("GT")) {
                // TASSEL 5 tells the two cases apart to say whether GT is misplaced or absent.
                if (format.contains("GT")) {
                    "VCF site $site does not have GT as its first FORMAT field: [$format]"
                } else {
                    "VCF site $site has no GT field in FORMAT: [$format]"
                }
            }
        }
        genotypes.add(readCalls(line, tabs[columns.fixed - 1] + 1, numberOfTaxa, alleles, site))
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

    /** The genotypes of one site, from the tab-separated taxon columns starting at [start]. */
    private fun readCalls(line: String, start: Int, numberOfTaxa: Int, alleles: ByteArray, site: String): ByteArray {
        val calls = ByteArray(numberOfTaxa)
        var taxon = 0
        var at = start
        while (true) {
            var end = at
            while (end < line.length && line[end] != '\t') end++
            require(taxon < numberOfTaxa) { "VCF site $site has more than $numberOfTaxa genotypes" }
            calls[taxon++] = call(line, at, end, alleles, site)
            if (end >= line.length) break
            at = end + 1
        }
        require(taxon == numberOfTaxa) { "VCF site $site has $taxon genotypes, expected $numberOfTaxa" }
        return calls
    }

    /**
     * One taxon's call, from the `GT` field of the column `line[from until to]`.
     *
     * The separator between the two allele indices is ignored, so `0|1` and `0/1` give the same genotype. TASSEL 5
     * also reads only the characters at offsets 0 and 2, which is why an allele index above 9 is not supported.
     */
    private fun call(line: String, from: Int, to: Int, alleles: ByteArray, site: String): Byte {
        require(from < to) { "VCF site $site has an empty genotype column" }
        // TASSEL 5 treats any column starting with a dot as missing, which covers ./. and .|. and a bare dot.
        if (line[from] == '.') return Genotype.UNKNOWN
        var end = from
        while (end < to && line[end] != ':') end++
        val first: Int
        val second: Int
        if (end - from == 1) {
            first = line[from] - '0'
            second = first
        } else {
            require(end - from >= 3) {
                // TASSEL 5 indexes character 2 unguarded here and dies with a string-index error instead.
                "VCF site $site has a genotype that is not diploid: [${line.substring(from, end)}]"
            }
            first = line[from] - '0'
            second = line[from + 2] - '0'
        }
        require(first <= alleles.size - 1 && second <= alleles.size - 1) {
            "VCF site $site has a genotype naming an allele it does not list: " +
                "[${line.substring(from, end)}] with ${alleles.size} variants"
        }
        // A half-missing call such as 0/. gives a negative index, which TASSEL 5 checks only after the bound above.
        return if (first < 0 || second < 0) Genotype.UNKNOWN else Genotype.of(alleles[first], alleles[second])
    }

    /**
     * Site indices in [PositionList] order. Sites already in order keep their file order, so a sorted file is read
     * unchanged and sites sharing a position keep the order they were written in.
     */
    private fun sortedOrder(positions: List<Position>): List<Int> {
        val order = positions.indices.toList()
        val sorted = (1 until positions.size).all { positions[it - 1] <= positions[it] }
        return if (sorted) order else order.sortedBy { positions[it] }
    }

    private fun excerpt(line: String): String =
        line.take(EXCERPT_LENGTH) + if (line.length > EXCERPT_LENGTH) "..." else ""

    private const val EXCERPT_LENGTH = 25
    private const val PROGRESS_SITES = 4096
}
