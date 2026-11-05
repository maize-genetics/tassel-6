package net.maizegenetics.tassel.io

import kotlinx.io.Sink
import kotlinx.io.Source
import kotlinx.io.readLine
import kotlinx.io.writeString
import net.maizegenetics.tassel.genotype.Allele
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
 * HapMap genotype files: eleven fixed columns, then one column per taxon, and one line per site.
 *
 * Genotypes are nucleotides, written either as one IUPAC character per taxon (`A`, `R`, `N`) or as two allele
 * characters (`AA`, `AG`, `NN`); [read] tells the two apart the way TASSEL 5 does, from the width of the first site's
 * line. Compression is not handled here: hand [read] an already-decompressed [Source].
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromHapMap
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToHapmap
 */
object HapMap {
    /** The columns before the first taxon: `rs#`, `alleles`, `chrom`, `pos`, and seven more. */
    val FIXED_COLUMNS: List<String> = listOf(
        "rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID",
        "QCcode",
    )

    /**
     * Reads a HapMap file.
     *
     * Lines starting with `##` are skipped, the next line names the taxa, and every line after that is one site.
     * Sites that are not in [PositionList] order are sorted, with their genotypes, as TASSEL 5's
     * `BuilderFromHapMap.sortPositions()` does.
     *
     * @throws IllegalArgumentException if the file is empty, a header is malformed, a position is not an integer, a
     *   genotype has no nucleotide meaning, or a line has the wrong number of genotypes.
     */
    suspend fun read(source: Source, progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
        var line = source.readLine() ?: throw IllegalArgumentException("HapMap input is empty")
        while (line.startsWith("##")) {
            // TASSEL 5 reads taxa annotations from ##SAMPLE lines; Taxon carries only a name, so they are skipped.
            line = source.readLine() ?: throw IllegalArgumentException("HapMap input has no header line")
        }
        val taxa = readTaxa(line)

        val firstSite = source.readLine()
            ?: return InMemoryGenotypeTable.wrap(taxa, PositionList(emptyList()), ByteArray(0))
        val oneLetter = isOneLetter(firstSite, taxa.size)

        val positions = ArrayList<Position>()
        val genotypes = ArrayList<ByteArray>()
        val chromosomes = HashMap<String, Chromosome>()
        var siteLine: String? = firstSite
        while (siteLine != null) {
            readSite(siteLine, taxa.size, oneLetter, chromosomes, positions, genotypes)
            if (positions.size % PROGRESS_SITES == 0) progress.report(positions.size, 0, "Reading HapMap")
            siteLine = source.readLine()
        }

        require(positions.size.toLong() * taxa.size <= Int.MAX_VALUE) {
            "HapMap file is too large to hold in memory: ${positions.size} sites x ${taxa.size} taxa"
        }
        val order = sortedOrder(positions)
        val combined = ByteArray(positions.size * taxa.size)
        order.forEachIndexed { site, from -> genotypes[from].copyInto(combined, site * taxa.size) }
        val positionList = PositionList(order.map { positions[it] })
        progress.report(positions.size, positions.size, "Reading HapMap")
        return InMemoryGenotypeTable.wrap(taxa, positionList, combined)
    }

    /**
     * Writes [table] in the format TASSEL 5 exports: no `##` lines, `\n` line endings, `+` in the strand column, and
     * `NA` in the columns HapMap has no value for. The `alleles` column is recomputed from the genotypes, most
     * frequent allele first, so it need not match the file the table was read from.
     *
     * @param diploid writes two allele characters per genotype (`AG`) instead of one IUPAC character (`R`).
     * @param delimiter must be a tab or a space, as in TASSEL 5.
     */
    suspend fun write(
        table: GenotypeTable,
        sink: Sink,
        diploid: Boolean = false,
        delimiter: Char = '\t',
        progress: ProgressReporter = ProgressReporter.None,
    ) {
        require(delimiter == '\t' || delimiter == ' ') { "HapMap delimiter must be a tab or a space: [$delimiter]" }
        val line = StringBuilder()
        FIXED_COLUMNS.forEach { line.append(it).append(delimiter) }
        table.taxa.forEachIndexed { index, taxon ->
            if (index != 0) line.append(delimiter)
            line.append(taxon.name.trim())
        }
        sink.writeString(line.append('\n').toString())

        val numberOfSites = table.numberOfSites
        for (block in table.siteBlocks()) {
            for (site in block.sites) {
                val position = table.positions[site]
                val siteGenotypes = block.siteGenotypes(site)
                line.clear()
                line.append(position.name).append(delimiter)
                line.append(allelesColumn(block.siteStats(site).alleles)).append(delimiter)
                line.append(position.chromosome.name).append(delimiter)
                line.append(position.position.toString()).append(delimiter)
                // TASSEL 5 writes the strand as + for every site and has no value for the last six columns.
                line.append('+')
                repeat(6) { line.append(delimiter).append("NA") }
                for (genotype in siteGenotypes) {
                    line.append(delimiter)
                    if (diploid) {
                        line.append(Nucleotide.alleleString(Genotype.firstAllele(genotype)))
                        line.append(Nucleotide.alleleString(Genotype.secondAllele(genotype)))
                    } else {
                        line.append(Nucleotide.iupac(genotype))
                    }
                }
                sink.writeString(line.append('\n').toString())
            }
            progress.report(block.sites.last + 1, numberOfSites, "Writing HapMap")
        }
        sink.flush()
    }

    /** `A/C`, or `NA` when the site has no called allele, as TASSEL 5's HapMap export writes it. */
    private fun allelesColumn(alleles: List<Byte>): String =
        if (alleles.isEmpty()) "NA" else alleles.joinToString("/") { Nucleotide.alleleString(it) }

    private fun readTaxa(header: String): TaxaList {
        val columns = splitOnWhitespace(header)
        require(columns.size > FIXED_COLUMNS.size) {
            "HapMap header must have ${FIXED_COLUMNS.size} fixed columns and at least one taxon: [$header]"
        }
        return TaxaList.of(columns.subList(FIXED_COLUMNS.size, columns.size))
    }

    /**
     * True when each genotype is one IUPAC character, false when it is two allele characters, decided as TASSEL 5
     * does: from the average number of characters per taxon on the first site's line.
     */
    private fun isOneLetter(firstSite: String, numberOfTaxa: Int): Boolean {
        val start = whitespaceIndex(firstSite, FIXED_COLUMNS.size)
        require(start >= 0) { "HapMap site line has fewer than ${FIXED_COLUMNS.size} columns: ${excerpt(firstSite)}" }
        // The genotypes are one character plus a delimiter each, except the last, hence the + 1.
        val average = (firstSite.length - start + 1).toDouble() / numberOfTaxa
        return when {
            average > 1.99 && average < 2.01 -> true
            average > 2.99 && average < 3.01 -> false
            else -> throw IllegalArgumentException(
                "HapMap genotypes must be one or two characters each, and the first site must have " +
                    "$numberOfTaxa of them: ${excerpt(firstSite)}",
            )
        }
    }

    private fun readSite(
        line: String,
        numberOfTaxa: Int,
        oneLetter: Boolean,
        chromosomes: HashMap<String, Chromosome>,
        positions: MutableList<Position>,
        genotypes: MutableList<ByteArray>,
    ) {
        val tabs = IntArray(FIXED_COLUMNS.size)
        var found = 0
        var index = 0
        while (found < tabs.size && index < line.length) {
            if (line[index] == '\t') tabs[found++] = index
            index++
        }
        require(found == tabs.size) {
            "HapMap site line must have ${FIXED_COLUMNS.size} tab-separated columns before the genotypes: " +
                excerpt(line)
        }

        val snpId = line.substring(0, tabs[SNP_ID_COLUMN])
        val variants = line.substring(tabs[VARIANT_COLUMN - 1] + 1, tabs[VARIANT_COLUMN])
        val chromosome = line.substring(tabs[CHROMOSOME_COLUMN - 1] + 1, tabs[CHROMOSOME_COLUMN])
        val text = line.substring(tabs[POSITION_COLUMN - 1] + 1, tabs[POSITION_COLUMN])
        val position = text.toIntOrNull()
            ?: throw IllegalArgumentException("HapMap position must be an integer: [${text.trim()}] in $snpId")
        positions.add(
            Position(
                chromosome = chromosomes.getOrPut(chromosome) { Chromosome.of(chromosome) },
                position = position,
                snpId = snpId,
                knownVariants = if (variants.isEmpty()) emptyList() else variants.split('/'),
            ),
        )

        val start = tabs[FIXED_COLUMNS.size - 1] + 1
        val step = if (oneLetter) 2 else 3
        val site = ByteArray(numberOfTaxa)
        var taxon = 0
        var at = start
        while (at < line.length) {
            require(taxon < numberOfTaxa) { "HapMap site $snpId has more than $numberOfTaxa genotypes" }
            val genotype = if (oneLetter) {
                Nucleotide.genotype(line[at])
            } else {
                require(at + 1 < line.length) { "HapMap site $snpId ends inside a genotype" }
                // TASSEL 5 reads the two characters in reverse: the second is the first allele of the genotype.
                Genotype.of(
                    Genotype.secondAllele(Nucleotide.genotype(line[at + 1])),
                    Genotype.secondAllele(Nucleotide.genotype(line[at])),
                )
            }
            require(genotype != UNDEFINED_GENOTYPE) {
                "HapMap site $snpId has an illegal genotype: [${line.substring(at, at + step - 1)}]"
            }
            site[taxon++] = genotype
            at += step
        }
        require(taxon == numberOfTaxa) { "HapMap site $snpId has $taxon genotypes, expected $numberOfTaxa" }
        genotypes.add(site)
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

    /** Index of the character after the [count]th whitespace character, or -1 if the line has fewer than [count]. */
    private fun whitespaceIndex(line: String, count: Int): Int {
        var found = 0
        for (index in line.indices) {
            if (line[index].isFieldSeparator() && ++found == count) return index + 1
        }
        return -1
    }

    private fun splitOnWhitespace(line: String): List<String> {
        val parts = ArrayList<String>()
        var start = 0
        for (index in line.indices) {
            if (line[index].isFieldSeparator()) {
                parts.add(line.substring(start, index))
                start = index + 1
            }
        }
        parts.add(line.substring(start))
        while (parts.isNotEmpty() && parts.last().isEmpty()) parts.removeAt(parts.size - 1)
        return parts
    }

    /** The characters Java's `\s` matches, which is how TASSEL 5 splits HapMap header lines. */
    private fun Char.isFieldSeparator(): Boolean = this == ' ' || this in '\t'..'\r'

    private fun excerpt(line: String): String =
        line.take(EXCERPT_LENGTH) + if (line.length > EXCERPT_LENGTH) "..." else ""

    private const val SNP_ID_COLUMN = 0
    private const val VARIANT_COLUMN = 1
    private const val CHROMOSOME_COLUMN = 2
    private const val POSITION_COLUMN = 3
    private const val EXCERPT_LENGTH = 25

    /**
     * Sites between progress reports. Each report is a suspension point (ADR-0019), so this is also
     * how often a single-threaded caller gets the thread back. The tutorial file has 3,093 sites, so
     * anything much larger reports once, after the parse, which is too late to be worth reporting.
     */
    private const val PROGRESS_SITES = 256
    private val UNDEFINED_GENOTYPE: Byte = Genotype.of(Allele.UNDEFINED, Allele.UNDEFINED)
}
