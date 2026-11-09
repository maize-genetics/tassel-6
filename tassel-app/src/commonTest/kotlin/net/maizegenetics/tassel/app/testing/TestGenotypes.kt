package net.maizegenetics.tassel.app.testing

import kotlinx.io.Buffer
import kotlinx.io.Source
import kotlinx.io.writeString
import net.maizegenetics.tassel.app.data.GenotypeFile
import net.maizegenetics.tassel.app.data.GenotypeFiles
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/** Fixtures shared by the `tassel-app` tests. No test here reads a file, so none of them needs `testdata/`. */

/**
 * A table of [numberOfTaxa] by [numberOfSites] whose call at a cell is a function of the cell, so a test can tell
 * which cell a genotype came from.
 */
fun testTable(numberOfTaxa: Int, numberOfSites: Int): GenotypeTable = InMemoryGenotypeTable.of(
    taxa = TaxaList.of(List(numberOfTaxa) { "taxon$it" }),
    positions = PositionList(List(numberOfSites) { Position(Chromosome.of("1"), it + 1) }),
) { taxon, site -> cellGenotype(taxon, site) }

/** The call [testTable] puts at a cell: one of the four bases, homozygous, picked so neighbors differ. */
fun cellGenotype(taxon: Int, site: Int): Byte = HOMOZYGOTES[(taxon + 2 * site) % HOMOZYGOTES.size]

private val HOMOZYGOTES = listOf("AA", "CC", "GG", "TT").map { Nucleotide.genotype(it) }

/**
 * A HapMap file with these [taxa] and one site per entry of [calls], each entry holding one IUPAC character per
 * taxon. Site `s` is named `S$s` and sits at position `s + 1` of chromosome 1.
 */
fun hapMapText(taxa: List<String>, calls: List<String>): String {
    calls.forEachIndexed { site, line ->
        require(line.length == taxa.size) { "site $site has ${line.length} calls for ${taxa.size} taxa" }
    }
    val header = (HapMap.FIXED_COLUMNS + taxa).joinToString("\t")
    val sites = calls.mapIndexed { site, line ->
        val fixed = listOf("S$site", "A/C", "1", "${site + 1}", "+", "NA", "NA", "NA", "NA", "NA", "NA")
        (fixed + line.map { it.toString() }).joinToString("\t")
    }
    return (listOf(header) + sites).joinToString("\n") + "\n"
}

/**
 * A [GenotypeFile] over text held in memory. The browser's tutorial request works this way: no path, text handed to
 * the reader through a `kotlinx.io.Buffer`. A picked file is bytes; see [ByteGenotypeFile].
 *
 * @param beforeOpen runs before the bytes are handed over, which lets a test hold a load open or fail it.
 */
class TextGenotypeFile(
    override val name: String,
    private val text: String,
    private val beforeOpen: suspend () -> Unit = {},
) : GenotypeFile {
    /** How many times [open] has been called. */
    var opened: Int = 0
        private set

    override suspend fun open(): Source {
        beforeOpen()
        opened++
        return Buffer().apply { writeString(text) }
    }
}

/**
 * A [GenotypeFile] over bytes held in memory.
 *
 * [TextGenotypeFile] writes UTF-8, which cannot hold a gzip member. Dispatch tests use this for `.vcf.gz` and
 * `.hmp.txt.gz`.
 */
class ByteGenotypeFile(override val name: String, private val bytes: ByteArray) : GenotypeFile {
    override suspend fun open(): Source = Buffer().apply { write(bytes) }
}

/** A [GenotypeFiles] that hands out files a test prepared, and records what the app asked for. */
class FakeGenotypeFiles(
    private val tutorial: GenotypeFile = TextGenotypeFile("tutorial.hmp.txt", SMALL_HAPMAP),
    private val chosen: GenotypeFile? = null,
) : GenotypeFiles {
    /** How many times the app asked for the tutorial dataset. */
    var tutorialRequests: Int = 0
        private set

    /** How many times the app asked the user for a file. */
    var chooseRequests: Int = 0
        private set

    override suspend fun tutorial(): GenotypeFile = tutorial.also { tutorialRequests++ }

    override suspend fun choose(): List<GenotypeFile>? {
        chooseRequests++
        return chosen?.let { listOf(it) }
    }
}

/**
 * A [GenotypeFile] that is only a name, and optionally an [index]. Region-step tests never open it: pairing and the
 * fake indexed opener run first.
 */
class NamedGenotypeFile(override val name: String, override val index: GenotypeFile? = null) : GenotypeFile {
    override suspend fun open(): Source = error("$name was opened")
}

/** Three taxa, two sites: enough to tell that a file was parsed and not invented. */
val SMALL_HAPMAP: String = hapMapText(listOf("B73", "Mo17", "CML247"), listOf("ACM", "GGT"))
