package net.maizegenetics.tassel.app.data

import net.maizegenetics.tassel.app.logAppMessage
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeReadOptions
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/** The two phenotype files bundled with the tutorial genotype, in tree order. */
val TUTORIAL_PHENOTYPE_FILE_NAMES: List<String> = listOf("mdp_traits.txt", "mdp_population_structure.txt")

/** The tags a TASSEL phenotype file may start with: `<Phenotype>` and the four legacy ones, case ignored. */
private val PHENOTYPE_TAGS = listOf("<phenotype", "<trait", "<covariate", "<factor", "<header")

/** Bytes read to find the first tag. Enough for any comment block a file starts with. */
private const val SNIFF_BYTES = 4096

/** Suffixes of files that are never phenotypes here, so picking a big `.vcf.gz` costs no read to sniff. */
private val BINARY_SUFFIXES = listOf(".gz", ".bcf", ".csi", ".tbi")

/**
 * Whether the first content line of [head] is a phenotype tag. Blank lines, `#` lines, and a byte-order mark before
 * it are skipped, as [PhenotypeFiles.read] skips them. This is `tassel-cli`'s `sniffInput` on bytes instead of a path.
 */
internal fun startsWithPhenotypeTag(head: ByteArray): Boolean {
    val firstLine = head.decodeToString()
        .removePrefix("\uFEFF")
        .lineSequence()
        .map { it.trim() }
        .firstOrNull { it.isNotEmpty() && !it.startsWith("#") }
        ?: return false
    val lower = firstLine.lowercase()
    return PHENOTYPE_TAGS.any { lower.startsWith(it) }
}

/**
 * Whether this file is a phenotype, decided from its content and never its name (except to skip compressed and
 * index files, which cannot be one). Only the first [SNIFF_BYTES] are read: positionally when the file can seek, so
 * the browser does not copy a large genotype file just to look at its first line.
 */
internal suspend fun GenotypeFile.isPhenotype(): Boolean {
    if (BINARY_SUFFIXES.any { name.endsWith(it, ignoreCase = true) }) return false
    val random = openRandomAccess()
    val head = if (random != null) {
        random.use {
            val bytes = ByteArray(minOf(SNIFF_BYTES.toLong(), it.size).toInt())
            bytes.copyOf(maxOf(it.readAt(0, bytes, 0, bytes.size), 0))
        }
    } else {
        open().use { source ->
            source.request(SNIFF_BYTES.toLong())
            val bytes = ByteArray(SNIFF_BYTES)
            val count = source.peek().use { it.readAtMostTo(bytes, 0, SNIFF_BYTES) }
            bytes.copyOf(maxOf(count, 0))
        }
    }
    return startsWithPhenotypeTag(head)
}

/**
 * Reads this file through [PhenotypeFiles.read], naming the table after the file. The reader's warnings go to the
 * app log, since the UI has nowhere to show them.
 *
 * @see net.maizegenetics.phenotype.PhenotypeBuilder.importPhenotypeFile
 */
suspend fun GenotypeFile.readPhenotype(progress: ProgressReporter = ProgressReporter.None): PhenotypeTable {
    val read = open().use { PhenotypeFiles.read(it, PhenotypeReadOptions(name = name), progress) }
    read.warnings.forEach { logAppMessage("$name: $it") }
    return read.table
}
