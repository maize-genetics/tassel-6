package net.maizegenetics.tassel.app.export

import kotlinx.io.Sink
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/**
 * The file name a phenotype dataset is exported under: the dataset name, with a trailing `.txt` replaced by `.txt`.
 * TASSEL 5 and `tassel-cli` write phenotypes as `.txt`, and a derived dataset such as
 * `mdp_traits.txt + filter-phenotype` has no extension to strip.
 */
fun phenotypeFileName(baseName: String): String =
    (if (baseName.endsWith(TEXT_EXTENSION, ignoreCase = true)) baseName.dropLast(TEXT_EXTENSION.length) else baseName) +
        TEXT_EXTENSION

private const val TEXT_EXTENSION = ".txt"

/** Writes [table] as `<Phenotype>` text, the bytes [PhenotypeTable.writeDelimited] produces. */
fun writePhenotype(table: PhenotypeTable, sink: Sink) {
    sink.appendText { table.writeDelimited(it) }
}
