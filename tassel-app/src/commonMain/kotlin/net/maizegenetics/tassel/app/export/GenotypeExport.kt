package net.maizegenetics.tassel.app.export

import kotlinx.io.Sink
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.Bcf
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.Vcf

/**
 * A genotype export. The app asks for one of these and then writes that format. It does not infer the format
 * from a file name.
 *
 * @property label the menu text.
 * @property extension the suffix appended to the dataset name.
 */
enum class GenotypeExportFormat(val label: String, val extension: String) {
    HapMap("HapMap", ".hmp.txt"),
    Vcf("VCF", ".vcf"),
    VcfGz("VCF (bgzip)", ".vcf.gz"),
    Bcf("BCF", ".bcf"),
}

/**
 * The file name a genotype dataset is exported under.
 *
 * A trailing `.hmp.txt`, `.vcf.gz`, `.vcf`, or `.bcf` is removed, ignoring case, and [format]'s extension is
 * appended. `.vcf.gz` is matched before `.vcf`. A name that ends with none of those stays as it is, which is
 * how a filtered dataset named `tutorial.hmp.txt + filter-sites` keeps that stem. The writers take the sink
 * directly, which is why a genotype export is not a result pane.
 */
fun genotypeFileName(baseName: String, format: GenotypeExportFormat = GenotypeExportFormat.HapMap): String {
    val suffix = GENOTYPE_EXTENSIONS.firstOrNull { baseName.endsWith(it, ignoreCase = true) }
    val stem = if (suffix == null) baseName else baseName.dropLast(suffix.length)
    return stem + format.extension
}

/** Longer suffixes first, so `.vcf.gz` is not left as `.vcf` with `.gz` still attached. */
private val GENOTYPE_EXTENSIONS = listOf(".vcf.gz", ".hmp.txt", ".vcf", ".bcf")

/**
 * Writes [table] in [format].
 *
 * HapMap, plain VCF, bgzipped VCF, and BCF all go through the core writers. Bgzipped VCF and BCF use the
 * fixed DEFLATE level, not stored blocks.
 */
suspend fun writeGenotype(table: GenotypeTable, format: GenotypeExportFormat, sink: Sink) {
    when (format) {
        GenotypeExportFormat.HapMap -> HapMap.write(table, sink)
        GenotypeExportFormat.Vcf -> Vcf.write(table, sink)
        GenotypeExportFormat.VcfGz -> Vcf.writeBgzf(table, sink)
        GenotypeExportFormat.Bcf -> Bcf.write(table, sink)
    }
}
