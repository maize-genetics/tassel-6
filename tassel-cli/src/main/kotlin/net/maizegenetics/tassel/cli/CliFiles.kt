package net.maizegenetics.tassel.cli

import kotlinx.io.Sink
import kotlinx.io.asSink
import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
import kotlinx.io.writeString
import net.maizegenetics.tassel.analysis.GenotypeSummary
import net.maizegenetics.tassel.analysis.PrincipalComponentsResult
import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.io.decompressing
import net.maizegenetics.tassel.io.isBgzf
import net.maizegenetics.tassel.io.openRandomAccessSource
import net.maizegenetics.tassel.io.parseVcfRegion
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable

/** TASSEL 5's extension for a HapMap export. */
const val HAPMAP_EXTENSION: String = ".hmp.txt"

/** TASSEL 5's extension for a table, matrix, or phenotype export. */
const val TEXT_EXTENSION: String = ".txt"

/** Plain-text VCF, which `--input` reads. */
const val VCF_EXTENSION: String = ".vcf"

/** bgzip-compressed VCF. The bytes are detected by content; this suffix only selects the VCF reader. */
const val VCF_GZIP_EXTENSION: String = ".vcf.gz"

/**
 * Reads a genotype file. `tassel-core` takes a `kotlinx.io.Source` and knows nothing about paths
 * (ADR-0017), so opening the file is the CLI's job.
 *
 * The format comes from the extension, as in TASSEL 5: [VCF_EXTENSION] and [VCF_GZIP_EXTENSION] are
 * read as VCF, and everything else, including a gzipped HapMap, as HapMap. `.vcf.gz` does not end
 * with `.vcf`, so the gzipped suffix is matched on its own. A `.vcf.gz` whose first member is BGZF
 * is read from a random-access source, so inflate can use more than one worker. Plain gzip stays on
 * the sequential source. Gzip itself is detected from the first bytes. A missing BGZF EOF block is a
 * warning on standard error.
 *
 * [region] reads one region of a bgzipped VCF through [IndexedVcf]. [index] overrides the sibling
 * lookup, which is `<path>.csi` and then `<path>.tbi` (ADR-0021). A `.csi` that is present is the
 * index even when a `.tbi` is beside it. Without [region], [index] is ignored; the command rejects
 * that combination before calling this.
 *
 * @throws kotlinx.io.IOException if the file cannot be read.
 * @throws IllegalArgumentException if it is not a file, the bytes are not a format the reader accepts,
 *   or [region] is set and the file, the index, or the region cannot be read.
 */
suspend fun readGenotypes(
    path: String,
    progress: ProgressReporter = ProgressReporter.None,
    region: String? = null,
    index: String? = null,
): GenotypeTable {
    if (region != null) return readIndexedRegion(path, region, index, progress)
    val file = Path(path)
    val metadata = SystemFileSystem.metadataOrNull(file)
    require(metadata != null) { "no such file: $path" }
    require(metadata.isRegularFile) { "not a file: $path" }
    CliLog.info("reading $path")
    if (path.endsWith(VCF_GZIP_EXTENSION)) {
        openRandomAccessSource(path).use { random ->
            if (random.isBgzf()) {
                val read = Vcf.read(random, progress)
                read.warnings.forEach(CliLog::warn)
                return read.table
            }
        }
    }
    return SystemFileSystem.source(file).buffered().use { raw ->
        val opened = raw.decompressing()
        try {
            val table = if (isVcfPath(path)) Vcf.read(opened.source, progress) else HapMap.read(opened.source, progress)
            opened.warnings.forEach(CliLog::warn)
            table
        } finally {
            if (opened.source !== raw) opened.source.close()
        }
    }
}

/** `.vcf.gz` does not end with `.vcf`, so the gzipped suffix has to be tested on its own. */
private fun isVcfPath(path: String): Boolean = path.endsWith(VCF_GZIP_EXTENSION) || path.endsWith(VCF_EXTENSION)

/**
 * Reads [region] from the bgzipped VCF at [path].
 *
 * The index is [index] when the user named one, otherwise the `.csi` beside [path], then the `.tbi`.
 * A header that lists `##contig` ids must name every index contig. The tutorial header lists none, so
 * that check does not apply to it. Wheat's `chr6D` is in the header and absent from the CSI, so an
 * index is allowed to be a subset.
 */
private suspend fun readIndexedRegion(
    path: String,
    region: String,
    index: String?,
    progress: ProgressReporter,
): GenotypeTable {
    val parsed = parseVcfRegion(region)
    if (!path.endsWith(VCF_GZIP_EXTENSION)) {
        throw IllegalArgumentException("--region needs a bgzipped VCF (.vcf.gz); $path is not one")
    }
    val file = Path(path)
    val metadata = SystemFileSystem.metadataOrNull(file)
    require(metadata != null) { "no such file: $path" }
    require(metadata.isRegularFile) { "not a file: $path" }
    CliLog.info("reading $path region $region")
    return openRandomAccessSource(path).use { data ->
        if (!data.isBgzf()) {
            throw IllegalArgumentException("--region needs a bgzipped VCF; $path is not bgzip compressed")
        }
        val indexPath = resolveIndex(path, index)
        openRandomAccessSource(indexPath).use { indexSource ->
            val opened = try {
                IndexedVcf.open(data, indexSource)
            } catch (e: IllegalArgumentException) {
                throw IllegalArgumentException("cannot read index $indexPath: ${e.message ?: "invalid index"}", e)
            }
            opened.use { vcf ->
                rejectMismatchedIndex(vcf, path, indexPath)
                val read = vcf.read(parsed, progress)
                read.warnings.forEach(CliLog::warn)
                read.table
            }
        }
    }
}

/**
 * The index for [vcfPath]: [explicit] when the user passed `--index`, otherwise `<vcfPath>.csi`,
 * then `<vcfPath>.tbi`. A `.csi` that exists is used even if it cannot be read.
 */
private fun resolveIndex(vcfPath: String, explicit: String?): String {
    if (explicit != null) return explicit
    val csi = "$vcfPath.csi"
    if (isRegularFile(csi)) return csi
    val tbi = "$vcfPath.tbi"
    if (isRegularFile(tbi)) return tbi
    throw IllegalArgumentException(
        "no index found for $vcfPath. Run `bcftools index -c $vcfPath`. " +
            "`tabix -p vcf` only works for contigs under 2^29.",
    )
}

private fun isRegularFile(path: String): Boolean {
    val metadata = SystemFileSystem.metadataOrNull(Path(path))
    return metadata != null && metadata.isRegularFile
}

private fun rejectMismatchedIndex(vcf: IndexedVcf, dataPath: String, indexPath: String) {
    val headerIds = vcf.headerContigIds
    if (headerIds.isEmpty()) return
    val header = headerIds.toSet()
    val missing = vcf.contigs.map { it.name }.filter { it !in header }
    if (missing.isEmpty()) return
    throw IllegalArgumentException(
        "index $indexPath does not match $dataPath: " +
            "index contigs ${missing.joinToString(", ")} are not in the VCF header " +
            "(${headerIds.joinToString(", ")})",
    )
}

/**
 * One file an operation's result turns into. TASSEL 5 writes a multi-dataset result as several
 * files named from one `-export` list; here they are named from one `--output` base.
 *
 * @property part file name suffix for one of several results, or `null` when there is only one.
 * @property extension TASSEL 5's extension for this kind of file.
 */
internal class OutputPart(val part: String?, val extension: String, val write: suspend (Sink) -> Unit)

/**
 * The files [result] becomes. An operation's output type decides them, so adding an operation whose
 * output is one of these types needs no change here.
 *
 * @throws com.github.ajalt.clikt.core.CliktError if nothing knows how to write this type.
 */
internal fun outputParts(result: Any, progress: ProgressReporter = ProgressReporter.None): List<OutputPart> =
    when (result) {
        is GenotypeTable -> listOf(OutputPart(null, HAPMAP_EXTENSION) { HapMap.write(result, it, progress = progress) })
        is DistanceMatrix -> listOf(OutputPart(null, TEXT_EXTENSION) { result.writeDelimited(SinkAppendable(it)) })
        is ResultTable -> listOf(OutputPart(null, TEXT_EXTENSION) { result.writeDelimited(SinkAppendable(it)) })
        is PhenotypeTable -> listOf(OutputPart(null, TEXT_EXTENSION) { result.writeDelimited(SinkAppendable(it)) })
        is GenotypeSummary -> listOfNotNull(
            result.overview?.let { table -> tablePart("overview", table) },
            result.alleles?.let { table -> tablePart("alleles", table) },
            result.sites?.let { table -> tablePart("sites", table) },
            result.taxa?.let { table -> tablePart("taxa", table) },
        )
        is PrincipalComponentsResult -> listOfNotNull(
            // TASSEL 5 exports the components themselves as "pcs"; see testdata/golden/MANIFEST [pca].
            OutputPart("pcs", TEXT_EXTENSION) { result.principalComponents.writeDelimited(SinkAppendable(it)) },
            result.eigenvalues?.let { table -> tablePart("eigenvalues", table) },
            result.eigenvectors?.let { table -> tablePart("eigenvectors", table) },
        )
        else -> unsupported("no writer for a ${result::class.simpleName} result; add one to outputParts()")
    }

private fun tablePart(part: String, table: ResultTable) =
    OutputPart(part, TEXT_EXTENSION) { table.writeDelimited(SinkAppendable(it)) }

/**
 * Writes [result] to files under the [output] base name, or to standard output when it is `null`.
 *
 * The base name works like TASSEL 5's `-export`: the extension for the kind of file is appended
 * unless the base already ends with it, and a result with several parts puts the part name in
 * between (`--output run` gives `run.overview.txt`, `run.sites.txt`, and so on).
 */
internal suspend fun writeResult(result: Any, output: String?, progress: ProgressReporter = ProgressReporter.None) {
    val parts = outputParts(result, progress)
    if (output == null) {
        val sink = System.out.asSink().buffered()
        parts.forEach { it.write(sink) }
        sink.flush()
        return
    }
    for (part in parts) {
        val path = Path(outputPath(output, part.part, part.extension))
        path.parent?.let { parent ->
            if (!SystemFileSystem.exists(parent)) SystemFileSystem.createDirectories(parent)
        }
        SystemFileSystem.sink(path).buffered().use { part.write(it) }
        CliLog.info("wrote $path")
    }
}

/** The file name for one [part] of a result written under the base name [output]. */
internal fun outputPath(output: String, part: String?, extension: String): String {
    val base = output.removeSuffix(extension)
    return if (part == null) base + extension else "$base.$part$extension"
}

/** Lets the core's `Appendable` writers stream into a `kotlinx.io.Sink` without buffering the whole file. */
private class SinkAppendable(private val sink: Sink) : Appendable {
    override fun append(value: Char): Appendable = apply { sink.writeString(value.toString()) }

    override fun append(value: CharSequence?): Appendable = apply { sink.writeString(value.toString()) }

    override fun append(value: CharSequence?, startIndex: Int, endIndex: Int): Appendable =
        apply { sink.writeString(value.toString(), startIndex, endIndex) }
}
