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
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable

/** TASSEL 5's extension for a HapMap export. */
const val HAPMAP_EXTENSION: String = ".hmp.txt"

/** TASSEL 5's extension for a table, matrix, or phenotype export. */
const val TEXT_EXTENSION: String = ".txt"

/** Plain-text VCF, which `--input` reads. */
const val VCF_EXTENSION: String = ".vcf"

/** bgzip-compressed VCF, which nothing can decompress yet. */
const val VCF_GZIP_EXTENSION: String = ".vcf.gz"

/**
 * Reads a genotype file. `tassel-core` takes a `kotlinx.io.Source` and knows nothing about paths
 * (ADR-0017), so opening the file is the CLI's job.
 *
 * The format comes from the extension, as in TASSEL 5: `$VCF_EXTENSION` is read as VCF and
 * everything else as HapMap.
 *
 * @throws kotlinx.io.IOException if the file cannot be read.
 * @throws IllegalArgumentException if it is not a file, or is in a format nothing can read yet.
 */
suspend fun readGenotypes(path: String, progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
    val file = Path(path)
    val metadata = SystemFileSystem.metadataOrNull(file)
    require(metadata != null) { "no such file: $path" }
    require(metadata.isRegularFile) { "not a file: $path" }
    require(!path.endsWith(VCF_GZIP_EXTENSION)) {
        "bgzip-compressed VCF is not supported yet (tassel-6#14 reads plain $VCF_EXTENSION only); " +
            "decompress it first: $path"
    }
    CliLog.info("reading $path")
    return SystemFileSystem.source(file).buffered().use { source ->
        if (path.endsWith(VCF_EXTENSION)) Vcf.read(source, progress) else HapMap.read(source, progress)
    }
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
