package net.maizegenetics.tassel.app.data

import net.maizegenetics.tassel.app.logAppMessage
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.VcfRegion
import net.maizegenetics.tassel.io.isBgzf
import net.maizegenetics.tassel.io.parseVcfRegion
import net.maizegenetics.tassel.io.resolveContigName
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.time.TimeSource

/**
 * What [resolveOpen] decided to do with one chooser result.
 *
 * [WholeFile] is HapMap, gzipped HapMap, or a plain VCF: the reader [GenotypeFile.readGenotypes] already has.
 * [Indexed] is a bgzipped VCF plus the index that belongs with it. [Rejected] is a missing or mismatched index,
 * reported in the UI before any read.
 */
internal sealed interface OpenPlan {
    /** Read [file] from start to end. */
    data class WholeFile(val file: GenotypeFile) : OpenPlan

    /** Open [data] with [index] and ask for a region before reading sites. */
    data class Indexed(val data: GenotypeFile, val index: GenotypeFile) : OpenPlan

    /** [name] is the file the message is about. [message] is shown as a failure and nothing is read. */
    data class Rejected(val name: String, val message: String) : OpenPlan
}

/**
 * Pairs a bgzipped VCF with its index, or keeps a single genotype file on the whole-file path.
 *
 * Names are case-sensitive and the list order does not matter. A [GenotypeFile.index] already set on the `.vcf.gz`
 * wins, even when that index has some other name: desktop's second dialog and the CLI's `--index` both work that
 * way. Otherwise the index must be named `<data>.csi` or `<data>.tbi`. When both are picked, `.csi` wins and the
 * `.tbi` is ignored (ADR-0021). Anything else is a mismatch, and a `.vcf.gz` with no index is not read whole.
 */
internal fun resolveOpen(files: List<GenotypeFile>): OpenPlan {
    if (files.isEmpty()) return OpenPlan.Rejected("", "No file was selected.")
    val dataFiles = files.filter { it.name.endsWith(VCF_GZIP_EXTENSION) }
    if (dataFiles.isEmpty()) {
        if (files.size == 1) return OpenPlan.WholeFile(files[0])
        return OpenPlan.Rejected(files[0].name, "Open one genotype file, or a .vcf.gz with its index.")
    }
    if (dataFiles.size > 1) return OpenPlan.Rejected(dataFiles[0].name, "Open one .vcf.gz with its index.")
    val data = dataFiles[0]
    val attached = data.index
    if (attached != null) {
        val extras = files.filter { it !== data && it !== attached && !it.isCompanionIndex(data.name) }
        if (extras.isNotEmpty()) return mismatch(data, extras[0])
        return OpenPlan.Indexed(data, attached)
    }
    val others = files.filter { it !== data }
    val csiName = data.name + CSI_EXTENSION
    val tbiName = data.name + TBI_EXTENSION
    val unrelated = others.filter { it.name != csiName && it.name != tbiName }
    if (unrelated.isNotEmpty()) return mismatch(data, unrelated[0])
    val csi = others.filter { it.name == csiName }
    val tbi = others.filter { it.name == tbiName }
    if (csi.size > 1 || tbi.size > 1) {
        return OpenPlan.Rejected(data.name, "More than one index was picked for ${data.name}.")
    }
    val index = csi.firstOrNull() ?: tbi.firstOrNull()
        ?: return OpenPlan.Rejected(data.name, missingIndex(data.name))
    return OpenPlan.Indexed(data, index)
}

/**
 * The message to show for [text], or null when it names one region of [contigNames].
 *
 * Blank text has not been parsed yet. Anything else goes through [parseVcfRegion] and [resolveContigName], which are
 * the same checks [IndexedVcf.read] makes, so the prompt can refuse a region before the read starts.
 */
internal fun validateRegion(text: String, contigNames: List<String>): String? {
    if (text.isBlank()) return "Enter a region, or choose a contig."
    val region = try {
        parseVcfRegion(text.trim())
    } catch (failure: IllegalArgumentException) {
        return failure.message ?: "VCF region is malformed"
    }
    return try {
        resolveContigName(region.contig, contigNames)
        null
    } catch (failure: IllegalArgumentException) {
        failure.message
    }
}

/**
 * The CLI's header check (ADR-0021): an index may name a subset of the header's `##contig` ids, and may not name
 * one the header does not. An empty header list, which is the tutorial VCF, skips the check.
 *
 * Returns the failure message, or null when the index is acceptable.
 */
internal fun indexMismatch(
    dataName: String,
    indexName: String,
    indexContigs: List<String>,
    headerContigIds: List<String>,
): String? {
    if (headerContigIds.isEmpty()) return null
    val header = headerContigIds.toSet()
    val missing = indexContigs.filter { it !in header }
    if (missing.isEmpty()) return null
    return "index $indexName does not match $dataName: " +
        "index contigs ${missing.joinToString(", ")} are not in the VCF header " +
        "(${headerContigIds.joinToString(", ")})"
}

/**
 * A bgzipped VCF opened with its index, ready for one region at a time.
 *
 * [openRegionSource] reads the header and the index. [read] reads one region. [close] closes the reader and both
 * sources. The view model holds one of these between the prompt and the load.
 */
internal interface RegionSource : AutoCloseable {
    /** File name of the VCF, without a directory. */
    val fileName: String

    /** Index contigs, in index order, with `##contig` lengths where the header has them. */
    val contigs: List<IndexedVcf.Contig>

    /** The sites in [region]. The caller reports [progress]. */
    suspend fun read(region: VcfRegion, progress: ProgressReporter): GenotypeTable
}

/**
 * Opens [data] with [index] for a region read.
 *
 * Both must be random-access sources. [data] must be BGZF. A bad index, or an index whose contigs are not in the
 * header, closes both sources and throws [IllegalArgumentException]. The byte counter is reset after the header
 * read, so [RegionSource.read] counts only the region's slices.
 */
internal suspend fun openRegionSource(data: GenotypeFile, index: GenotypeFile): RegionSource {
    val raw = data.openRandomAccess()
        ?: throw IllegalArgumentException("${data.name} cannot be read by region")
    val counted = CountingRandomAccessSource(raw)
    val indexSource = try {
        index.openRandomAccess()
    } catch (failure: Throwable) {
        counted.close()
        throw failure
    }
    if (indexSource == null) {
        counted.close()
        throw IllegalArgumentException("${index.name} cannot be read by region")
    }
    var handedOff = false
    var opened: IndexedVcf? = null
    try {
        if (!counted.isBgzf()) throw IllegalArgumentException("${data.name} is not bgzip compressed")
        opened = try {
            IndexedVcf.open(counted, indexSource)
        } catch (failure: IllegalArgumentException) {
            throw IllegalArgumentException(
                "cannot read index ${index.name}: ${failure.message ?: "invalid index"}",
                failure,
            )
        }
        val mismatch = indexMismatch(
            data.name,
            index.name,
            opened.contigs.map { it.name },
            opened.headerContigIds,
        )
        if (mismatch != null) throw IllegalArgumentException(mismatch)
        counted.reset()
        val session = IndexedRegionSource(data.name, opened, counted, indexSource)
        handedOff = true
        return session
    } finally {
        if (!handedOff) {
            opened?.close()
            counted.close()
            indexSource.close()
        }
    }
}

/**
 * Counts bytes copied by [readAt], and the smallest position those bytes came from.
 *
 * [reset] drops the header read, which starts at offset 0, so a later region read can show that it never touched
 * the front of a file past 2 GiB.
 */
internal class CountingRandomAccessSource(private val inner: RandomAccessSource) : RandomAccessSource {
    /** Bytes returned by [readAt] since the last [reset]. */
    var bytes: Long = 0
        private set

    /** Smallest [readAt] position that returned a byte since the last [reset], or null when none has. */
    var minimumOffset: Long? = null
        private set

    /** Forgets [bytes] and [minimumOffset]. Does not move the source. */
    fun reset() {
        bytes = 0
        minimumOffset = null
    }

    override val size: Long get() = inner.size

    override suspend fun readAt(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        val count = inner.readAt(position, destination, destinationOffset, length)
        if (count > 0) {
            bytes += count
            val seen = minimumOffset
            if (seen == null || position < seen) minimumOffset = position
        }
        return count
    }

    override fun close() {
        inner.close()
    }
}

private class IndexedRegionSource(
    override val fileName: String,
    private val indexed: IndexedVcf,
    private val data: CountingRandomAccessSource,
    private val index: RandomAccessSource,
) : RegionSource {
    override val contigs: List<IndexedVcf.Contig> get() = indexed.contigs

    override suspend fun read(region: VcfRegion, progress: ProgressReporter): GenotypeTable {
        val started = TimeSource.Monotonic.markNow()
        val table = indexed.read(region, progress).table
        val elapsed = started.elapsedNow().inWholeMilliseconds
        val offset = data.minimumOffset?.toString() ?: "none"
        logAppMessage(
            "region $fileName ${region.describe()} $elapsed ms " +
                "bytes=${data.bytes} minOffset=$offset fileSize=${data.size}",
        )
        return table
    }

    override fun close() {
        indexed.close()
        data.close()
        index.close()
    }
}

private fun GenotypeFile.isCompanionIndex(dataName: String): Boolean =
    name == dataName + CSI_EXTENSION || name == dataName + TBI_EXTENSION

private fun mismatch(data: GenotypeFile, other: GenotypeFile): OpenPlan.Rejected = OpenPlan.Rejected(
    data.name,
    "${data.name} and ${other.name} do not belong together. " +
        "The index must be named ${data.name}.csi or ${data.name}.tbi.",
)

private fun missingIndex(name: String): String =
    "No index found for $name. Pick the .csi or .tbi file that belongs with it."

/** `chr3B`, `chr3B:1000-2000`, `chr3B:1000-`, or `chr3B:1000`, the forms [parseVcfRegion] accepts. */
private fun VcfRegion.describe(): String = when {
    start == null -> contig
    end == null -> "$contig:$start-"
    start == end -> "$contig:$start"
    else -> "$contig:$start-$end"
}

private const val CSI_EXTENSION: String = ".csi"
private const val TBI_EXTENSION: String = ".tbi"
