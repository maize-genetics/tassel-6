package net.maizegenetics.tassel.app.data

import kotlinx.io.Source
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.Bcf
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.io.decompressing
import net.maizegenetics.tassel.io.isBcf22
import net.maizegenetics.tassel.io.isBgzf
import net.maizegenetics.tassel.ops.ProgressReporter

/** File name of the tutorial dataset bundled with both targets. */
const val TUTORIAL_FILE_NAME: String = "mdp_genotype.hmp.txt"

/**
 * A genotype file the app can read: a display [name] and a way to get at its bytes.
 *
 * [open] suspends and yields bytes rather than a path because the browser has no file system: a file the user picked
 * there is read asynchronously out of the page. `tassel-core` readers take a `kotlinx.io.Source` and know nothing
 * about paths for the same reason (ADR-0017).
 */
interface GenotypeFile {
    /** File name, without a directory. */
    val name: String

    /** Opens the bytes. The caller closes the returned source. */
    suspend fun open(): Source

    /**
     * Opens the same bytes for positional reads, or null when this file is only a stream.
     *
     * A path on disk returns one, so a bgzipped VCF can use [Vcf.read] of a [RandomAccessSource] and inflate on
     * more than one worker. A file the browser user picked returns one over that `File`, so a region read can seek
     * without copying the file into memory. Bytes already copied into memory return null and stay on
     * [Source.decompressing].
     */
    suspend fun openRandomAccess(): RandomAccessSource? = null

    /**
     * The index that belongs with this file, when the platform already resolved one.
     *
     * Desktop sets this from the `.csi` or `.tbi` beside a `.vcf.gz`, or from the `.csi` beside a `.bcf`. The browser
     * leaves it null: the page cannot see neighboring files, so [resolveOpen] pairs the picked names instead. A set
     * index is used even when its name is not `<data>.csi` or `<data>.tbi`, which is the desktop sibling and the
     * CLI's `--index`.
     */
    val index: GenotypeFile?
        get() = null
}

/**
 * Where the app gets genotype files. Desktop opens the file system through a native dialog; the browser reads an
 * `<input type="file">` and an HTTP request. Both sit behind this one interface so the shared UI never learns which
 * target it is on. A bgzipped VCF's or BCF's index is a second [GenotypeFile], not a path: `commonMain` has no file
 * system (ADR-0017, ADR-0022).
 */
interface GenotypeFiles {
    /** The tutorial dataset shipped with the app ([TUTORIAL_FILE_NAME]). */
    suspend fun tutorial(): GenotypeFile

    /**
     * The phenotype files that ship with the tutorial genotype ([TUTORIAL_PHENOTYPE_FILE_NAMES]), opened beside it.
     * Empty by default, for a target or a test with none.
     */
    suspend fun tutorialPhenotypes(): List<GenotypeFile> = emptyList()

    /**
     * Asks the user for files, or returns null when they cancel.
     *
     * One file is the usual open. A `.vcf.gz` or `.bcf` may arrive with [GenotypeFile.index] set, or as one of several
     * files the browser user picked together. [resolveOpen] decides which is the data file and which is the index.
     */
    suspend fun choose(): List<GenotypeFile>?
}

/**
 * Reads this file as BCF, VCF, or HapMap, reporting progress.
 *
 * [BCF_EXTENSION] is always read with [Bcf.read]. A BGZF file with any other name is also BCF when [isBcf22] is
 * true, so a misnamed BCF does not reach [Vcf.read] or [HapMap.read] (ADR-0022). [VCF_GZIP_EXTENSION] and
 * [VCF_EXTENSION] are VCF, and everything else, including a gzipped HapMap, is HapMap. `.vcf.gz` does not end with
 * `.vcf`, so that suffix is matched on its own. Gzip is detected from the first bytes, not the name.
 *
 * A BGZF file that [openRandomAccess] can open is read through [Bcf.read] or [Vcf.read] of that source, so desktop
 * inflate can use more than one worker. Raw BCF, plain gzip, and every in-memory file go through [Source]. A
 * missing BGZF EOF block stays a successful read; the app has nowhere to show the warning the CLI prints.
 *
 * Parsing suspends at every progress report (ADR-0019), so a reporter that yields keeps a single-threaded caller
 * responsive. Which context the parse runs on is still the caller's choice; see `AppViewModel`.
 *
 * @see net.maizegenetics.analysis.data.FileLoadPlugin
 */
suspend fun GenotypeFile.readGenotypes(progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
    val bcfName = name.endsWith(BCF_EXTENSION)
    openRandomAccess()?.use { random ->
        if (random.isBgzf()) {
            // The name forces BCF even when the magic is wrong, so HapMap bytes named `.bcf` stay an error from
            // Bcf.read. isBcf22 covers a BGZF BCF that was not named `.bcf`.
            if (bcfName || random.isBcf22()) return Bcf.read(random, progress).table
            if (name.endsWith(VCF_GZIP_EXTENSION)) return Vcf.read(random, progress).table
        }
    }
    return open().use { raw ->
        if (bcfName) return Bcf.read(raw, progress).table
        val opened = raw.decompressing()
        try {
            if (isVcfName(name)) Vcf.read(opened.source, progress) else HapMap.read(opened.source, progress)
        } finally {
            if (opened.source !== raw) opened.source.close()
        }
    }
}

/** `.vcf.gz` does not end with `.vcf`, so the gzipped suffix has to be tested on its own. */
private fun isVcfName(name: String): Boolean = name.endsWith(VCF_GZIP_EXTENSION) || name.endsWith(VCF_EXTENSION)

private const val VCF_EXTENSION: String = ".vcf"

/** `.vcf.gz` does not end with `.vcf`. Desktop lookup and [resolveOpen] share this suffix. */
internal const val VCF_GZIP_EXTENSION: String = ".vcf.gz"

/** BCF 2.2, bgzipped or raw. Desktop lookup and [resolveOpen] share this suffix. A `.bcf` is never HapMap. */
internal const val BCF_EXTENSION: String = ".bcf"
