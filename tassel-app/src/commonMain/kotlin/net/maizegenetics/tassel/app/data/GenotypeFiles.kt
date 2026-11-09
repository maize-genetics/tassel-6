package net.maizegenetics.tassel.app.data

import kotlinx.io.Source
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.io.decompressing
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
     * more than one worker. Bytes already in memory, which is every file the browser opens, return null and stay
     * on [Source.decompressing].
     */
    suspend fun openRandomAccess(): RandomAccessSource? = null
}

/**
 * Where the app gets genotype files. Desktop opens the file system through a native dialog; the browser reads an
 * `<input type="file">` and an HTTP request. Both sit behind this one interface so the shared UI, and the analysis
 * screens of tassel-6#22, never learn which target they are on.
 */
interface GenotypeFiles {
    /** The tutorial dataset shipped with the app ([TUTORIAL_FILE_NAME]). */
    suspend fun tutorial(): GenotypeFile

    /** Asks the user for a file, or returns null when they cancel. */
    suspend fun choose(): GenotypeFile?
}

/**
 * Reads this file as VCF or HapMap, reporting progress.
 *
 * The name chooses the reader, as TASSEL 5's "Make Best Guess" does for these suffixes: [VCF_GZIP_EXTENSION] and
 * [VCF_EXTENSION] are VCF, and everything else, including a gzipped HapMap, is HapMap. `.vcf.gz` does not end with
 * `.vcf`, so that suffix is matched on its own. Gzip is detected from the first bytes, not the name.
 *
 * A `.vcf.gz` that [openRandomAccess] can open and whose first member is BGZF is read through [Vcf.read] of that
 * source, so desktop inflate can use more than one worker. Plain gzip, and every in-memory file, goes through
 * [Source.decompressing]. A missing BGZF EOF block stays a successful read; the app has nowhere to show the warning
 * the CLI prints.
 *
 * Parsing suspends at every progress report (ADR-0019), so a reporter that yields keeps a single-threaded caller
 * responsive. Which context the parse runs on is still the caller's choice; see `AppViewModel`.
 *
 * @see net.maizegenetics.analysis.data.FileLoadPlugin
 */
suspend fun GenotypeFile.readGenotypes(progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
    if (name.endsWith(VCF_GZIP_EXTENSION)) {
        openRandomAccess()?.use { random ->
            if (random.isBgzf()) return Vcf.read(random, progress).table
        }
    }
    return open().use { raw ->
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
private const val VCF_GZIP_EXTENSION: String = ".vcf.gz"
