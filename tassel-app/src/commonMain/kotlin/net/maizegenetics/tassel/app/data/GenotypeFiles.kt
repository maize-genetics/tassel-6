package net.maizegenetics.tassel.app.data

import kotlinx.io.Source
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
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
 * Reads this file as HapMap, reporting progress.
 *
 * Parsing suspends at every progress report (ADR-0019), so a reporter that yields keeps a single-threaded caller
 * responsive. Which context the parse runs on is still the caller's choice; see `AppViewModel`.
 */
suspend fun GenotypeFile.readGenotypes(progress: ProgressReporter = ProgressReporter.None): GenotypeTable =
    open().use { HapMap.read(it, progress) }
