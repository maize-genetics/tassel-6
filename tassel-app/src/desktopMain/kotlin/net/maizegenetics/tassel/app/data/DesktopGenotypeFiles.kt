package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.swing.Swing
import kotlinx.coroutines.withContext
import kotlinx.io.Source
import kotlinx.io.asSource
import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.openRandomAccessSource
import java.awt.FileDialog
import java.awt.Frame

/**
 * Desktop file loading: a native [FileDialog] for opening, and the tutorial dataset read off the classpath, where
 * `tassel-app/build.gradle.kts` puts it.
 *
 * Reading a path is the desktop half of [GenotypeFiles]; the browser half cannot do it at all, which is why the
 * interface deals in sources rather than paths.
 *
 * @param owner the window the dialog is modal to, or null for an unowned dialog.
 */
class DesktopGenotypeFiles(private val owner: Frame? = null) : GenotypeFiles {
    override suspend fun tutorial(): GenotypeFile = ClasspathGenotypeFile(TUTORIAL_FILE_NAME)

    override suspend fun choose(): List<GenotypeFile>? = withContext(Dispatchers.Swing) {
        // FileDialog is the native chooser on macOS and Windows, which JFileChooser is not. It blocks until the
        // user answers, so it runs on the AWT event thread rather than the Compose one.
        // Unfiltered. A `*.hmp.txt` filter hid VCF, and FilenameFilter does nothing on macOS. The reader
        // picks HapMap or VCF from the file name.
        val picked = pick("Open genotype file") ?: return@withContext null
        if (!picked.name.endsWith(VCF_GZIP_EXTENSION)) return@withContext listOf(picked)
        // A sibling index opens the region step. Without one, the gzipped VCF is read whole, same as a plain VCF.
        val sibling = siblingIndex(picked.path.toString()) ?: return@withContext listOf(picked)
        listOf(picked.withIndex(PathGenotypeFile(Path(sibling))))
    }

    /** One native dialog. Null when the user cancels. */
    private fun pick(title: String): PathGenotypeFile? {
        val dialog = FileDialog(owner, title, FileDialog.LOAD)
        dialog.isVisible = true
        val directory = dialog.directory
        val file = dialog.file
        return if (directory == null || file == null) null else PathGenotypeFile(Path(directory, file))
    }
}

/**
 * The index beside [vcfPath]: `<vcfPath>.csi`, then `<vcfPath>.tbi`.
 *
 * A `.csi` that exists is used even when a `.tbi` is beside it, and even when the `.csi` cannot be read. Null when
 * neither file is there. This is the desktop half of ADR-0021; the CLI keeps its own copy because index discovery
 * stays out of `tassel-core`.
 */
internal fun siblingIndex(vcfPath: String): String? {
    val csi = vcfPath + CSI_EXTENSION
    if (isRegularFile(csi)) return csi
    val tbi = vcfPath + TBI_EXTENSION
    if (isRegularFile(tbi)) return tbi
    return null
}

/** A file at [path]. Tests use this to read `testdata/` without a dialog. */
internal fun pathGenotypeFile(path: String): GenotypeFile = PathGenotypeFile(Path(path))

private fun isRegularFile(path: String): Boolean {
    val metadata = SystemFileSystem.metadataOrNull(Path(path))
    return metadata != null && metadata.isRegularFile
}

/** A file on disk. [openRandomAccess] is what lets a bgzipped VCF inflate on more than one worker. */
private class PathGenotypeFile(val path: Path, override val index: GenotypeFile? = null) : GenotypeFile {
    override val name: String get() = path.name

    /** The same file with [index] attached. The dialog result is copied this way; the path is unchanged. */
    fun withIndex(index: GenotypeFile): PathGenotypeFile = PathGenotypeFile(path, index)

    override suspend fun open(): Source = SystemFileSystem.source(path).buffered()

    override suspend fun openRandomAccess(): RandomAccessSource = openRandomAccessSource(path.toString())
}

/** A file packaged with the app, read straight off the classpath so it streams rather than loading whole. */
private class ClasspathGenotypeFile(override val name: String) : GenotypeFile {
    override suspend fun open(): Source {
        val stream = DesktopGenotypeFiles::class.java.getResourceAsStream("/$name")
            ?: error("$name is missing from the app resources; check the bundleTutorialData task")
        return stream.asSource().buffered()
    }
}

private const val CSI_EXTENSION: String = ".csi"
private const val TBI_EXTENSION: String = ".tbi"
