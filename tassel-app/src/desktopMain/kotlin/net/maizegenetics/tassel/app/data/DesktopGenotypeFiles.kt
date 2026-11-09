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

    override suspend fun choose(): GenotypeFile? = withContext(Dispatchers.Swing) {
        // FileDialog is the native chooser on macOS and Windows, which JFileChooser is not. It blocks until the
        // user answers, so it runs on the AWT event thread rather than the Compose one.
        // Unfiltered. A `*.hmp.txt` filter hid VCF, and FilenameFilter does nothing on macOS. The reader
        // picks HapMap or VCF from the file name.
        val dialog = FileDialog(owner, "Open genotype file", FileDialog.LOAD)
        dialog.isVisible = true
        val directory = dialog.directory
        val file = dialog.file
        if (directory == null || file == null) null else PathGenotypeFile(Path(directory, file))
    }
}

/** A file at [path]. Tests use this to read `testdata/` without a dialog. */
internal fun pathGenotypeFile(path: String): GenotypeFile = PathGenotypeFile(Path(path))

/** A file on disk. [openRandomAccess] is what lets a bgzipped VCF inflate on more than one worker. */
private class PathGenotypeFile(private val path: Path) : GenotypeFile {
    override val name: String get() = path.name

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
