package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.swing.Swing
import kotlinx.coroutines.withContext
import kotlinx.io.Source
import kotlinx.io.asSource
import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
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
        val dialog = FileDialog(owner, "Open genotype file", FileDialog.LOAD)
        dialog.file = "*.hmp.txt"
        dialog.isVisible = true
        val directory = dialog.directory
        val file = dialog.file
        if (directory == null || file == null) null else PathGenotypeFile(Path(directory, file))
    }
}

/** A file on disk. */
private class PathGenotypeFile(private val path: Path) : GenotypeFile {
    override val name: String get() = path.name

    override suspend fun open(): Source = SystemFileSystem.source(path).buffered()
}

/** A file packaged with the app, read straight off the classpath so it streams rather than loading whole. */
private class ClasspathGenotypeFile(override val name: String) : GenotypeFile {
    override suspend fun open(): Source {
        val stream = DesktopGenotypeFiles::class.java.getResourceAsStream("/$name")
            ?: error("$name is missing from the app resources; check the bundleTutorialData task")
        return stream.asSource().buffered()
    }
}
