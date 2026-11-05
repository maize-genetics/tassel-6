package net.maizegenetics.tassel.app.export

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.swing.Swing
import kotlinx.coroutines.withContext
import kotlinx.io.Sink
import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
import java.awt.FileDialog
import java.awt.Frame

/**
 * Desktop export: a native save dialog, then the file on disk.
 *
 * Writing happens on a background dispatcher, because a kinship matrix is megabytes and the Compose thread should
 * not be the one spelling them out.
 *
 * @param owner the window the dialog is modal to, or null for an unowned dialog.
 */
class DesktopResultFiles(private val owner: Frame? = null) : ResultFiles {
    override suspend fun save(fileName: String, write: suspend (Sink) -> Unit): Boolean {
        val path = withContext(Dispatchers.Swing) {
            // FileDialog is the native chooser on macOS and Windows, as in DesktopGenotypeFiles.
            val dialog = FileDialog(owner, "Save result", FileDialog.SAVE)
            dialog.file = fileName
            dialog.isVisible = true
            val directory = dialog.directory
            val chosen = dialog.file
            if (directory == null || chosen == null) null else Path(directory, chosen)
        } ?: return false
        withContext(Dispatchers.IO) {
            SystemFileSystem.sink(path).buffered().use { write(it) }
        }
        return true
    }
}
