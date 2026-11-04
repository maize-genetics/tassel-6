package net.maizegenetics.tassel.app

import androidx.compose.ui.Alignment
import androidx.compose.ui.unit.DpSize
import androidx.compose.ui.unit.dp
import androidx.compose.ui.window.Window
import androidx.compose.ui.window.WindowPosition
import androidx.compose.ui.window.application
import androidx.compose.ui.window.rememberWindowState
import kotlinx.coroutines.Dispatchers
import net.maizegenetics.tassel.app.data.DesktopGenotypeFiles
import net.maizegenetics.tassel.tasselVersion

fun main() = application {
    val state = rememberWindowState(
        size = DpSize(1280.dp, 800.dp),
        position = WindowPosition(Alignment.Center),
    )
    Window(onCloseRequest = ::exitApplication, state = state, title = "TASSEL ${tasselVersion()}") {
        // Desktop has threads, so HapMap parsing stays off the Compose thread.
        App(files = DesktopGenotypeFiles(window), readContext = Dispatchers.Default)
    }
}
