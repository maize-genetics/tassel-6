package net.maizegenetics.tassel.app

import androidx.compose.ui.Alignment
import androidx.compose.ui.unit.DpSize
import androidx.compose.ui.unit.dp
import androidx.compose.ui.window.Window
import androidx.compose.ui.window.WindowPosition
import androidx.compose.ui.window.application
import androidx.compose.ui.window.rememberWindowState
import kotlinx.coroutines.Dispatchers
import net.maizegenetics.tassel.app.analysis.LocalAnalysisEngine
import net.maizegenetics.tassel.app.data.DesktopGenotypeFiles
import net.maizegenetics.tassel.app.export.DesktopResultFiles
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.tasselVersion

fun main() = application {
    val state = rememberWindowState(
        size = DpSize(1280.dp, 800.dp),
        position = WindowPosition(Alignment.Center),
    )
    Window(onCloseRequest = ::exitApplication, state = state, title = "TASSEL ${tasselVersion()}") {
        // Desktop has threads, so HapMap parsing and the analyses stay off the Compose thread.
        App(
            files = DesktopGenotypeFiles(window),
            engine = LocalAnalysisEngine(LocalEngine(TasselOperations.registry, Dispatchers.Default)),
            results = DesktopResultFiles(window),
            readContext = Dispatchers.Default,
        )
    }
}
