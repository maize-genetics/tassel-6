package net.maizegenetics.tassel.app

import androidx.compose.ui.ExperimentalComposeUiApi
import androidx.compose.ui.window.ComposeViewport
import net.maizegenetics.tassel.app.analysis.LocalAnalysisEngine
import net.maizegenetics.tassel.app.data.BrowserGenotypeFiles
import net.maizegenetics.tassel.app.export.BrowserResultFiles
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.ops.TasselOperations
import kotlin.coroutines.EmptyCoroutineContext

@OptIn(ExperimentalComposeUiApi::class)
fun main() {
    ComposeViewport {
        // The browser has one thread and this app has no worker, so parsing and the analyses share it with the
        // UI; the empty context keeps them on that thread rather than pretending otherwise. Dispatchers.Default
        // would be the same thread with an extra hop. See tassel-6#54 and the PR for tassel-6#22.
        App(
            files = BrowserGenotypeFiles(),
            engine = LocalAnalysisEngine(LocalEngine(TasselOperations.registry, EmptyCoroutineContext)),
            results = BrowserResultFiles(),
        )
    }
}
