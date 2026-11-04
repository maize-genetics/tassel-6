package net.maizegenetics.tassel.app

import androidx.compose.ui.ExperimentalComposeUiApi
import androidx.compose.ui.window.ComposeViewport
import net.maizegenetics.tassel.app.data.BrowserGenotypeFiles

@OptIn(ExperimentalComposeUiApi::class)
fun main() {
    ComposeViewport {
        // The browser has one thread and this app has no worker, so parsing shares it with the UI; the default
        // context keeps the parse on that thread rather than pretending otherwise.
        App(files = BrowserGenotypeFiles())
    }
}
