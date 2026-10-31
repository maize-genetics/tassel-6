package net.maizegenetics.tassel.app

import androidx.compose.ui.window.Window
import androidx.compose.ui.window.application
import net.maizegenetics.tassel.tasselVersion

fun main() = application {
    Window(onCloseRequest = ::exitApplication, title = "TASSEL ${tasselVersion()}") {
        App()
    }
}
