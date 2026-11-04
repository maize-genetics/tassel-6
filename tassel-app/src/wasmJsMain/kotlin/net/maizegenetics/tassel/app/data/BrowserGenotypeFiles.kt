package net.maizegenetics.tassel.app.data

import kotlinx.browser.document
import kotlinx.coroutines.suspendCancellableCoroutine
import kotlinx.io.Buffer
import kotlinx.io.Source
import kotlinx.io.writeString
import org.w3c.dom.HTMLInputElement
import org.w3c.files.File
import org.w3c.files.FileReader
import org.w3c.xhr.TEXT
import org.w3c.xhr.XMLHttpRequest
import org.w3c.xhr.XMLHttpRequestResponseType
import kotlin.coroutines.resume
import kotlin.coroutines.resumeWithException
import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsString

/**
 * Browser file loading. There is no file system here, so both paths end in bytes handed to the reader through a
 * `kotlinx.io.Buffer`: the tutorial dataset comes over HTTP from the page the app is served from, and a file the
 * user picks comes from an `<input type="file">` read by a `FileReader`.
 *
 * HapMap is a text format, so both reads ask for text and write it into the buffer as UTF-8. That crosses the
 * JavaScript boundary once per file instead of once per byte, which an `ArrayBuffer` copied element by element would
 * do, and it needs no dependency beyond the browser APIs Compose for Wasm already puts on the classpath.
 */
class BrowserGenotypeFiles : GenotypeFiles {
    override suspend fun tutorial(): GenotypeFile = TextGenotypeFile(TUTORIAL_FILE_NAME) { get(TUTORIAL_FILE_NAME) }

    override suspend fun choose(): GenotypeFile? {
        val file = pickFile() ?: return null
        return TextGenotypeFile(file.name) { readAsText(file) }
    }
}

/** A file whose bytes come from text fetched when it is opened. */
private class TextGenotypeFile(override val name: String, private val text: suspend () -> String) : GenotypeFile {
    override suspend fun open(): Source = Buffer().apply { writeString(text()) }
}

/**
 * Shows a file chooser and returns what the user picked, or null when they cancel.
 *
 * The input element is added to the document because Safari ignores `click()` on a detached one, hidden because the
 * app's own UI is the chooser the user sees, and removed again either way.
 */
private suspend fun pickFile(): File? = suspendCancellableCoroutine { continuation ->
    val input = document.createElement("input") as HTMLInputElement
    input.type = "file"
    input.accept = ".txt,.hmp.txt,text/plain"
    input.style.display = "none"
    document.body?.appendChild(input)

    var answered = false
    fun answer(file: File?) {
        if (answered) return
        answered = true
        input.remove()
        continuation.resume(file)
    }
    input.addEventListener("change") { answer(input.files?.item(0)) }
    // Fired when the chooser is dismissed without a file. Without it the coroutine would never resume.
    input.addEventListener("cancel") { answer(null) }
    continuation.invokeOnCancellation { input.remove() }
    input.click()
}

/** Reads a picked file as text. */
@OptIn(ExperimentalWasmJsInterop::class)
private suspend fun readAsText(file: File): String = suspendCancellableCoroutine { continuation ->
    val reader = FileReader()
    // readAsText leaves the text in result, typed as JsAny? because a FileReader can also produce an ArrayBuffer.
    reader.onload = { continuation.resume((reader.result as JsString).toString()) }
    reader.onerror = { continuation.resumeWithException(IllegalStateException("could not read ${file.name}")) }
    continuation.invokeOnCancellation { reader.abort() }
    reader.readAsText(file)
}

/** Fetches [url] relative to the page as text. */
private suspend fun get(url: String): String = suspendCancellableCoroutine { continuation ->
    val request = XMLHttpRequest()
    request.open("GET", url)
    request.responseType = XMLHttpRequestResponseType.TEXT
    request.onload = {
        if (request.status.toInt() == HTTP_OK) {
            continuation.resume(request.responseText)
        } else {
            continuation.resumeWithException(
                IllegalStateException("GET $url returned ${request.status} ${request.statusText}"),
            )
        }
    }
    request.onerror = { continuation.resumeWithException(IllegalStateException("GET $url failed")) }
    continuation.invokeOnCancellation { request.abort() }
    request.send()
}

private const val HTTP_OK = 200
