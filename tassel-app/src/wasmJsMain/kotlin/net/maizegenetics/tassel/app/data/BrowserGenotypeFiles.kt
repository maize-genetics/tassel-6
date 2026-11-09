package net.maizegenetics.tassel.app.data

import kotlinx.browser.document
import kotlinx.coroutines.suspendCancellableCoroutine
import kotlinx.io.Buffer
import kotlinx.io.Source
import kotlinx.io.writeString
import net.maizegenetics.tassel.io.arrayBufferToByteArray
import org.w3c.dom.HTMLInputElement
import org.w3c.files.File
import org.w3c.files.FileReader
import org.w3c.xhr.TEXT
import org.w3c.xhr.XMLHttpRequest
import org.w3c.xhr.XMLHttpRequestResponseType
import kotlin.coroutines.resume
import kotlin.coroutines.resumeWithException
import kotlin.js.ExperimentalWasmJsInterop

/**
 * Browser file loading. There is no file system here, so both paths end in bytes handed to the reader through a
 * `kotlinx.io.Buffer`: the tutorial dataset comes over HTTP from the page the app is served from, and a file the
 * user picks comes from an `<input type="file">` read by a `FileReader`.
 *
 * A picked file is read as an `ArrayBuffer`. HapMap is text, but gzip is not, and `readAsText` replaces every byte
 * that is not valid UTF-8. The tutorial stays text: it is the uncompressed HapMap the page serves, and the tutorial
 * button is not a file the user picked.
 */
class BrowserGenotypeFiles : GenotypeFiles {
    override suspend fun tutorial(): GenotypeFile = TextGenotypeFile(TUTORIAL_FILE_NAME) { get(TUTORIAL_FILE_NAME) }

    override suspend fun choose(): GenotypeFile? {
        val file = pickFile() ?: return null
        return BytesGenotypeFile(file.name) { readAsBytes(file) }
    }
}

/** A file whose bytes come from text fetched when it is opened. The bundled tutorial uses this. */
private class TextGenotypeFile(override val name: String, private val text: suspend () -> String) : GenotypeFile {
    override suspend fun open(): Source = Buffer().apply { writeString(text()) }
}

/** A file whose bytes are copied whole, so a gzip member arrives unchanged. */
private class BytesGenotypeFile(override val name: String, private val bytes: suspend () -> ByteArray) : GenotypeFile {
    override suspend fun open(): Source = Buffer().apply { write(bytes()) }
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
    input.accept = ".vcf,.vcf.gz,.hmp.txt,.hmp.txt.gz,.txt"
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

/** Reads a picked file as bytes. `readAsText` would replace anything that is not UTF-8, including every gzip file. */
@OptIn(ExperimentalWasmJsInterop::class)
private suspend fun readAsBytes(file: File): ByteArray = suspendCancellableCoroutine { continuation ->
    val reader = FileReader()
    // readAsArrayBuffer leaves the buffer in result, typed as JsAny? because a FileReader can also produce text.
    reader.onload = {
        val buffer = reader.result
        if (buffer == null) {
            continuation.resumeWithException(IllegalStateException("could not read ${file.name}"))
        } else {
            continuation.resume(arrayBufferToByteArray(buffer))
        }
    }
    reader.onerror = { continuation.resumeWithException(IllegalStateException("could not read ${file.name}")) }
    continuation.invokeOnCancellation { reader.abort() }
    reader.readAsArrayBuffer(file)
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
