package net.maizegenetics.tassel.app.data

import kotlinx.browser.document
import kotlinx.coroutines.suspendCancellableCoroutine
import kotlinx.io.Buffer
import kotlinx.io.Source
import kotlinx.io.writeString
import net.maizegenetics.tassel.io.RandomAccessSource
import org.w3c.dom.HTMLInputElement
import org.w3c.files.File
import org.w3c.xhr.TEXT
import org.w3c.xhr.XMLHttpRequest
import org.w3c.xhr.XMLHttpRequestResponseType
import kotlin.coroutines.resume
import kotlin.coroutines.resumeWithException

/**
 * Browser file loading. There is no file system here, so both paths end in bytes handed to the reader through a
 * `kotlinx.io.Buffer`: the tutorial dataset comes over HTTP from the page the app is served from, and a file the
 * user picks comes from an `<input type="file">` read by a `FileReader`.
 *
 * A picked file is read as an `ArrayBuffer` when the whole file is wanted. HapMap is text, but gzip is not, and
 * `readAsText` replaces every byte that is not valid UTF-8. A `.vcf.gz` opened for a region is not copied whole:
 * [BlobGenotypeFile] seeks through [BlobRandomAccessSource]. The tutorial stays text: it is the uncompressed HapMap
 * the page serves, and the tutorial button is not a file the user picked.
 *
 * The chooser has `multiple` set. Pairing a VCF with its index is [resolveOpen]'s job, because this page cannot see
 * the files next to the one the user picked.
 */
class BrowserGenotypeFiles : GenotypeFiles {
    override suspend fun tutorial(): GenotypeFile = TextGenotypeFile(TUTORIAL_FILE_NAME) { get(TUTORIAL_FILE_NAME) }

    override suspend fun choose(): List<GenotypeFile>? {
        val files = pickFiles() ?: return null
        return files.map { BlobGenotypeFile(it) }
    }
}

/** A file whose bytes come from text fetched when it is opened. The bundled tutorial uses this. */
private class TextGenotypeFile(override val name: String, private val text: suspend () -> String) : GenotypeFile {
    override suspend fun open(): Source = Buffer().apply { writeString(text()) }
}

/**
 * A file the user picked. [open] copies it whole, which is what HapMap and a plain VCF need. [openRandomAccess]
 * seeks, which is what a region of a bgzipped VCF needs.
 */
private class BlobGenotypeFile(private val file: File) : GenotypeFile {
    override val name: String get() = file.name

    override suspend fun open(): Source = Buffer().apply { write(readBlobBytes(file, file.name)) }

    override suspend fun openRandomAccess(): RandomAccessSource = BlobRandomAccessSource(file)
}

/**
 * Shows a file chooser and returns what the user picked, or null when they cancel.
 *
 * `multiple` is set so a bgzipped VCF and its index can be picked together. The input element is added to the
 * document because Safari ignores `click()` on a detached one, hidden because the app's own UI is the chooser the
 * user sees, and removed again either way. The [File] objects are copied out before the element is removed.
 */
private suspend fun pickFiles(): List<File>? = suspendCancellableCoroutine { continuation ->
    val input = document.createElement("input") as HTMLInputElement
    input.type = "file"
    input.multiple = true
    input.accept = ".vcf,.vcf.gz,.hmp.txt,.hmp.txt.gz,.txt,.csi,.tbi"
    input.style.display = "none"
    document.body?.appendChild(input)

    var answered = false
    fun answer(files: List<File>?) {
        if (answered) return
        answered = true
        input.remove()
        continuation.resume(files)
    }
    input.addEventListener("change") {
        val selected = input.files
        val picked = if (selected == null || selected.length == 0) {
            emptyList()
        } else {
            List(selected.length) { index -> selected.item(index) }.filterNotNull()
        }
        answer(picked.ifEmpty { null })
    }
    // Fired when the chooser is dismissed without a file. Without it the coroutine would never resume.
    input.addEventListener("cancel") { answer(null) }
    continuation.invokeOnCancellation { input.remove() }
    input.click()
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
