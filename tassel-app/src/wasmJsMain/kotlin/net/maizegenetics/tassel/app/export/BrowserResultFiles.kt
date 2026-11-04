package net.maizegenetics.tassel.app.export

import kotlinx.browser.document
import kotlinx.io.Buffer
import kotlinx.io.Sink
import kotlinx.io.readString
import org.w3c.dom.HTMLAnchorElement
import org.w3c.dom.url.URL
import org.w3c.files.Blob
import org.w3c.files.BlobPropertyBag
import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsArray
import kotlin.js.set

/**
 * Browser export: the file is built in memory and handed to the page as a download, because there is no file
 * system here and no way to stream into one.
 *
 * The whole file is buffered before the download starts, which is the cost of having no file system; for the
 * results this PoC produces (a 281 x 281 kinship matrix is about 1 MB of text) that is acceptable.
 */
class BrowserResultFiles : ResultFiles {
    @OptIn(ExperimentalWasmJsInterop::class)
    override suspend fun save(fileName: String, write: (Sink) -> Unit): Boolean {
        val buffer = Buffer()
        write(buffer)
        val parts = JsArray<JsAny?>()
        parts[0] = buffer.readString().toJsString()
        val blob = Blob(parts, BlobPropertyBag(type = "text/plain;charset=utf-8"))
        val url = URL.createObjectURL(blob)
        try {
            val anchor = document.createElement("a") as HTMLAnchorElement
            anchor.href = url
            anchor.download = fileName
            // Safari ignores click() on a detached element, the same reason the file input is attached.
            anchor.style.display = "none"
            document.body?.appendChild(anchor)
            anchor.click()
            anchor.remove()
        } finally {
            URL.revokeObjectURL(url)
        }
        return true
    }
}
