package net.maizegenetics.tassel.spike.gram

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import kotlinx.serialization.json.Json
import kotlinx.serialization.json.JsonElement
import kotlinx.serialization.json.JsonObject
import kotlinx.serialization.json.JsonPrimitive
import kotlinx.serialization.json.booleanOrNull
import kotlinx.serialization.json.contentOrNull
import kotlinx.serialization.json.jsonObject
import kotlinx.serialization.json.jsonPrimitive
import java.net.InetSocketAddress
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.Executors

/**
 * Serves the wasm browser distribution and `gram_f32.wgsl`, and writes each POSTed JSON report
 * under the reports directory. JDK `HttpServer` only: no new dependency.
 */
fun main(args: Array<String>) {
    val distribution = Path.of(args[0]).toAbsolutePath().normalize()
    val shader = Path.of(args[1]).toAbsolutePath().normalize()
    val reports = Path.of(args[2]).toAbsolutePath().normalize()
    val port = args.getOrNull(3)?.toInt() ?: 8765
    Files.createDirectories(reports)
    println("host: $machineHost")
    val server = HttpServer.create(InetSocketAddress("127.0.0.1", port), 0)
    server.createContext("/") { exchange -> handle(exchange, distribution, shader, reports) }
    server.executor = Executors.newFixedThreadPool(4)
    server.start()
    println("GPU Gram spike at http://127.0.0.1:$port/")
    println("reports: $reports")
    Thread.currentThread().join()
}

private fun handle(exchange: HttpExchange, distribution: Path, shader: Path, reports: Path) {
    try {
        val path = exchange.requestURI.path
        if (exchange.requestMethod == "GET" && path == "/host") {
            val body = machineHost.toByteArray()
            exchange.responseHeaders.add("Content-Type", "text/plain; charset=utf-8")
            exchange.responseHeaders.add("Cache-Control", "no-store")
            exchange.sendResponseHeaders(200, body.size.toLong())
            exchange.responseBody.use { it.write(body) }
            return
        }
        if (exchange.requestMethod == "POST" && path == "/report") {
            val saved = saveReport(reports, exchange.requestBody.readBytes())
            val ok = """{"ok":true,"file":"${saved.primary}"}""".toByteArray()
            exchange.responseHeaders.add("Content-Type", "application/json")
            exchange.responseHeaders.add("Cache-Control", "no-store")
            exchange.sendResponseHeaders(200, ok.size.toLong())
            exchange.responseBody.use { it.write(ok) }
            println("wrote ${saved.logged}")
            return
        }
        if (exchange.requestMethod != "GET" && exchange.requestMethod != "HEAD") {
            exchange.sendResponseHeaders(405, -1)
            exchange.close()
            return
        }
        val file = if (path == "/gram_f32.wgsl") {
            shader
        } else {
            val relative = if (path == "/" || path.isEmpty()) "index.html" else path.removePrefix("/")
            val resolved = distribution.resolve(relative).normalize()
            if (!resolved.startsWith(distribution) || !Files.isRegularFile(resolved)) {
                exchange.sendResponseHeaders(404, -1)
                exchange.close()
                return
            }
            resolved
        }
        val bytes = Files.readAllBytes(file)
        exchange.responseHeaders.add("Content-Type", contentType(file))
        exchange.responseHeaders.add("Cache-Control", "no-store")
        exchange.sendResponseHeaders(200, bytes.size.toLong())
        if (exchange.requestMethod != "HEAD") {
            exchange.responseBody.use { it.write(bytes) }
        } else {
            exchange.close()
        }
    } catch (error: Exception) {
        val message = (error.message ?: "error").toByteArray()
        exchange.sendResponseHeaders(500, message.size.toLong())
        exchange.responseBody.use { it.write(message) }
    }
}

private val reportJson: Json = Json { prettyPrint = true }

private val machineHost: String = readMachineHost()

private fun readMachineHost(): String {
    val model = try {
        val process = ProcessBuilder("sysctl", "-n", "machdep.cpu.brand_string")
            .redirectErrorStream(true)
            .start()
        val text = process.inputStream.readAllBytes().decodeToString().trim()
        if (process.waitFor() == 0 && text.isNotEmpty()) text else System.getProperty("os.arch")
    } catch (error: Exception) {
        System.err.println("host lookup failed: ${error.message}")
        System.getProperty("os.arch")
    }
    val os = listOf("os.name", "os.version", "os.arch").joinToString(" ") { System.getProperty(it) }
    return "$model; $os"
}

private class SavedReport(val primary: String, val logged: String)

private fun saveReport(reports: Path, body: ByteArray): SavedReport {
    val text = ensureHost(body.toString(Charsets.UTF_8))
    val bytes = text.toByteArray()
    val name = "browser-${System.currentTimeMillis()}.json"
    Files.write(reports.resolve(name), bytes)
    val root = Json.parseToJsonElement(text).jsonObject
    val complete = root["complete"]?.jsonPrimitive?.booleanOrNull == true
    val engine = engineName(root["userAgent"]?.jsonPrimitive?.contentOrNull ?: "")
    val logged = if (complete) {
        Files.write(reports.resolve("$engine.json"), bytes)
        "$name and $engine.json (${bytes.size} bytes)"
    } else {
        "$name (${bytes.size} bytes)"
    }
    return SavedReport(name, logged)
}

private fun ensureHost(text: String): String {
    val root = Json.parseToJsonElement(text).jsonObject
    val host = root["host"]?.jsonPrimitive?.contentOrNull
    if (!host.isNullOrBlank()) return text
    val withHost = linkedMapOf<String, JsonElement>("host" to JsonPrimitive(machineHost))
    withHost.putAll(root)
    return reportJson.encodeToString(JsonObject(withHost))
}

private fun engineName(userAgent: String): String = when {
    "Firefox/" in userAgent -> "firefox"
    "Chrome/" in userAgent || "Chromium/" in userAgent -> "chrome"
    "Safari/" in userAgent -> "safari"
    else -> "browser"
}

private fun contentType(file: Path): String = when (file.fileName.toString().substringAfterLast('.', "")) {
    "html" -> "text/html; charset=utf-8"
    "js" -> "text/javascript; charset=utf-8"
    "wasm" -> "application/wasm"
    "json" -> "application/json"
    "wgsl" -> "text/plain; charset=utf-8"
    "map" -> "application/json"
    else -> "application/octet-stream"
}
