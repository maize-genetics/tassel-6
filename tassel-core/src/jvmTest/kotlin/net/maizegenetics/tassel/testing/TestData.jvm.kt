package net.maizegenetics.tassel.testing

import java.io.File

internal actual fun environmentVariable(name: String): String? = System.getenv(name)

internal actual fun readFile(absolutePath: String): ByteArray = File(absolutePath).readBytes()
