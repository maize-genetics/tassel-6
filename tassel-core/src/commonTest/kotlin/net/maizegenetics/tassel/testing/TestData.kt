package net.maizegenetics.tassel.testing

/**
 * Reads files under the repository's `testdata/` directory from `commonTest`, on every target.
 *
 * Gradle passes the absolute path of `testdata/` to each test runtime in the [ENV_VAR]
 * environment variable (see `tassel-core/build.gradle.kts`).
 */
object TestData {
    const val ENV_VAR = "TASSEL_TESTDATA_DIR"

    val root: String by lazy {
        environmentVariable(ENV_VAR)?.trimEnd('/')
            ?: error("$ENV_VAR is not set; run tests through Gradle so it points at testdata/")
    }

    /** Bytes of [path], relative to `testdata/` (for example `tutorial/mdp_genotype.hmp.txt`). */
    fun readBytes(path: String): ByteArray = readFile("$root/$path")

    /** Text of [path], relative to `testdata/`, decoded as UTF-8 with line endings unchanged. */
    fun readText(path: String): String = readBytes(path).decodeToString()
}

internal expect fun environmentVariable(name: String): String?

internal expect fun readFile(absolutePath: String): ByteArray
