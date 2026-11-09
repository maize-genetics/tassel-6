package net.maizegenetics.tassel.testing

/**
 * Parsed `testdata/htslib/MANIFEST`, written by `scripts/generate-htslib-fixtures.sh`
 * (ADR-0021, tassel-6#80).
 *
 * The format matches `testdata/golden/MANIFEST`: `[section]` headers, `key = value` lines,
 * and `#` comments. `command` and `file` may repeat. A `file` value is a path relative to
 * `testdata/`.
 */
class HtslibManifest(val sections: Map<String, Section>) {
    class Section(
        val name: String,
        val properties: Map<String, String>,
        val commands: List<String>,
        val files: List<String>,
    ) {
        operator fun get(key: String): String = properties[key] ?: error("[$name] has no '$key'")
    }

    fun section(name: String): Section = sections[name] ?: error("MANIFEST has no [$name] section")

    /** Every output file, relative to `testdata/`, in manifest order. */
    val files: List<String> by lazy { sections.values.flatMap { it.files } }

    companion object {
        const val PATH = "htslib/MANIFEST"

        fun load(): HtslibManifest = parse(TestData.readText(PATH))

        fun parse(text: String): HtslibManifest {
            val sections = linkedMapOf<String, Section>()
            var name: String? = null
            var properties = linkedMapOf<String, String>()
            var commands = mutableListOf<String>()
            var files = mutableListOf<String>()

            fun finish() {
                val current = name ?: return
                require(current !in sections) { "MANIFEST repeats section [$current]" }
                sections[current] = Section(current, properties, commands, files)
            }

            for ((index, raw) in text.lines().withIndex()) {
                val line = raw.trim()
                if (line.isEmpty() || line.startsWith("#")) continue
                if (line.startsWith("[") && line.endsWith("]")) {
                    finish()
                    name = line.substring(1, line.length - 1)
                    properties = linkedMapOf()
                    commands = mutableListOf()
                    files = mutableListOf()
                    continue
                }
                val separator = line.indexOf(" = ")
                require(name != null && separator > 0) {
                    "MANIFEST line ${index + 1} is not 'key = value': $raw"
                }
                val key = line.substring(0, separator)
                val value = line.substring(separator + 3)
                when (key) {
                    "command" -> commands += value
                    "file" -> files += value
                    else -> {
                        require(key !in properties) { "MANIFEST section [$name] repeats '$key'" }
                        properties[key] = value
                    }
                }
            }
            finish()
            return HtslibManifest(sections)
        }
    }
}
