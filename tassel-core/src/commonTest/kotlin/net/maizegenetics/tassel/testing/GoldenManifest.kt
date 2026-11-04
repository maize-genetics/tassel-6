package net.maizegenetics.tassel.testing

/**
 * Parsed `testdata/golden/MANIFEST`, written by `scripts/generate-goldens.sh` (ADR-0009).
 *
 * The format is `[section]` headers followed by `key = value` lines; `#` starts a comment line.
 * A key that starts with `golden/` is an output file and its value is that file's [Tolerance].
 * `command` may repeat within a section.
 */
class GoldenManifest(val sections: Map<String, Section>) {
    class Section(
        val name: String,
        val properties: Map<String, String>,
        val commands: List<String>,
        val outputs: Map<String, Tolerance>,
    ) {
        operator fun get(key: String): String = properties[key] ?: error("[$name] has no '$key'")
    }

    fun section(name: String): Section = sections[name] ?: error("MANIFEST has no [$name] section")

    /** All golden output files, relative to `testdata/`, with their tolerances. */
    val outputs: Map<String, Tolerance> by lazy { buildMap { sections.values.forEach { putAll(it.outputs) } } }

    fun tolerance(path: String): Tolerance = outputs[path] ?: error("MANIFEST does not list $path")

    companion object {
        const val PATH = "golden/MANIFEST"

        fun load(): GoldenManifest = parse(TestData.readText(PATH))

        fun parse(text: String): GoldenManifest {
            val sections = linkedMapOf<String, Section>()
            var name: String? = null
            var properties = linkedMapOf<String, String>()
            var commands = mutableListOf<String>()
            var outputs = linkedMapOf<String, Tolerance>()

            fun finish() {
                val current = name ?: return
                require(current !in sections) { "MANIFEST repeats section [$current]" }
                sections[current] = Section(current, properties, commands, outputs)
            }

            for ((index, raw) in text.lines().withIndex()) {
                val line = raw.trim()
                if (line.isEmpty() || line.startsWith("#")) continue
                if (line.startsWith("[") && line.endsWith("]")) {
                    finish()
                    name = line.substring(1, line.length - 1)
                    properties = linkedMapOf()
                    commands = mutableListOf()
                    outputs = linkedMapOf()
                    continue
                }
                val separator = line.indexOf(" = ")
                require(name != null && separator > 0) { "MANIFEST line ${index + 1} is not 'key = value': $raw" }
                val key = line.substring(0, separator)
                val value = line.substring(separator + 3)
                when {
                    key == "command" -> commands += value
                    key.startsWith("golden/") -> outputs[key] = Tolerance.parse(value)
                    else -> properties[key] = value
                }
            }
            finish()
            return GoldenManifest(sections)
        }
    }
}

/** How a golden file is compared; see the legend at the top of `testdata/golden/MANIFEST`. */
sealed interface Tolerance {
    /** Whether the data rows may be compared in any order, which is the `rows-unordered` keyword. */
    val rowsUnordered: Boolean

    /** Identical text after normalizing line endings. */
    data class Exact(override val rowsUnordered: Boolean = false) : Tolerance

    /** Headers, row labels, and layout identical; numbers are checked in the `.full.txt` twin. */
    data class Format(override val rowsUnordered: Boolean = false) : Tolerance

    /** Every number within [value] absolutely; [upToSign] allows negating whole columns. */
    data class Absolute(
        val value: Double,
        val upToSign: Boolean = false,
        override val rowsUnordered: Boolean = false,
    ) : Tolerance

    /** Every number within [value] relative to the golden; [upToSign] allows negating whole columns. */
    data class Relative(
        val value: Double,
        val upToSign: Boolean = false,
        override val rowsUnordered: Boolean = false,
    ) : Tolerance

    companion object {
        /** Keywords that may follow the tolerance itself, in any order. */
        private val MODIFIERS = setOf("up-to-sign", "rows-unordered")

        fun parse(text: String): Tolerance {
            val words = text.split(' ').filter { it.isNotEmpty() }
            val modifiers = words.takeLastWhile { it in MODIFIERS }
            val terms = words.dropLast(modifiers.size)
            val upToSign = "up-to-sign" in modifiers
            val rowsUnordered = "rows-unordered" in modifiers
            fun number(): Double = terms.getOrNull(1)?.toDoubleOrNull() ?: error("Tolerance needs a number: '$text'")
            val tolerance = when (terms.firstOrNull()) {
                "exact" -> Exact(rowsUnordered)
                "format" -> Format(rowsUnordered)
                "abs" -> Absolute(number(), upToSign, rowsUnordered)
                "rel" -> Relative(number(), upToSign, rowsUnordered)
                else -> error("Unknown tolerance: '$text'")
            }
            val expectedTerms = if (tolerance is Absolute || tolerance is Relative) 2 else 1
            require(
                terms.size == expectedTerms &&
                    modifiers.size == modifiers.distinct().size &&
                    (!upToSign || expectedTerms == 2),
            ) { "Malformed tolerance: '$text'" }
            return tolerance
        }
    }
}
