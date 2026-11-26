package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

/**
 * The legacy writers against TASSEL 5's own GLM reports (#188): records parsed from the full-precision twin of each
 * report must write the rounded report in [NumberForm.Tassel] and the full-precision twin itself in
 * [NumberForm.Full], byte for byte.
 */
class AssociationReportGoldenTest {
    private val provenance = AnalysisProvenance(operation = "glm", version = "test")

    private fun text(path: String) = normalizeLineEndings(goldenText("golden/$path"))

    private fun check(directory: String, name: String): Int {
        val statsFull = text("$directory/$name.stats.full.txt")
        val effectsFull = text("$directory/$name.effects.full.txt")
        val layout = detectLayout(statsFull.substringBefore('\n'))
        val parsed = parseGlm(statsFull, effectsFull)
        val table = AssociationTable.of(parsed.records, parsed.effects, provenance)

        fun stats(form: NumberForm) = buildString {
            table.writeGlmStatistics(this, layout.layout, form, layout.permutations, layout.appendAddDom)
        }

        fun effects(form: NumberForm) = buildString { table.writeGlmEffects(this, layout.layout, form) }
        assertSame(text("$directory/$name.stats.txt"), stats(NumberForm.Tassel), "$directory/$name stats, rounded")
        assertSame(statsFull, stats(NumberForm.Full), "$directory/$name stats, full")
        assertSame(
            text("$directory/$name.effects.txt"),
            effects(NumberForm.Tassel),
            "$directory/$name effects, rounded",
        )
        assertSame(effectsFull, effects(NumberForm.Full), "$directory/$name effects, full")
        return table.size + table.effectCount
    }

    /** Fails at the first cell that differs, not with two whole reports. */
    private fun assertSame(expected: String, actual: String, label: String) {
        if (expected == actual) return
        val want = expected.split('\n')
        val got = actual.split('\n')
        for (line in want.indices) {
            if (want[line] == got.getOrNull(line)) continue
            val wantCells = want[line].split('\t')
            val gotCells = got.getOrNull(line)?.split('\t').orEmpty()
            val column = wantCells.indices.firstOrNull { wantCells[it] != gotCells.getOrNull(it) } ?: 0
            val expectedCell = wantCells[column]
            fail(
                "$label, line ${line + 1} column ${column + 1}: expected '$expectedCell' but was '${gotCells.getOrNull(
                    column,
                )}'",
            )
        }
        fail("$label: ${want.size} lines expected, ${got.size} written")
    }

    @Test
    fun reproducesEveryGlmGolden() {
        var cells = 0
        for ((directory, names) in RUNS) for (name in names) cells += check(directory, name)
        assertTrue(cells > 100_000, "only $cells rows were checked")
    }

    @Test
    fun usesBothLayouts() {
        val full = text("glm-refprob/mdp_traits.stats.full.txt")
        assertEquals(GlmLayout.ReferenceProbability, detectLayout(full.substringBefore('\n')).layout)
        val add = text("glm-options/append-add-dom.stats.full.txt")
        val layout = detectLayout(add.substringBefore('\n'))
        assertEquals(GlmLayout.Genotype, layout.layout)
        assertTrue(layout.appendAddDom)
    }

    companion object {
        val RUNS = listOf(
            "glm-tutorial" to listOf("mdp_traits"),
            "glm-tutorial-maf" to listOf("mdp_traits"),
            "glm-no-covariates" to listOf("mdp_traits"),
            "glm-refprob" to listOf("mdp_traits"),
            "glm-options" to listOf("append-add-dom", "biallelic-only", "max-p-0.01", "min-class-size-5"),
            "glm-synthetic" to listOf(
                "append-add-dom",
                "biallelic-only",
                "collinear",
                "factor-only",
                "min-class-size-3",
                "model",
            ),
        )
    }
}
