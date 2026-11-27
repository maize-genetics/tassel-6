package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenManifest
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.math.abs
import kotlin.test.fail

/**
 * The departures from TASSEL 5 in the GLM goldens, read from the `[glm-departures]` section of `testdata/golden/MANIFEST` (ADR-0029). The
 * full-precision reports are compared with the golden at its MANIFEST tolerance after the cells that file names are
 * set aside: a `skip` cell is replaced by the golden's text, and a cell within a `bound` of TASSEL 5's own rounding
 * noise likewise. A `skip` rule that sets aside no cell that differs is stale and fails.
 */
internal object GlmDepartures {
    private class Skip(val id: String, val file: String, val rows: String, val columns: Set<String>)

    private val skips = ArrayList<Skip>()
    private var epsilonMultiple = 0.0
    private var estimateBound = 0.0

    init {
        val section = goldenManifest.section("glm-departures")
        for ((key, value) in section.properties) {
            val fields = value.split('|').map { it.trim() }
            when {
                key.startsWith("skip.") -> skips += Skip(fields[0], fields[1], fields[2], fields[3].split(',').toSet())
                key == "bound.cancellation_eps_multiple" -> epsilonMultiple = fields[0].toDouble()
                key == "bound.estimate_absolute" -> estimateBound = fields[0].toDouble()
                key == "description" -> Unit
                else -> error("Unknown departure rule '$key'")
            }
        }
    }

    /**
     * D4's bound, `bound.cancellation_eps_multiple`: the multiple of eps of the error sum of squares within which a sum of
     * squares that is zero in exact arithmetic is rounding noise. `GlmFastTest` builds its kernel-against-kernel bound
     * from it, so the number has one source.
     */
    val cancellationEpsMultiple: Double get() = epsilonMultiple

    /** The ids of the rules, for the departure table. */
    val ids: Set<String> get() = skips.map { it.id }.toSet()

    fun assertFullMatches(label: String, kind: String, actual: String) {
        val file = "$label.$kind.full.txt"
        val path = "golden/$file"
        val golden = normalizeLineEndings(goldenText(path)).split('\n').map { it.split('\t') }
        val rows = normalizeLineEndings(actual).split('\n').map { it.split('\t') }
        val header = golden[0]
        val applying = skips.filter { it.file == file }
        val exercised = HashSet<Skip>()
        val adjusted = rows.mapIndexed { r, row ->
            if (r == 0 || r >= golden.size || row.size != golden[r].size) return@mapIndexed row
            val marker = golden[r].getOrNull(1)
            val noise = rowNoise(header, golden[r])
            row.mapIndexed { c, cell ->
                val expected = golden[r][c]
                if (cell == expected) return@mapIndexed cell
                val name = header[c]
                val rule = applying.firstOrNull {
                    name in it.columns && (it.rows == "*" || it.rows == "marker=$marker")
                }
                if (rule != null) {
                    if (differs(expected, cell)) exercised += rule
                    return@mapIndexed expected
                }
                if (withinNoise(name, expected, cell, noise)) expected else cell
            }
        }
        for (rule in applying) {
            if (rule !in exercised) fail("Departure ${rule.id} in $file names no cell that differs: the rule is stale")
        }
        assertMatchesGolden(path, adjusted.joinToString("\n") { it.joinToString("\t") })
    }

    private fun differs(expected: String, actual: String): Boolean {
        val e = expected.toDoubleOrNull()
        val a = actual.toDoubleOrNull()
        if (e == null || a == null) return expected != actual
        if (e.isNaN() && a.isNaN()) return false
        return abs(e - a) > 1e-8 * abs(e)
    }

    /** The error degrees of freedom and error sum of squares of a statistics row. */
    private class RowNoise(val errorDf: Double, val errorSs: Double)

    private fun rowNoise(header: List<String>, row: List<String>): RowNoise? {
        val df = header.indexOf("error_df")
        val ms = header.indexOf("error_MS")
        if (df < 0 || ms < 0 || row.size != header.size) return null
        val errorDf = row[df].toDoubleOrNull() ?: return null
        val errorMs = row[ms].toDoubleOrNull() ?: return null
        return RowNoise(errorDf, errorDf * errorMs)
    }

    /**
     * Whether [actual] is within the rounding noise of TASSEL 5's [expected] (D4): a few eps of the error sum of squares,
     * as a difference of residual sums of squares carries, for the statistics made from one, and an absolute bound for a
     * class estimate that is itself nearly zero.
     */
    private fun withinNoise(column: String, expected: String, actual: String, noise: RowNoise?): Boolean {
        val e = expected.toDoubleOrNull() ?: return false
        val a = actual.toDoubleOrNull() ?: return false
        val delta = abs(e - a)
        val bound = epsilonMultiple * EPSILON
        return when (column) {
            // The noise of a sum of squares in F units is bound * SSE / MSE, which is bound * error df.
            "marker_F", "add_F" -> noise != null && delta <= bound * noise.errorDf
            // dom_F is a difference of two such differences, so it carries twice the noise.
            "dom_F" -> noise != null && delta <= 2 * bound * noise.errorDf
            "marker_Rsq" -> delta <= bound
            "marker_MS", "model_MS" -> noise != null && delta <= bound * noise.errorSs
            "Estimate" -> abs(e) < 1e-3 && delta <= estimateBound
            else -> false
        }
    }

    private const val EPSILON = 2.220446049250313E-16
}
