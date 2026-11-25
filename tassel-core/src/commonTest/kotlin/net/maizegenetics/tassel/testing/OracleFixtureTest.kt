package net.maizegenetics.tassel.testing

import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Smoke test for `testdata/oracle/` (tassel-6#183): the MANIFEST parses with the golden syntax, records
 * its tools, and one NIST fixture and one mpmath fixture read on every target. The statistics kernels
 * that consume these files (#184 to #186) compare against them; this test only proves they load.
 */
class OracleFixtureTest {
    @Test
    fun manifestRecordsToolsAndListsEveryFixture() {
        val manifest = OracleManifest.load()
        val tools = manifest.section("tools")
        for (key in listOf("python_image", "r_image")) {
            assertTrue(tools[key].contains("@sha256:"), "$key is pinned by digest: ${tools[key]}")
        }
        assertEquals("1.3.0", tools["mpmath_version"])
        assertTrue(tools["r_version"].startsWith("R version 4.4.2"), tools["r_version"])
        assertEquals(Tolerance.Relative(1e-10), manifest.tolerance("oracle/nist/Norris.certified.txt"))
        assertEquals(Tolerance.Exact(), manifest.tolerance("oracle/nist/Norris.data.tsv"))
        for ((path, _) in manifest.outputs) {
            assertTrue(TestData.readText(path).isNotEmpty(), "$path is empty")
        }
        assertEquals(11 * 2 + 1 + 4, manifest.outputs.size)
    }

    @Test
    fun readsANistFixture() {
        val certified = keyValues(TestData.readOracleText("nist/Norris.certified.txt"))
        assertEquals("36", certified.getValue("observations"))
        assertEquals("-0.262323073774029 0.232818234301152", certified.getValue("B0"))
        assertEquals("0.884796396144373", certified.getValue("residual_sd"))
        val rows = TestData.readOracleText("nist/Norris.data.tsv").lines().filter { it.isNotEmpty() }
        assertEquals("y\tx", rows.first())
        assertEquals(36, rows.size - 1)
        assertEquals("0.1\t0.2", rows[1])
    }

    @Test
    fun readsAnMpmathFixture() {
        val rows = TestData.readOracleText("tails/log-sf.tsv").lines().filter { it.isNotEmpty() }
        assertEquals(listOf("dist", "df1", "df2", "x", "log_sf"), rows.first().split('\t'))
        val first = rows[1].split('\t')
        assertEquals(listOf("chisq", "1", "-"), first.take(3))
        // Chi-square with 1 df at x = 0.8108148784511222: P(X > x) = 0.3678..., so log p is about -1.
        assertTrue(abs(first[4].toDouble() + 1.0) < 1e-8, rows[1])
        val deepest = rows.drop(1).minOf { it.split('\t')[4].toDouble() }
        assertTrue(deepest < -49_000.0, "deepest log p is $deepest")
    }

    private fun keyValues(text: String): Map<String, String> =
        text.lines().filter { " = " in it }.associate { it.substringBefore(" = ") to it.substringAfter(" = ") }
}
