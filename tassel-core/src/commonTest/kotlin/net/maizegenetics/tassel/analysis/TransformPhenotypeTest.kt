package net.maizegenetics.tassel.analysis

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenTolerance
import kotlin.math.abs
import kotlin.math.ln
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * [transformPhenotype] and the `transform-phenotype` spec: TASSEL 5's `TransformDataPlugin` exports at a relative
 * tolerance, and hand-computed cases for what TASSEL 5 gets wrong or leaves undefined (ADR-0023).
 *
 * @see net.maizegenetics.analysis.numericaltransform.TransformDataPlugin
 */
class TransformPhenotypeTest {
    private val nan = Float.NaN

    private fun table(vararg attributes: net.maizegenetics.tassel.phenotype.PhenotypeAttribute): PhenotypeTable {
        val rows = attributes.first().size
        return PhenotypeTable(listOf(PhenotypeTable.taxa(List(rows) { "t$it" })) + attributes, "in")
    }

    private fun data(name: String, vararg values: Float) = NumericAttribute(name, AttributeType.data, values)

    private fun PhenotypeTable.floats(name: String) = (attribute(name) as NumericAttribute).toFloatArray().toList()

    private suspend fun run(table: PhenotypeTable, transform: PhenotypeTransform) = transformPhenotype(table, transform)

    private fun assertClose(expected: List<Float>, actual: List<Float>, relative: Double = 1e-6) {
        assertEquals(expected.size, actual.size)
        expected.zip(actual).forEachIndexed { index, (e, a) ->
            if (e.isNaN()) {
                assertTrue(a.isNaN(), "row $index: expected missing, actual $a")
            } else {
                assertTrue(abs(e - a) <= relative * abs(e), "row $index: expected $e, actual $a")
            }
        }
    }

    /** `(x - mean) / sd` with the n - 1 standard deviation, computed independently of the code under test. */
    private fun zScores(values: List<Double>): List<Double> {
        val mean = values.sum() / values.size
        val sd = sqrt(values.sumOf { (it - mean) * (it - mean) } / (values.size - 1))
        return values.map { (it - mean) / sd }
    }

    // --- TASSEL 5 goldens ---------------------------------------------------------------------

    @Test
    fun theGoldensAreComparedAsPhenotypesWithARelativeTolerance() {
        GOLDENS.forEach { (golden, _) ->
            assertEquals(Tolerance.Phenotype(relative = 1e-6, absolute = 1e-6), goldenTolerance(golden))
        }
    }

    @Test
    fun matchesTassel5sTransformsOfTheTutorialTraits() = runTest {
        val traits = PhenotypeFiles.read(Buffer().apply { write(TestData.readBytes("tutorial/mdp_traits.txt")) }).table
        GOLDENS.filter { it.value.second == "tutorial/mdp_traits.txt" }.forEach { (golden, spec) ->
            assertMatchesGolden(golden, run(traits, spec.first).table.toDelimitedString())
        }
    }

    @Test
    fun matchesTassel5StandardizingWithinReplicates() = runTest {
        val input = PhenotypeFiles.read(
            Buffer().apply {
                write(TestData.readBytes("phenotype/earht_with_rep.txt"))
            },
        ).table
        val (golden, spec) = GOLDENS.entries.single { it.value.second == "phenotype/earht_with_rep.txt" }
        assertMatchesGolden(golden, run(input, spec.first).table.toDelimitedString())
    }

    // --- Math ---------------------------------------------------------------------------------

    @Test
    fun takesNaturalAndBaseTwoAndBaseTenLogarithmsInDoubleAndNarrowsToFloat() = runTest {
        val input = table(data("x", 1f, 8f, 100f, nan))
        val natural = run(input, PhenotypeTransform(log = LogBase.natural)).table.floats("x")
        assertEquals(listOf(0f, ln(8.0).toFloat(), ln(100.0).toFloat()), natural.take(3))
        assertTrue(natural[3].isNaN())
        assertEquals(
            listOf(0f, 3f, 6.643856f),
            run(input, PhenotypeTransform(log = LogBase.base_2)).table.floats("x").take(3),
        )
        assertEquals(
            listOf(0f, 0.90309f, 2f),
            run(input, PhenotypeTransform(log = LogBase.base_10)).table.floats("x").take(3),
        )
    }

    @Test
    fun powerWinsOverLog() = runTest {
        val input = table(data("x", 2f, 3f))
        val result = run(input, PhenotypeTransform(log = LogBase.natural, power = 2.0))
        assertEquals(listOf(4f, 9f), result.table.floats("x"))
        assertEquals(emptyList(), result.warnings)
    }

    @Test
    fun logOfZeroAndOfANegativeValueIsMissingWithAWarning() = runTest {
        val result = run(table(data("x", 1f, 0f, -1f, 10f)), PhenotypeTransform(log = LogBase.natural))
        val values = result.table.floats("x")
        assertEquals(0f, values[0])
        assertTrue(values[1].isNaN() && values[2].isNaN())
        assertEquals(ln(10.0).toFloat(), values[3])
        assertEquals(listOf("Trait 'x': 2 values had no finite result, so they are missing"), result.warnings)
    }

    @Test
    fun aPowerWithNoRealResultOrOneTooBigForAFloatIsMissing() = runTest {
        val root = run(table(data("x", 4f, -4f)), PhenotypeTransform(power = 0.5))
        assertEquals(listOf(2f), root.table.floats("x").take(1))
        assertTrue(root.table.floats("x")[1].isNaN())
        assertEquals(1, root.warnings.size)
        val huge = run(table(data("x", 1e30f, 2f)), PhenotypeTransform(power = 2.0))
        assertTrue(huge.table.floats("x")[0].isNaN())
        assertEquals(4f, huge.table.floats("x")[1])
        assertContains(huge.warnings.single(), "1 values")
    }

    @Test
    fun missingInputsStayMissingWithoutAWarning() = runTest {
        val input = table(data("x", nan, 2f, nan))
        // Math.pow(NaN, 0) is 1; a missing value must stay missing.
        val result = run(input, PhenotypeTransform(power = 0.0))
        val values = result.table.floats("x")
        assertTrue(values[0].isNaN() && values[2].isNaN())
        assertEquals(1f, values[1])
        assertEquals(emptyList(), result.warnings)
        val standardized = run(table(data("x", nan, 1f, 3f)), PhenotypeTransform(standardize = true))
        assertTrue(standardized.table.floats("x")[0].isNaN())
        assertEquals(emptyList(), standardized.warnings)
    }

    @Test
    fun standardizesWithTheNMinusOneStandardDeviationSkippingMissingValues() = runTest {
        val result = run(table(data("x", 1f, 2f, nan, 3f, 4f)), PhenotypeTransform(standardize = true))
        val sd = sqrt(5.0 / 3.0)
        val z = listOf(-1.5, -0.5, 0.5, 1.5).map { (it / sd).toFloat() }
        assertClose(listOf(z[0], z[1], nan, z[2], z[3]), result.table.floats("x"))
        assertEquals(emptyList(), result.warnings)
    }

    @Test
    fun standardizesTheLogOrPowerResultNotTheOriginal() = runTest {
        val result = run(table(data("x", 1f, 10f, 100f)), PhenotypeTransform(log = LogBase.base_10, standardize = true))
        assertEquals(listOf(-1f, 0f, 1f), result.table.floats("x"))
    }

    @Test
    fun standardizingOneValueOrAConstantTraitGivesMissingWithAWarning() = runTest {
        val one = run(table(data("x", 5f, nan)), PhenotypeTransform(standardize = true))
        assertTrue(one.table.floats("x").all { it.isNaN() })
        assertEquals(listOf("Trait 'x': 1 values had no finite result, so they are missing"), one.warnings)

        val constant = run(table(data("x", 2f, 2f, 2f)), PhenotypeTransform(standardize = true))
        assertTrue(constant.table.floats("x").all { it.isNaN() })
        assertContains(constant.warnings.single(), "3 values")
    }

    @Test
    fun aLogThatFailsIsNotCountedAgainWhenStandardizing() = runTest {
        val result =
            run(table(data("x", 0f, 1f, 10f, 100f)), PhenotypeTransform(log = LogBase.natural, standardize = true))
        assertTrue(result.table.floats("x")[0].isNaN())
        assertContains(result.warnings.single(), "1 values")
    }

    // --- Groups -------------------------------------------------------------------------------

    @Test
    fun standardizesWithinGroupsOfTwoFactors() = runTest {
        val input = table(
            FactorAttribute.of("f1", listOf("a", "a", "b", "b", "a", "b")),
            FactorAttribute.of("f2", listOf("x", "x", "x", "y", "y", "y")),
            data("v", 1f, 3f, 5f, 2f, 7f, 6f),
        )
        val result = run(input, PhenotypeTransform(factors = listOf("f1", "f2"), standardize = true))
        val h = (1 / sqrt(2.0)).toFloat()
        // (a,x) = rows 0 and 1, (b,y) = rows 3 and 5; (b,x) and (a,y) have one row each.
        assertClose(listOf(-h, h, nan, -h, nan, h), result.table.floats("v"))
        assertEquals(listOf("Trait 'v': 2 values had no finite result, so they are missing"), result.warnings)
        // One factor groups differently: a = rows 0, 1, 4 and b = rows 2, 3, 5.
        val byOne = run(input, PhenotypeTransform(factors = listOf("f1"), standardize = true)).table.floats("v")
        val a = zScores(listOf(1.0, 3.0, 7.0))
        val b = zScores(listOf(5.0, 2.0, 6.0))
        assertClose(listOf(a[0], a[1], b[0], b[1], a[2], b[2]).map { it.toFloat() }, byOne)
    }

    @Test
    fun aRowWithAMissingFactorBelongsToNoGroupSoItIsMissing() = runTest {
        val input = table(
            FactorAttribute.of("f", listOf("a", null, "a", "a")),
            data("v", 1f, 2f, 3f, 5f),
        )
        val result = run(input, PhenotypeTransform(factors = listOf("f"), standardize = true))
        val values = result.table.floats("v")
        assertTrue(values[1].isNaN())
        assertTrue(values.filterIndexed { i, _ -> i != 1 }.none { it.isNaN() })
        assertEquals(listOf("Trait 'v': 1 have a missing factor value, so they are missing"), result.warnings)
    }

    // --- Selection ----------------------------------------------------------------------------

    @Test
    fun transformsEveryDataAttributeByDefaultAndOnlyThose() = runTest {
        val factor = FactorAttribute.of("f", listOf("a", "b"))
        val covariate = NumericAttribute("q", AttributeType.covariate, floatArrayOf(4f, 9f))
        val input = table(factor, data("x", 4f, 9f), covariate, data("y", 16f, 25f))
        val all = run(input, PhenotypeTransform(power = 0.5)).table
        assertEquals(input.attributes.map { it.name }, all.attributes.map { it.name })
        assertEquals(listOf(2f, 3f), all.floats("x"))
        assertEquals(listOf(4f, 5f), all.floats("y"))
        assertEquals(covariate, all.attribute("q"))
        assertEquals(factor, all.attribute("f"))

        val named = run(input, PhenotypeTransform(traits = listOf("y"), power = 0.5)).table
        assertEquals(listOf(4f, 9f), named.floats("x"))
        assertEquals(listOf(4f, 5f), named.floats("y"))
    }

    @Test
    fun rejectsTraitsThatAreNotDataAndFactorsThatAreNotFactors() = runTest {
        val input = table(
            FactorAttribute.of("f", listOf("a", "b")),
            NumericAttribute("q", AttributeType.covariate, floatArrayOf(1f, 2f)),
            data("x", 1f, 2f),
        )
        assertContains(
            assertFailsWith<IllegalArgumentException> {
                run(input, PhenotypeTransform(traits = listOf("q"), power = 2.0))
            }.message!!,
            "only data attributes",
        )
        assertFailsWith<IllegalArgumentException> { run(input, PhenotypeTransform(traits = listOf("f"), power = 2.0)) }
        assertFailsWith<IllegalArgumentException> {
            run(input, PhenotypeTransform(traits = listOf("nope"), power = 2.0))
        }
        assertContains(
            assertFailsWith<IllegalArgumentException> {
                run(input, PhenotypeTransform(factors = listOf("x"), standardize = true))
            }.message!!,
            "not a factor",
        )
    }

    @Test
    fun aTransformNeedsSomethingToDoAndFactorsNeedStandardize() {
        assertFailsWith<IllegalArgumentException> { PhenotypeTransform() }
        assertFailsWith<IllegalArgumentException> { PhenotypeTransform(factors = listOf("f"), power = 2.0) }
        assertFailsWith<IllegalArgumentException> { PhenotypeTransform(power = Double.NaN) }
    }

    // --- Spec ---------------------------------------------------------------------------------

    @Test
    fun specDeclaresTassel5sNamesAndRunsThroughTheEngine() = runTest {
        assertEquals("TransformDataPlugin", transformPhenotypeOperation.legacyPlugin)
        listOf("traits", "factors", "log", "power", "standardize").forEach {
            assertTrue(
                transformPhenotypeOperation.parameter(it) != null,
                it,
            )
        }
        assertEquals("factor", transformPhenotypeOperation.parameter("factors")?.legacy)
        assertEquals("traits", transformPhenotypeOperation.parameter("traits")?.legacy)
        assertEquals(false, transformPhenotypeOperation.parameter("standardize")?.default)

        val input = table(data("x", 1f, 100f), data("y", 1f, 10f))
        val engine = LocalEngine(OperationRegistry(transformPhenotypeOperation))
        val handle = engine.run(
            RunRequest(
                "transform-phenotype",
                engine.put(input),
                mapOf("traits" to "x", "log" to "base_10", "standardize" to "false"),
            ),
        )
        val result = engine.get(handle) as PhenotypeTransformResult
        assertEquals(listOf(0f, 2f), result.table.floats("x"))
        assertEquals(listOf(1f, 10f), result.table.floats("y"))
    }

    private companion object {
        /** Golden to its input and the TASSEL 5 arguments it came from. */
        val GOLDENS = mapOf(
            "golden/phenotype-transform/mdp_traits.log-natural.txt" to
                (PhenotypeTransform(log = LogBase.natural) to "tutorial/mdp_traits.txt"),
            "golden/phenotype-transform/mdp_traits.log-base_10.txt" to
                (PhenotypeTransform(log = LogBase.base_10) to "tutorial/mdp_traits.txt"),
            "golden/phenotype-transform/mdp_traits.power-0.5.txt" to
                (PhenotypeTransform(power = 0.5) to "tutorial/mdp_traits.txt"),
            "golden/phenotype-transform/mdp_traits.standardize.txt" to
                (PhenotypeTransform(standardize = true) to "tutorial/mdp_traits.txt"),
            "golden/phenotype-transform/earht_with_rep.standardize-rep.txt" to
                (PhenotypeTransform(factors = listOf("rep"), standardize = true) to "phenotype/earht_with_rep.txt"),
        )
    }
}
