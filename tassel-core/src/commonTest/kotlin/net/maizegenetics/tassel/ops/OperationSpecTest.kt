package net.maizegenetics.tassel.ops

import kotlinx.coroutines.test.runTest
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

class OperationSpecTest {
    @Test
    fun parametersKeepDeclarationOrderAndMetadata() {
        assertEquals(
            listOf("factor", "repeats", "limit", "negate", "rounding", "report"),
            scaleValues.parameters.map {
                it.name
            },
        )
        val factor = scaleValues.parameter("factor")!!
        assertEquals("scaleFactor", factor.legacy)
        assertEquals("Multiplier", factor.help)
        assertEquals(1.0, factor.default)
        assertEquals(false, factor.required)
    }

    @Test
    fun defaultsFillUnsuppliedValues() {
        val values = scaleValues.bind()
        assertEquals(
            mapOf(
                "factor" to 1.0,
                "repeats" to 1,
                "limit" to null,
                "negate" to false,
                "rounding" to Rounding.NONE,
                "report" to null,
            ),
            values.toMap(),
        )
        assertEquals(values, scaleValues.bind(mapOf("factor" to null, "repeats" to "")))
    }

    @Test
    fun textValuesAreParsed() {
        val values = scaleValues.bind(
            mapOf(
                "factor" to "2.5",
                "repeats" to "3",
                "limit" to "1e2",
                "negate" to "TRUE",
                "rounding" to "FLOOR",
                "report" to "out/report.txt",
            ),
        )
        assertEquals(listOf(2.5, 3, 100, true, Rounding.FLOOR, "out/report.txt"), values.toMap().values.toList())
    }

    @Test
    fun typedValuesAreAccepted() {
        val values = scaleValues.bind(
            mapOf(
                "factor" to 2,
                "repeats" to 2,
                "negate" to true,
                "rounding" to Rounding.CEILING,
            ),
        )
        assertEquals(2.0, values["factor"])
        assertEquals(Rounding.CEILING, values["rounding"])
    }

    @Test
    fun rangesAreInclusive() {
        for (ok in listOf(0.0, 10.0)) assertEquals(ok, scaleValues.bind(mapOf("factor" to ok))["factor"])
        for (ok in listOf("1", "5")) assertEquals(ok.toInt(), scaleValues.bind(mapOf("repeats" to ok))["repeats"])
        val problems = problemsOf { scaleValues.bind(mapOf("factor" to -0.25, "repeats" to "6", "limit" to 1001)) }
        assertEquals(listOf("factor", "repeats", "limit"), problems.map { it.parameter })
        assertTrue(problems.all { "outside range" in it.message }, problems.toString())
    }

    @Test
    fun badTypesAndTextAreRejected() {
        val problems = problemsOf {
            scaleValues.bind(
                mapOf(
                    "factor" to "abc",
                    "repeats" to "2.5",
                    "negate" to "yes",
                    "rounding" to "floor",
                    "limit" to 3.0,
                    "report" to 7,
                ),
            )
        }
        assertEquals(
            setOf("factor", "repeats", "negate", "rounding", "limit", "report"),
            problems.map {
                it.parameter
            }.toSet(),
        )
        assertTrue(problems.first { it.parameter == "rounding" }.message.contains("one of NONE, FLOOR, CEILING"))
    }

    @Test
    fun unknownNamesAreRejectedWithLegacyHint() {
        val e = assertFailsWith<InvalidParametersException> {
            scaleValues.bind(mapOf("bogus" to "1", "scaleFactor" to "2.0"))
        }
        assertEquals("scale-values", e.operationId)
        assertEquals(listOf("bogus", "scaleFactor"), e.problems.map { it.parameter })
        assertTrue("TASSEL 5 name of 'factor'" in e.problems[1].message, e.message)
    }

    @Test
    fun requiredParametersMustBeSupplied() {
        val path = readLines.parameter("path")!!
        assertTrue(path.required)
        assertNull(path.default)
        assertEquals("required parameter is missing", problemsOf { readLines.bind() }.single().message)
        assertEquals("a.txt", readLines.bind(mapOf("path" to "a.txt"))["path"])
    }

    @Test
    fun runReadsBoundValues() = runTest {
        val reports = mutableListOf<Progress>()
        val output = scaleValues.run(
            listOf(1.0, 2.0, 3.0),
            mapOf("factor" to "2", "repeats" to "2", "limit" to "2", "negate" to "true", "report" to "r.txt"),
        ) { reports += it }
        assertEquals(listOf(-4.0, -8.0), output)
        assertEquals(listOf(Progress(1, 2, "r.txt"), Progress(2, 2, "r.txt")), reports)
        assertEquals("read in.txt", readLines.run(Unit, mapOf("path" to "in.txt")))
    }

    @Test
    fun valuesBelongToTheirOperation() = runTest {
        assertFailsWith<IllegalArgumentException> { scaleValues.run(emptyList(), readLines.bind(mapOf("path" to "x"))) }
        assertFailsWith<IllegalArgumentException> {
            readLines.bind(mapOf("path" to "x"))[scaleValues.parameter("factor")!!]
        }
        assertFailsWith<IllegalArgumentException> { scaleValues.bind()["nope"] }
    }

    @Test
    fun invalidDeclarationsFail() {
        fun declare(block: OperationScope<Unit, Unit>.() -> Unit) = operation<Unit, Unit>(id = "bad", block = block)
        val finish: OperationScope<Unit, Unit>.() -> Unit = { run { _, _ -> } }

        assertEquals("ok", operation<Unit, Unit>(id = "ok", block = finish).id)
        assertFailsWith<IllegalArgumentException> { operation<Unit, Unit>(id = "Not_Kebab", block = finish) }
        assertFailsWith<IllegalArgumentException> { declare { } }
        assertFailsWith<IllegalArgumentException> {
            declare {
                finish()
                finish()
            }
        }
        assertFailsWith<IllegalArgumentException> { declare { double(default = 0.6, range = 0.0..0.5) } }
        assertFailsWith<IllegalArgumentException> { declare { int(range = 5..1) } }
        assertFailsWith<IllegalArgumentException> {
            declare {
                val a by int(default = 1, legacy = "same")
                val b by int(default = 2, legacy = "same")
                finish()
            }
        }
        assertFailsWith<IllegalArgumentException> {
            declare {
                @Suppress("ktlint:standard:property-naming", "LocalVariableName")
                val Upper by boolean()
                finish()
            }
        }
    }

    /** The ADR-0002 sketch, with a stand-in for `GenotypeTable`, compiles and runs as written. */
    @Test
    fun adrSketchWorks() = runTest {
        class Table(val maf: List<Double>)
        fun filterSites(input: Table, minMaf: Double, minCount: Int, progress: ProgressReporter): Table {
            progress.report(1, 1)
            return Table(input.maf.filter { it >= minMaf && minCount >= 0 })
        }

        val filterSitesOp = operation<Table, Table>(
            id = "filter-sites",
            legacyPlugin = "FilterSiteBuilderPlugin",
        ) {
            val minMaf by double(
                default = 0.0,
                range = 0.0..0.5,
                legacy = "siteMinAlleleFreq",
                help = "Minimum minor allele frequency",
            )
            val minCount by int(default = 0, legacy = "siteMinCount", help = "Minimum number of non-missing taxa")
            run { input, progress -> filterSites(input, minMaf(), minCount(), progress) }
        }

        val output = filterSitesOp.run(Table(listOf(0.01, 0.2, 0.4)), mapOf("minMaf" to "0.1"))
        assertEquals(listOf(0.2, 0.4), output.maf)
        assertEquals("minMaf", filterSitesOp.parameterByLegacyName("siteMinAlleleFreq")?.name)
    }

    private fun problemsOf(block: () -> Unit): List<ParameterProblem> =
        assertFailsWith<InvalidParametersException>(block = block).problems
}
