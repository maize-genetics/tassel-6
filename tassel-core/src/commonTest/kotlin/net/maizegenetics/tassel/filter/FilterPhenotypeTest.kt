package net.maizegenetics.tassel.filter

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertSame
import kotlin.test.assertTrue

/**
 * [filterPhenotype] and the `filter-phenotype` spec, against hand-written expected tables.
 *
 * @see net.maizegenetics.phenotype.PhenotypeBuilder
 * @see net.maizegenetics.analysis.filter.FilterTraitsPlugin
 */
class FilterPhenotypeTest {
    private val nan = Float.NaN

    private fun table(
        taxa: List<String>,
        env: List<String?>,
        yield: List<Float>,
        height: List<Float>,
        q: List<Float>,
    ) = PhenotypeTable(
        listOf(
            PhenotypeTable.taxa(taxa),
            FactorAttribute.of("env", env),
            NumericAttribute("yield", AttributeType.data, yield.toFloatArray()),
            NumericAttribute("height", AttributeType.data, height.toFloatArray()),
            NumericAttribute("q", AttributeType.covariate, q.toFloatArray()),
        ),
        "t",
    )

    /** Taxa A and B are replicated; row 3 has no env, row 2 has no yield, row 5 has no q. */
    private val input = table(
        taxa = listOf("A", "A", "B", "B", "C", "D"),
        env = listOf("e1", "e2", "e1", null, "e2", "e1"),
        yield = listOf(1f, 2f, nan, 4f, 5f, 6f),
        height = listOf(10f, 20f, 30f, 40f, 50f, 60f),
        q = listOf(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, nan),
    )

    private suspend fun filter(f: PhenotypeFilter) = filterPhenotype(input, f)

    private fun rows(t: PhenotypeTable, attribute: String) =
        (0 until t.numberOfObservations).map { t.attribute(attribute).format(it) }

    @Test
    fun defaultFilterReturnsTheInputItself() = runTest {
        assertSame(input, filter(PhenotypeFilter()))
        assertSame(input, filter(PhenotypeFilter(taxa = listOf("A", "B", "C", "D"))))
        assertSame(
            input,
            filter(
                PhenotypeFilter(
                    toData = listOf("yield"),
                    keepAttributes = input.attributes.map {
                        it.name
                    },
                ),
            ),
        )
    }

    @Test
    fun keepsAndDropsTaxaInTableOrderWithReplicates() = runTest {
        val kept = filter(PhenotypeFilter(taxa = listOf("C", "A")))
        assertEquals(listOf("A", "A", "C"), rows(kept, "Taxa"))
        assertEquals(listOf("10.0", "20.0", "50.0"), rows(kept, "height"))
        val dropped = filter(PhenotypeFilter(taxa = listOf("A", "B"), includeTaxa = false))
        assertEquals(listOf("C", "D"), rows(dropped, "Taxa"))
        assertEquals(listOf("e2", "e1"), rows(dropped, "env"))
    }

    @Test
    fun keepsAndDropsFactorLevelsAndTreatsMissingAsNoLevel() = runTest {
        val kept = filter(PhenotypeFilter(factor = "env", levels = listOf("e1")))
        assertEquals(listOf("A", "B", "D"), rows(kept, "Taxa"))
        assertEquals(listOf("e1"), (kept.attribute("env") as FactorAttribute).levels)
        // The row with a missing env is not at e1, so dropping e1 keeps it.
        val dropped = filter(PhenotypeFilter(factor = "env", levels = listOf("e1"), includeLevels = false))
        assertEquals(listOf("A", "B", "C"), rows(dropped, "Taxa"))
        assertEquals(listOf("e2", "?", "e2"), rows(dropped, "env"))
    }

    @Test
    fun rangeIsInclusiveAndMissingFailsIt() = runTest {
        assertEquals(
            listOf("2.0", "4.0", "5.0"),
            rows(filter(PhenotypeFilter(attribute = "yield", min = 2.0, max = 5.0)), "yield"),
        )
        assertEquals(listOf("1.0", "2.0"), rows(filter(PhenotypeFilter(attribute = "yield", max = 2.0)), "yield"))
        assertEquals(
            listOf("4.0", "5.0", "6.0"),
            rows(filter(PhenotypeFilter(attribute = "yield", min = 3.0)), "yield"),
        )
        // Bounds compare as Float, so 0.2 keeps the stored 0.2f.
        assertEquals(listOf("0.2", "0.3"), rows(filter(PhenotypeFilter(attribute = "q", min = 0.2, max = 0.3)), "q"))
    }

    @Test
    fun requireNonMissingNamedOrAll() = runTest {
        assertEquals(
            listOf("A", "A", "B", "C", "D"),
            rows(filter(PhenotypeFilter(requireNonMissing = listOf("yield"))), "Taxa"),
        )
        assertEquals(
            listOf("A", "A", "C"),
            rows(filter(PhenotypeFilter(requireNonMissing = listOf("yield", "env", "q"))), "Taxa"),
        )
        val all = filter(PhenotypeFilter(requireAllNonMissing = true))
        assertEquals(listOf("A", "A", "C"), rows(all, "Taxa"))
        assertEquals(listOf("10.0", "20.0", "50.0"), rows(all, "height"))
    }

    @Test
    fun keepsAndDropsAttributesInTableOrder() = runTest {
        val kept = filter(PhenotypeFilter(keepAttributes = listOf("q", "yield")))
        assertEquals(listOf("Taxa", "yield", "q"), kept.attributes.map { it.name })
        assertEquals(6, kept.numberOfObservations)
        val dropped = filter(PhenotypeFilter(dropAttributes = listOf("env", "height")))
        assertEquals(listOf("Taxa", "yield", "q"), dropped.attributes.map { it.name })
        val both =
            filter(PhenotypeFilter(keepAttributes = listOf("env", "yield", "q"), dropAttributes = listOf("yield")))
        assertEquals(listOf("Taxa", "env", "q"), both.attributes.map { it.name })
    }

    @Test
    fun keepTypesAlwaysKeepsTaxa() = runTest {
        assertEquals(
            listOf("Taxa", "yield", "height"),
            filter(PhenotypeFilter(keepTypes = setOf(AttributeType.data))).attributes.map {
                it.name
            },
        )
        assertEquals(
            listOf("Taxa", "env", "q"),
            filter(PhenotypeFilter(keepTypes = setOf(AttributeType.factor, AttributeType.covariate))).attributes.map {
                it.name
            },
        )
        assertEquals(listOf("Taxa"), filter(PhenotypeFilter(keepTypes = emptySet())).attributes.map { it.name })
    }

    @Test
    fun excludeLastDropsTheLastNonTaxaAttribute() = runTest {
        assertEquals(
            listOf("Taxa", "env", "yield", "height"),
            filter(PhenotypeFilter(excludeLast = true)).attributes.map {
                it.name
            },
        )
        val taxaOnly = PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("A"))))
        assertFailsWith<IllegalArgumentException> { filterPhenotype(taxaOnly, PhenotypeFilter(excludeLast = true)) }
    }

    @Test
    fun retypesNumericAttributesWithoutChangingValues() = runTest {
        val result = filter(PhenotypeFilter(toCovariate = listOf("yield"), toData = listOf("q")))
        assertEquals(AttributeType.covariate, result.attribute("yield").type)
        assertEquals(AttributeType.data, result.attribute("q").type)
        assertEquals(AttributeType.data, result.attribute("height").type)
        assertEquals(rows(input, "yield"), rows(result, "yield"))
        assertEquals(input.attributes.map { it.name }, result.attributes.map { it.name })
        assertEquals(6, result.numberOfObservations)
    }

    @Test
    fun retypingAnAttributeThatIsAlreadyThatTypeIsNotAChange() = runTest {
        assertSame(input, filter(PhenotypeFilter(toCovariate = listOf("q"), toData = listOf("yield"))))
    }

    @Test
    fun filtersCombine() = runTest {
        val result = filter(
            PhenotypeFilter(
                taxa = listOf("D"),
                includeTaxa = false,
                factor = "env",
                levels = listOf("e1"),
                attribute = "height",
                min = 15.0,
                requireNonMissing = listOf("yield"),
                dropAttributes = listOf("env"),
                toCovariate = listOf("yield"),
            ),
        )
        // Not D, at e1, height at least 15, yield present: only B's first row (30) fails yield, A's first fails height.
        assertEquals(0, result.numberOfObservations)
        val looser = filter(
            PhenotypeFilter(
                taxa = listOf("D"),
                includeTaxa = false,
                attribute = "height",
                min = 15.0,
                requireNonMissing = listOf("yield"),
                dropAttributes = listOf("env"),
                toCovariate = listOf("yield"),
            ),
        )
        assertEquals(listOf("A", "B", "C"), rows(looser, "Taxa"))
        assertEquals(listOf("Taxa", "yield", "height", "q"), looser.attributes.map { it.name })
        assertEquals(AttributeType.covariate, looser.attribute("yield").type)
    }

    @Test
    fun observationFiltersCanNameAttributesTheCallDrops() = runTest {
        val result = filter(PhenotypeFilter(factor = "env", levels = listOf("e2"), keepAttributes = listOf("height")))
        assertEquals(listOf("Taxa", "height"), result.attributes.map { it.name })
        assertEquals(listOf("20.0", "50.0"), rows(result, "height"))
    }

    @Test
    fun emptyingTheTableKeepsTheAttributes() = runTest {
        val result = filter(PhenotypeFilter(attribute = "height", min = 1000.0))
        assertEquals(0, result.numberOfObservations)
        assertEquals(input.attributes.map { it.name }, result.attributes.map { it.name })
        assertEquals(input.attributes.map { it.type }, result.attributes.map { it.type })
        assertEquals(
            PhenotypeTable.HEADER + "\n" + "taxa\tfactor\tdata\tdata\tcovariate\nTaxa\tenv\tyield\theight\tq\n",
            result.toDelimitedString(),
        )
    }

    @Test
    fun levelsUnusedByAnyObservationAreValid() = runTest {
        val withSpare = PhenotypeTable(
            listOf(PhenotypeTable.taxa(listOf("A", "B")), FactorAttribute("env", intArrayOf(0, 0), listOf("e1", "e2"))),
        )
        val result = filterPhenotype(withSpare, PhenotypeFilter(factor = "env", levels = listOf("e2")))
        assertEquals(0, result.numberOfObservations)
    }

    @Test
    fun unknownNamesAreErrorsThatListTheValidOnes() = runTest {
        suspend fun message(f: PhenotypeFilter): String =
            assertFailsWith<IllegalArgumentException> { filter(f) }.message.orEmpty()
        assertTrue("[A, B, C, D]" in message(PhenotypeFilter(taxa = listOf("Z"))))
        assertTrue("[e1, e2]" in message(PhenotypeFilter(factor = "env", levels = listOf("e3"))))
        assertTrue("Taxa, env, yield, height, q" in message(PhenotypeFilter(keepAttributes = listOf("nope"))))
        assertTrue("nope" in message(PhenotypeFilter(dropAttributes = listOf("nope"))))
        assertTrue("nope" in message(PhenotypeFilter(requireNonMissing = listOf("nope"))))
        assertTrue("nope" in message(PhenotypeFilter(attribute = "nope", min = 0.0)))
        assertTrue("nope" in message(PhenotypeFilter(toData = listOf("nope"))))
    }

    @Test
    fun rejectsMisusedAttributes() = runTest {
        suspend fun fails(f: PhenotypeFilter) = assertFailsWith<IllegalArgumentException> { filter(f) }
        fails(PhenotypeFilter(factor = "yield", levels = listOf("1")))
        fails(PhenotypeFilter(attribute = "env", min = 0.0))
        fails(PhenotypeFilter(toCovariate = listOf("env")))
        fails(PhenotypeFilter(toData = listOf("Taxa")))
        fails(PhenotypeFilter(dropAttributes = listOf("Taxa")))
    }

    @Test
    fun rejectsInconsistentFilters() {
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(factor = "env") }
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(levels = listOf("e1")) }
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(min = 1.0) }
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(attribute = "yield") }
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(attribute = "yield", min = 2.0, max = 1.0) }
        assertFailsWith<IllegalArgumentException> { PhenotypeFilter(toData = listOf("q"), toCovariate = listOf("q")) }
    }

    // The spec

    @Test
    fun specHasTassel5IdentityAndLegacyName() {
        assertEquals("filter-phenotype", filterPhenotypeOperation.id)
        assertEquals("FilterTraitsPlugin", filterPhenotypeOperation.legacyPlugin)
        assertEquals("excludeLast", filterPhenotypeOperation.parameterByLegacyName("excludeLastTrait")?.name)
        assertEquals(false, filterPhenotypeOperation.parameter("excludeLast")?.default)
        assertEquals(true, filterPhenotypeOperation.parameter("includeTaxa")?.default)
        assertEquals(true, filterPhenotypeOperation.parameter("includeLevels")?.default)
        assertEquals(false, filterPhenotypeOperation.parameter("requireAllNonMissing")?.default)
    }

    @Test
    fun specBindsTassel5TextAndRunsThroughTheEngine() = runTest {
        val engine = LocalEngine(OperationRegistry(filterPhenotypeOperation))
        val handle = engine.put(input)

        val noArguments = engine.run(RunRequest("filter-phenotype", handle))
        assertSame(input, engine.get(noArguments))

        val excluded = engine.run(RunRequest("filter-phenotype", handle, mapOf("excludeLast" to "true")))
        val result = engine.get(excluded) as PhenotypeTable
        assertEquals(listOf("Taxa", "env", "yield", "height"), result.attributes.map { it.name })

        val filtered = engine.run(
            RunRequest(
                "filter-phenotype",
                handle,
                mapOf(
                    "taxa" to "A,B",
                    "factor" to "env",
                    "levels" to "e1",
                    "attribute" to "height",
                    "max" to "30",
                    "keepTypes" to "data,covariate",
                    "toCovariate" to "yield",
                ),
            ),
        )
        val out = engine.get(filtered) as PhenotypeTable
        assertEquals(listOf("A", "B"), (0 until out.numberOfObservations).map { out.taxon(it) })
        assertEquals(listOf("Taxa", "yield", "height", "q"), out.attributes.map { it.name })
        assertEquals(AttributeType.covariate, out.attribute("yield").type)
        assertContentEquals(listOf("10.0", "30.0"), rows(out, "height"))
    }

    @Test
    fun specRejectsAnUnknownAttributeType() = runTest {
        assertFailsWith<IllegalArgumentException> {
            filterPhenotypeOperation.run(input, mapOf("keepTypes" to "trait"))
        }
    }
}
