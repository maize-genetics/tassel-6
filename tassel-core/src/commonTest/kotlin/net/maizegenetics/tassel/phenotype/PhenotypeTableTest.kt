package net.maizegenetics.tassel.phenotype

import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertNotSame
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * The phenotype model (ADR-0023) and the `<Phenotype>` text it writes, against the shape of
 * `testdata/golden/pca/pcs.txt`.
 *
 * @see net.maizegenetics.phenotype.PhenotypeUtils.write
 */
class PhenotypeTableTest {
    @Test
    fun writesTheThreeHeaderLinesAndOneRowPerObservation() {
        assertEquals(
            "<Phenotype>\n" +
                "taxa\tcovariate\tcovariate\n" +
                "Taxa\tPC1\tPC2\n" +
                "one\t0.5\t-1.25\n" +
                "two\t2.0\t3.5\n",
            pcs.toDelimitedString(),
        )
    }

    @Test
    fun writesFloatsInJavasSpellingOnEveryTarget() {
        val tiny =
            PhenotypeTable(
                listOf(PhenotypeTable.taxa(listOf("one")), NumericAttribute("PC1", floatArrayOf(3.1861675e-4f))),
            )
        assertTrue(tiny.toDelimitedString().endsWith("one\t3.1861675E-4\n"), tiny.toDelimitedString())
    }

    @Test
    fun writesMissingNumbersAsNaNAndMissingLevelsAsAQuestionMark() {
        assertEquals(
            "<Phenotype>\n" +
                "taxa\tfactor\tdata\tcovariate\n" +
                "Taxon\tenv\tEarHT\tQ1\n" +
                "B73\tIA\t64.75\t0.014\n" +
                "B73\tNY\tNaN\t0.014\n" +
                "Mo17\t?\t81.0\tNaN\n",
            trial.toDelimitedString(),
        )
    }

    /** The format the PCA golden is written in, read back off the golden itself. */
    @Test
    fun matchesTheShapeOfTheGoldenPhenotype() {
        val lines = normalizeLineEndings(goldenText("golden/pca/pcs.txt")).lines()
        assertEquals(PhenotypeTable.HEADER, lines[0])
        assertEquals(listOf(AttributeType.taxa.name) + List(5) { "covariate" }, lines[1].split("\t"))
        assertEquals(listOf(PhenotypeTable.TAXA) + (1..5).map { "PC$it" }, lines[2].split("\t"))
    }

    @Test
    fun knowsItsTaxaAndReplicates() {
        assertEquals(listOf("B73", "Mo17"), trial.taxa.names)
        assertTrue(trial.hasReplicatedTaxa)
        assertFalse(pcs.hasReplicatedTaxa)
        assertEquals("Mo17", trial.taxon(2))
        assertEquals("Taxon", trial.taxaAttribute.name)
        assertEquals(3, trial.numberOfObservations)
        assertEquals(4, trial.numberOfAttributes)
    }

    @Test
    fun findsAttributesByNameAndByType() {
        assertEquals(2, trial.indexOf("EarHT"))
        assertEquals(-1, trial.indexOf("earht"))
        assertEquals(AttributeType.covariate, trial.attribute("Q1").type)
        assertEquals(listOf("EarHT"), trial.attributes(AttributeType.data).map { it.name })
        val failure = assertFailsWith<IllegalArgumentException> { trial.attribute("Missing") }
        assertEquals("No attribute named 'Missing'; the attributes are [Taxon, env, EarHT, Q1]", failure.message)
    }

    @Test
    fun selectingRowsMakesANewTableAndLeavesTheOriginalAlone() {
        val selected = trial.selectRows(intArrayOf(2, 0, 0))
        assertEquals(listOf("Mo17", "B73", "B73"), selected.taxaAttribute.names)
        assertEquals(listOf("IA"), (selected.attribute("env") as FactorAttribute).levels)
        assertNull((selected.attribute("env") as FactorAttribute).level(0))
        assertEquals(3, trial.numberOfObservations)
        assertEquals(listOf("IA", "NY"), (trial.attribute("env") as FactorAttribute).levels)
        assertFailsWith<IllegalArgumentException> { trial.selectRows(intArrayOf(3)) }
    }

    @Test
    fun selectingAttributesKeepsTheTaxaColumnFirst() {
        assertEquals(
            listOf("Taxon", "Q1", "EarHT"),
            trial.selectAttributes(listOf("Q1", "EarHT")).attributes.map {
                it.name
            },
        )
        assertEquals(
            listOf("EarHT", "Taxon"),
            trial.selectAttributes(listOf("EarHT", "Taxon")).attributes.map {
                it.name
            },
        )
        assertFailsWith<IllegalArgumentException> { trial.selectAttributes(listOf("Q1", "Q1")) }
        assertFailsWith<IllegalArgumentException> { trial.selectAttributes(listOf("nope")) }
    }

    @Test
    fun rejectsTablesThatBreakTheModel() {
        val taxa = PhenotypeTable.taxa(listOf("one", "two"))
        val message = { block: () -> Unit -> assertFailsWith<IllegalArgumentException> { block() }.message }
        assertEquals("A phenotype needs at least one attribute", message { PhenotypeTable(emptyList()) })
        assertEquals(
            "Every attribute needs the same number of observations, found [2, 1]",
            message { PhenotypeTable(listOf(taxa, NumericAttribute("PC1", floatArrayOf(1f)))) },
        )
        assertEquals(
            "A phenotype needs exactly one taxa attribute, found 0",
            message { PhenotypeTable(listOf(NumericAttribute("PC1", floatArrayOf(1f)))) },
        )
        assertEquals(
            "A phenotype needs exactly one taxa attribute, found 2",
            message { PhenotypeTable(listOf(taxa, TaxaAttribute("Other", listOf("a", "b")))) },
        )
        assertEquals(
            "Attribute name 'x' appears more than once",
            message {
                PhenotypeTable(
                    listOf(
                        taxa,
                        NumericAttribute("x", floatArrayOf(1f, 2f)),
                        FactorAttribute.of("x", listOf("a", "b")),
                    ),
                )
            },
        )
    }

    @Test
    fun numericAttributesHoldFloatsWithNaNAsTheOnlyMissingValue() {
        val values =
            NumericAttribute(
                "t",
                AttributeType.data,
                floatArrayOf(1.5f, Float.NaN, Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY),
            )
        assertEquals(1.5, values.value(0))
        assertEquals(listOf(false, true, true, true), (0 until 4).map { values.isMissing(it) })
        assertEquals(3, values.missingCount)
        assertEquals(listOf("1.5", "NaN", "NaN", "NaN"), (0 until 4).map { values.format(it) })
        assertFailsWith<IllegalArgumentException> { NumericAttribute("t", AttributeType.factor, floatArrayOf(1f)) }
        assertFailsWith<IllegalArgumentException> { NumericAttribute("t", AttributeType.taxa, floatArrayOf(1f)) }
    }

    @Test
    fun numericAttributesCopyTheArrayTheyAreGiven() {
        val source = floatArrayOf(1f, 2f)
        val values = NumericAttribute("t", AttributeType.data, source)
        source[0] = 9f
        assertEquals(1f, values.floatValue(0))
        val copy = values.toFloatArray()
        copy[1] = 9f
        assertEquals(2f, values.floatValue(1))
        assertNotSame(copy, values.toFloatArray())
    }

    @Test
    fun factorLevelsAreSortedAsStringsAndMissingIsCodeMinusOne() {
        val rep = FactorAttribute.of("rep", listOf("2", "10", null, "2", "1"))
        assertEquals(listOf("1", "10", "2"), rep.levels)
        assertContentEquals(intArrayOf(2, 1, -1, 2, 0), IntArray(5) { rep.code(it) })
        assertEquals(listOf("2", "10", "?", "2", "1"), (0 until 5).map { rep.format(it) })
        assertEquals(1, rep.missingCount)
    }

    @Test
    fun factorAttributesRejectBadLevelsAndCodes() {
        val message = { block: () -> Unit -> assertFailsWith<IllegalArgumentException> { block() }.message }
        assertEquals(
            "Factor 'f' levels must be unique and sorted; 'b' comes before 'a'",
            message { FactorAttribute("f", intArrayOf(0), listOf("b", "a")) },
        )
        assertEquals(
            "Factor 'f' cannot have the level '?', which means missing",
            message { FactorAttribute("f", intArrayOf(0), listOf("?")) },
        )
        assertEquals(
            "Factor 'f' has code 1 at observation 0, outside -1 until 1",
            message { FactorAttribute("f", intArrayOf(1), listOf("a")) },
        )
        assertEquals("Factor 'f' has an empty level", message { FactorAttribute.of("f", listOf("")) })
    }

    @Test
    fun taxaAttributesRejectEmptyNames() {
        assertFailsWith<IllegalArgumentException> { TaxaAttribute("Taxa", listOf("a", "")) }
        assertFailsWith<IllegalArgumentException> { TaxaAttribute("", listOf("a")) }
    }

    @Test
    fun refusesToWriteWhatWouldNotReadBackAsItself() {
        val message = { table: PhenotypeTable ->
            assertFailsWith<IllegalArgumentException> { table.toDelimitedString() }.message
        }
        val taxa = PhenotypeTable.taxa(listOf("a"))
        assertEquals(
            "Attribute name 'Ear HT' contains whitespace, which the file format splits on",
            message(PhenotypeTable(listOf(taxa, NumericAttribute("Ear HT", floatArrayOf(1f))))),
        )
        assertEquals(
            "Taxon 'B 73' contains whitespace, which the file format splits on",
            message(PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("B 73"))))),
        )
        assertEquals(
            "Level 'NA' of factor 'env' would be read back as missing",
            message(PhenotypeTable(listOf(taxa, FactorAttribute.of("env", listOf("NA"))))),
        )
    }

    @Test
    fun tablesAreEqualWhenTheirAttributesAndNamesAre() {
        assertEquals(trial, trial.selectRows(intArrayOf(0, 1, 2)))
        assertEquals(trial.hashCode(), trial.selectRows(intArrayOf(0, 1, 2)).hashCode())
        assertFalse(trial == PhenotypeTable(trial.attributes, "Other"))
        assertFalse(trial == trial.selectRows(intArrayOf(0, 1)))
    }

    private val pcs = PhenotypeTable(
        listOf(
            PhenotypeTable.taxa(listOf("one", "two")),
            NumericAttribute("PC1", floatArrayOf(0.5f, 2.0f)),
            NumericAttribute("PC2", floatArrayOf(-1.25f, 3.5f)),
        ),
    )

    private val trial = PhenotypeTable(
        listOf(
            TaxaAttribute("Taxon", listOf("B73", "B73", "Mo17")),
            FactorAttribute.of("env", listOf("IA", "NY", null)),
            NumericAttribute("EarHT", AttributeType.data, floatArrayOf(64.75f, Float.NaN, 81f)),
            NumericAttribute("Q1", AttributeType.covariate, floatArrayOf(0.014f, 0.014f, Float.NaN)),
        ),
    )
}
