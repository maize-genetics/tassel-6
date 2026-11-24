package net.maizegenetics.tassel.join

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGolden
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * [joinPhenotypes] and the `join-phenotypes` spec (ADR-0024). The fixtures in `testdata/phenotype/joins/` are
 * hand-written, each with an expected file whose row and column order is part of the check, because TASSEL 5 sorts
 * its joins and these do not. TASSEL 5's `-intersect` and `-union` are cross-checks only, at `rows-unordered`, on
 * the tutorial inputs where both semantics agree.
 *
 * @see net.maizegenetics.phenotype.PhenotypeBuilder.mergeTwoPhenotypes
 */
class JoinPhenotypesTest {
    private suspend fun load(name: String): PhenotypeTable =
        PhenotypeFiles.read(Buffer().apply { write(TestData.readBytes("phenotype/joins/$name.txt")) }).table

    private fun expected(name: String): String = TestData.readText("phenotype/joins/$name.expected.txt")

    private suspend fun join(
        inputs: List<String>,
        type: JoinType,
        multiple: JoinMultiple = JoinMultiple.error,
        collision: JoinCollision = JoinCollision.error,
    ) = joinPhenotypes(inputs.map { load(it) }, JoinOptions(type, multiple, collision)).toDelimitedString()

    private suspend fun assertJoin(
        expected: String,
        vararg inputs: String,
        type: JoinType,
        multiple: JoinMultiple = JoinMultiple.error,
        collision: JoinCollision = JoinCollision.error,
    ) = assertEquals(expected(expected), join(inputs.toList(), type, multiple, collision), expected)

    private suspend fun failure(
        vararg inputs: String,
        type: JoinType = JoinType.intersect,
        multiple: JoinMultiple = JoinMultiple.error,
        collision: JoinCollision = JoinCollision.error,
    ): String = assertFailsWith<IllegalArgumentException> {
        join(inputs.toList(), type, multiple, collision)
    }.message.orEmpty()

    // --- Fixtures -----------------------------------------------------------------------------

    @Test
    fun replicatesJoinAQMatrixOneToMany() = runTest {
        assertJoin("reps-q.intersect", "reps", "q", type = JoinType.intersect)
        assertJoin("reps-q.union", "reps", "q", type = JoinType.union)
    }

    @Test
    fun swappingSidesChangesRowAndColumnOrderAndWhereRightOnlyRowsGo() = runTest {
        assertJoin("q-reps.intersect", "q", "reps", type = JoinType.intersect)
        assertJoin("q-reps.union", "q", "reps", type = JoinType.union)
    }

    @Test
    fun aSharedFactorIsPartOfTheKeyAndReplicatesStayOneToMany() = runTest {
        assertJoin("env-a-b.intersect", "env-a", "env-b", type = JoinType.intersect)
        assertJoin("env-a-b.union", "env-a", "env-b", type = JoinType.union)
    }

    @Test
    fun everySharedFactorIsPartOfTheKey() = runTest {
        assertJoin("env-a-c.intersect", "env-a", "env-c", type = JoinType.intersect)
        assertJoin("env-a-c.union", "env-a", "env-c", type = JoinType.union)
    }

    @Test
    fun aKeyRepeatedOnBothSidesFailsNamingTheKeysUnlessCartesian() = runTest {
        val message = failure("reps", "mm-b")
        assertContains(message, "B73")
        assertTrue("Mo17" !in message, message)
        assertContains(message, "cartesian")
        assertJoin(
            "reps-mm.cartesian.intersect",
            "reps",
            "mm-b",
            type = JoinType.intersect,
            multiple = JoinMultiple.cartesian,
        )
        assertJoin("reps-mm.cartesian.union", "reps", "mm-b", type = JoinType.union, multiple = JoinMultiple.cartesian)
    }

    @Test
    fun ambiguousKeysNameTheFactorLevelsAndAreCutAfterTen() = runTest {
        val taxa = (1..12).map { "t$it" }
        val twice = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(taxa + taxa),
                FactorAttribute.of("env", List(24) { "E1" }),
            ),
        )
        val message = assertFailsWith<IllegalArgumentException> {
            joinPhenotypes(listOf(twice, twice))
        }.message.orEmpty()
        assertContains(message, "t1 (env=E1)")
        assertContains(message, "t10 (env=E1)")
        assertTrue("t11" !in message, message)
        assertContains(message, "and 2 more")
    }

    @Test
    fun collisionsAreAnErrorByDefaultAndSuffixedOrCoalescedOnRequest() = runTest {
        assertContains(failure("col-a", "col-b"), "'yield'")
        assertJoin("col-a-b.suffix", "col-a", "col-b", type = JoinType.intersect, collision = JoinCollision.suffix)
        assertJoin("col-a-b.coalesce", "col-a", "col-b", type = JoinType.intersect, collision = JoinCollision.coalesce)
    }

    @Test
    fun coalesceFailsNamingTheKeyAndColumnWhenBothValuesDiffer() = runTest {
        val message = failure("col-a", "col-conflict", collision = JoinCollision.coalesce)
        assertContains(message, "'yield'")
        assertContains(message, "B73")
    }

    @Test
    fun aNameWithTwoTypesIsACollisionThatOnlySuffixCanResolve() = runTest {
        assertContains(failure("type-a", "type-b"), "'env'")
        assertContains(failure("type-a", "type-b", collision = JoinCollision.coalesce), "different tables")
        assertJoin("type-a-b.suffix", "type-a", "type-b", type = JoinType.intersect, collision = JoinCollision.suffix)
    }

    @Test
    fun aSuffixedNameThatIsAlreadyTakenFails() = runTest {
        val left = PhenotypeTable(
            listOf(PhenotypeTable.taxa(listOf("a")), num("yield", 1f), num("yield.1", 2f)),
        )
        val right = PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("a")), num("yield", 3f)))
        val message = assertFailsWith<IllegalArgumentException> {
            joinPhenotypes(listOf(left, right), JoinOptions(collision = JoinCollision.suffix))
        }.message.orEmpty()
        assertContains(message, "'yield.1'")
    }

    @Test
    fun aMissingKeyValueMatchesNothingNotEvenAnotherMissingOne() = runTest {
        assertJoin("miss.intersect", "miss-a", "miss-b", type = JoinType.intersect)
        assertJoin("miss.union", "miss-a", "miss-b", type = JoinType.union)
    }

    @Test
    fun threeInputsAreALeftFoldAndSuffixNumbersByOriginalPosition() = runTest {
        val inputs = arrayOf("fold-a", "fold-b", "fold-c")
        assertJoin("fold.suffix.intersect", *inputs, type = JoinType.intersect, collision = JoinCollision.suffix)
        assertJoin("fold.suffix.union", *inputs, type = JoinType.union, collision = JoinCollision.suffix)
        assertContains(failure(*inputs), "'yield'")
    }

    // --- Rules not tied to a fixture ----------------------------------------------------------

    private fun num(name: String, vararg values: Float, type: AttributeType = AttributeType.data) =
        NumericAttribute(name, type, values)

    @Test
    fun factorsMatchByLabelNotByLevelCode() = runTest {
        val left = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(listOf("a", "a")),
                FactorAttribute.of("env", listOf("E1", "E2")),
                num("x", 1f, 2f),
            ),
        )
        // Only E2 here, so its code is 0 where the left's E2 code is 1.
        val right = PhenotypeTable(
            listOf(PhenotypeTable.taxa(listOf("a")), FactorAttribute.of("env", listOf("E2")), num("y", 9f)),
        )
        val joined = joinPhenotypes(listOf(left, right))
        assertEquals(listOf("a"), joined.taxaAttribute.names)
        assertEquals(2f, (joined.attribute("x") as NumericAttribute).floatValue(0))
        assertEquals(9f, (joined.attribute("y") as NumericAttribute).floatValue(0))
    }

    @Test
    fun taxaMatchExactlyAndTheTaxaColumnNameIsIgnored() = runTest {
        val left =
            PhenotypeTable(
                listOf(
                    net.maizegenetics.tassel.phenotype.TaxaAttribute("Line", listOf("B73", "b73")),
                    num("x", 1f, 2f),
                ),
            )
        val right =
            PhenotypeTable(
                listOf(net.maizegenetics.tassel.phenotype.TaxaAttribute("Genotype", listOf("B73")), num("y", 5f)),
            )
        val joined = joinPhenotypes(listOf(left, right), JoinOptions(JoinType.union))
        assertEquals("Line", joined.taxaAttribute.name)
        assertEquals(listOf("B73", "b73"), joined.taxaAttribute.names)
        assertTrue((joined.attribute("y") as NumericAttribute).isMissing(1))
    }

    @Test
    fun covariatesAreNeverKeys() = runTest {
        val left =
            PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("a")), num("Q1", 0.5f, type = AttributeType.covariate)))
        val right =
            PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("a")), num("Q1", 0.7f, type = AttributeType.covariate)))
        assertContains(
            assertFailsWith<IllegalArgumentException> { joinPhenotypes(listOf(left, right)) }.message.orEmpty(),
            "'Q1'",
        )
        val joined = joinPhenotypes(listOf(left, right), JoinOptions(collision = JoinCollision.suffix))
        assertEquals(listOf("Taxa", "Q1.1", "Q1.2"), joined.attributes.map { it.name })
    }

    @Test
    fun aTaxaColumnNameThatIsAnotherKindOfColumnElsewhereFails() = runTest {
        val left = PhenotypeTable(listOf(PhenotypeTable.taxa(listOf("a")), num("x", 1f)))
        val right = PhenotypeTable(
            listOf(net.maizegenetics.tassel.phenotype.TaxaAttribute("Line", listOf("a")), num("Taxa", 2f)),
        )
        assertFailsWith<IllegalArgumentException> { joinPhenotypes(listOf(left, right)) }
    }

    @Test
    fun aJoinNeedsAtLeastTwoTables() = runTest {
        assertContains(
            assertFailsWith<IllegalArgumentException> { joinPhenotypes(listOf(load("reps"))) }.message.orEmpty(),
            "at least two",
        )
    }

    @Test
    fun theJoinReportsProgressPerFoldStep() = runTest {
        val steps = mutableListOf<Long>()
        joinPhenotypes(
            listOf("fold-a", "fold-b", "fold-c").map { load(it) },
            JoinOptions(collision = JoinCollision.suffix),
            net.maizegenetics.tassel.ops.ProgressReporter { steps += it.completed },
        )
        assertEquals(listOf(1L, 2L), steps)
    }

    // --- DataSet and spec ---------------------------------------------------------------------

    @Test
    fun aDataSetKeepsOrderFlattensAndSortsItemsByType() = runTest {
        val a = load("reps")
        val b = load("q")
        val set = DataSet(listOf(a, DataSet.of(b, a)))
        assertEquals(listOf(a, b, a), set.items)
        assertEquals(listOf(a, b, a), set.phenotypes())
        assertEquals(emptyList(), set.genotypes())
        assertEquals(3, set.size)
    }

    @Test
    fun specDeclaresTheKernelsParametersAndRunsThroughTheEngineWithADataSetInput() = runTest {
        assertEquals("join-phenotypes", joinPhenotypesOperation.id)
        assertEquals("IntersectionAlignmentPlugin", joinPhenotypesOperation.legacyPlugin)
        assertEquals("DataSet", joinPhenotypesOperation.descriptor.inputType)
        assertEquals(JoinType.intersect, joinPhenotypesOperation.parameter("type")?.default)
        assertEquals(JoinMultiple.error, joinPhenotypesOperation.parameter("multiple")?.default)
        assertEquals(JoinCollision.error, joinPhenotypesOperation.parameter("collision")?.default)

        val engine = LocalEngine(OperationRegistry(joinPhenotypesOperation))
        val input = engine.put(DataSet.of(load("reps"), load("q")))
        val union = engine.get(
            engine.run(RunRequest("join-phenotypes", input, mapOf("type" to "union"))),
        ) as PhenotypeTable
        assertEquals(expected("reps-q.union"), union.toDelimitedString())
        val inner = engine.get(engine.run(RunRequest("join-phenotypes", input))) as PhenotypeTable
        assertEquals(expected("reps-q.intersect"), inner.toDelimitedString())
    }

    @Test
    fun specRejectsAnInputThatIsNotAtLeastTwoPhenotypes() = runTest {
        val engine = LocalEngine(OperationRegistry(joinPhenotypesOperation))
        val one = engine.put(DataSet.of(load("reps")))
        val message = assertFailsWith<IllegalArgumentException> {
            engine.run(RunRequest("join-phenotypes", one))
        }.message.orEmpty()
        assertContains(message, "at least two phenotype tables")
        assertFailsWith<IllegalArgumentException> {
            engine.run(RunRequest("join-phenotypes", engine.put(load("reps"))))
        }
    }

    // --- TASSEL 5 cross-check -----------------------------------------------------------------

    @Test
    fun matchesTassel5sIntersectAndUnionOfTheTutorialTraitsAndPopulationStructure() = runTest {
        suspend fun tutorial(name: String) =
            PhenotypeFiles.read(Buffer().apply { write(TestData.readBytes("tutorial/$name.txt")) }).table
        val traits = tutorial("mdp_traits")
        val structure = tutorial("mdp_population_structure")
        assertMatchesGolden(
            "golden/phenotype-join/mdp_traits.q.intersect.txt",
            joinPhenotypes(listOf(traits, structure), JoinOptions(JoinType.intersect)).toDelimitedString(),
        )
        assertMatchesGolden(
            "golden/phenotype-join/mdp_traits.q.union.txt",
            joinPhenotypes(listOf(traits, structure), JoinOptions(JoinType.union)).toDelimitedString(),
        )
    }
}
