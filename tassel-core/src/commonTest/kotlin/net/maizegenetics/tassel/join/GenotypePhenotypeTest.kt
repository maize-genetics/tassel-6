package net.maizegenetics.tassel.join

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * [alignGenotypePhenotype] and the `join-genotype-phenotype` spec (ADR-0024), on hand-written tables. The genotype
 * taxa are deliberately not in phenotype order, as in TASSEL 5's `GenotypePhenotypeBuilderTest`.
 *
 * @see net.maizegenetics.phenotype.GenotypePhenotypeBuilder
 */
class GenotypePhenotypeTest {
    private val positions = PositionList((1..4).map { Position(Chromosome.of(1), it * 100) })

    /** Genotype taxa c, a, d, b; the call of taxon `t` at site `s` is `10 t + s`. */
    private class Counting(private val inner: GenotypeTable) : GenotypeTable {
        var cellsRead = 0
        var wholeSiteOrTaxonReads = 0
        override val taxa get() = inner.taxa
        override val positions get() = inner.positions

        override fun genotype(taxon: Int, site: Int): Byte {
            cellsRead++
            return inner.genotype(taxon, site)
        }

        override fun siteGenotypes(site: Int): ByteArray {
            wholeSiteOrTaxonReads++
            return super.siteGenotypes(site)
        }

        override fun taxonGenotypes(taxon: Int): ByteArray {
            wholeSiteOrTaxonReads++
            return super.taxonGenotypes(taxon)
        }
    }

    private fun genotypes(vararg names: String) =
        InMemoryGenotypeTable.of(TaxaList.of(*names), positions) { t, s -> (10 * t + s).toByte() }

    private val genotypeTable = genotypes("c", "a", "d", "b")

    /** Taxa a, b, b (a replicate), x (no genotype), c. */
    private val phenotype = PhenotypeTable(
        listOf(
            PhenotypeTable.taxa(listOf("a", "b", "b", "x", "c")),
            FactorAttribute.of("rep", listOf("1", "1", "2", "1", null)),
            NumericAttribute("yield", AttributeType.data, floatArrayOf(1f, 2f, 3f, 4f, Float.NaN)),
        ),
        "pheno",
    )

    private fun taxaOf(table: PhenotypeTable) = (0 until table.numberOfObservations).map { table.taxon(it) }

    @Test
    fun intersectDropsObservationsWithoutGenotypesAndMapsRowsInObservationOrder() {
        val joined = alignGenotypePhenotype(genotypeTable, phenotype, JoinType.intersect)
        assertEquals(listOf("a", "b", "b", "c"), taxaOf(joined.phenotype))
        assertEquals(4, joined.numberOfObservations)
        assertContentEquals(intArrayOf(1, 3, 3, 0), joined.genotypeRows())
        assertContentEquals(intArrayOf(1, 3, 0), joined.genotypeTaxonIndices())
        assertContentEquals(intArrayOf(0, 1, 2, 3), joined.observationsWithGenotype())
        // Columns and cells follow their rows.
        assertEquals(
            listOf("1", "1", "2", null),
            (0 until 4).map {
                (joined.phenotype.attribute("rep") as FactorAttribute).level(it)
            },
        )
        assertEquals(
            listOf(1f, 2f, 3f),
            (0 until 3).map {
                (joined.phenotype.attribute("yield") as NumericAttribute).floatValue(it)
            },
        )
        assertTrue((joined.phenotype.attribute("yield") as NumericAttribute).isMissing(3))
        assertEquals("pheno", joined.phenotype.name)
    }

    @Test
    fun unionKeepsUngenotypedObservationsAndAppendsGenotypeOnlyTaxaInGenotypeOrder() {
        val bigger = genotypes("c", "e", "a", "d", "b") // d and e are genotype-only
        val joined = alignGenotypePhenotype(bigger, phenotype, JoinType.union)
        assertEquals(listOf("a", "b", "b", "x", "c", "e", "d"), taxaOf(joined.phenotype))
        assertContentEquals(intArrayOf(2, 4, 4, -1, 0, 1, 3), joined.genotypeRows())
        assertEquals(7, joined.numberOfObservations)
        assertContentEquals(intArrayOf(2, 4, 0, 1, 3), joined.genotypeTaxonIndices())
        assertContentEquals(intArrayOf(0, 1, 2, 4, 5, 6), joined.observationsWithGenotype())
        val rep = joined.phenotype.attribute("rep") as FactorAttribute
        val yield = joined.phenotype.attribute("yield") as NumericAttribute
        assertEquals(listOf(null, null), listOf(5, 6).map { rep.level(it) })
        assertTrue(yield.isMissing(5) && yield.isMissing(6))
        assertEquals(4f, yield.floatValue(3)) // the ungenotyped observation keeps its values
        assertEquals(listOf("rep", "yield"), joined.phenotype.attributes.drop(1).map { it.name })
    }

    @Test
    fun unionWithNothingToAppendReturnsThePhenotypeUnchanged() {
        val joined = alignGenotypePhenotype(genotypes("a", "b", "c"), phenotype, JoinType.union)
        assertEquals(phenotype, joined.phenotype)
        assertContentEquals(intArrayOf(0, 1, 1, -1, 2), joined.genotypeRows())
    }

    @Test
    fun taxaMatchExactlyAndCaseSensitively() {
        val joined = alignGenotypePhenotype(genotypes("A", "a"), phenotype, JoinType.intersect)
        assertEquals(listOf("a"), taxaOf(joined.phenotype))
        assertContentEquals(intArrayOf(1), joined.genotypeRows())
    }

    @Test
    fun disjointTaxaGiveNoObservationsUnderIntersectAndAllRowsUnderUnion() {
        val other = genotypes("p", "q")
        val none = alignGenotypePhenotype(other, phenotype, JoinType.intersect)
        assertEquals(0, none.numberOfObservations)
        assertEquals(0, none.genotypeTaxonIndices().size)
        val all = alignGenotypePhenotype(other, phenotype, JoinType.union)
        assertEquals(listOf("a", "b", "b", "x", "c", "p", "q"), taxaOf(all.phenotype))
        assertContentEquals(intArrayOf(-1, -1, -1, -1, -1, 0, 1), all.genotypeRows())
    }

    @Test
    fun siteBlocksOnTheMappedIndicesReturnEachObservationsGenotype() {
        val joined = alignGenotypePhenotype(genotypeTable, phenotype, JoinType.intersect)
        val indices = joined.genotypeTaxonIndices() // a, b, c = rows 1, 3, 0
        val blocks = joined.genotypes.siteBlocks(blockSize = 3, taxa = indices).toList()
        assertEquals(listOf(0..2, 3..3), blocks.map { it.sites })
        // Observation i reads block taxon position of its genotype row, whatever the table's own order.
        for (observation in 0 until joined.numberOfObservations) {
            val position = indices.indexOf(joined.genotypeRow(observation))
            for (block in blocks) {
                for (site in block.sites) {
                    val expected = (
                        10 * genotypeTable.taxa.indexOf(
                            joined.phenotype.taxon(observation),
                        ) + site
                        ).toByte()
                    assertEquals(expected, block.genotype(position, site), "observation $observation site $site")
                }
            }
        }
    }

    @Test
    fun unionBlocksCoverGenotypeOnlyTaxaToo() {
        val joined = alignGenotypePhenotype(genotypeTable, phenotype, JoinType.union)
        val block = joined.genotypes.siteBlocks(taxa = joined.genotypeTaxonIndices()).single()
        // Taxa a, b, c (observed), then d (appended): table rows 1, 3, 0, 2.
        assertEquals(listOf<Byte>(10, 30, 0, 20), (0 until 4).map { block.genotype(it, 0) })
    }

    @Test
    fun alignmentAndMappedBlocksNeverMaterializeTheTaxaBySiteMatrix() {
        val stub = Counting(genotypeTable)
        val joined = alignGenotypePhenotype(stub, phenotype, JoinType.union)
        assertEquals(0, stub.cellsRead, "aligning reads no genotypes")
        assertEquals(0, stub.wholeSiteOrTaxonReads)
        joined.genotypeTaxonIndices()
        joined.observationsWithGenotype()
        assertEquals(0, stub.cellsRead)

        val indices = joined.genotypeTaxonIndices()
        joined.genotypes.siteBlocks(taxa = indices).forEach { block -> block.siteGenotypes(block.sites.first) }
        assertEquals(indices.size * stub.numberOfSites, stub.cellsRead, "one read per mapped taxon and site, once")
        assertEquals(0, stub.wholeSiteOrTaxonReads, "no whole-site or whole-taxon vector")
        assertTrue(indices.size < stub.numberOfTaxa || indices.size == 4)

        // Intersect maps fewer taxa than the table holds, and only those are read.
        val inner = Counting(genotypes("c", "a", "d", "b", "e", "f"))
        val aligned = alignGenotypePhenotype(inner, phenotype, JoinType.intersect)
        aligned.genotypes.siteBlocks(taxa = aligned.genotypeTaxonIndices()).count()
        assertEquals(3 * inner.numberOfSites, inner.cellsRead)
    }

    @Test
    fun theValueRejectsMappingsThatDisagreeWithThePhenotype() {
        assertFailsWith<IllegalArgumentException> { GenotypePhenotype(genotypeTable, phenotype, IntArray(2)) }
        assertFailsWith<IllegalArgumentException> {
            GenotypePhenotype(genotypeTable, phenotype, intArrayOf(1, 3, 3, -1, 9))
        }
        val message = assertFailsWith<IllegalArgumentException> {
            GenotypePhenotype(genotypeTable, phenotype, intArrayOf(0, 3, 3, -1, 0)) // observation 0 is a, row 0 is c
        }.message.orEmpty()
        assertContains(message, "'a'")
        // The caller's array is copied.
        val rows = intArrayOf(1, 3, 3, -1, 0)
        val value = GenotypePhenotype(genotypeTable, phenotype, rows)
        rows[0] = 0
        assertEquals(1, value.genotypeRow(0))
    }

    // --- Spec -----------------------------------------------------------------------------------

    @Test
    fun specDeclaresTypeAndRunsThroughTheEngineWithTheGenotypeAndPhenotypeInEitherOrder() = runTest {
        assertEquals("join-genotype-phenotype", joinGenotypePhenotypeOperation.id)
        assertEquals("DataSet", joinGenotypePhenotypeOperation.descriptor.inputType)
        assertEquals(JoinType.intersect, joinGenotypePhenotypeOperation.parameter("type")?.default)

        val engine = LocalEngine(OperationRegistry(joinGenotypePhenotypeOperation))
        for (items in listOf(listOf(genotypeTable, phenotype), listOf(phenotype, genotypeTable))) {
            val input = engine.put(DataSet(items))
            val inner = engine.get(engine.run(RunRequest("join-genotype-phenotype", input))) as GenotypePhenotype
            assertContentEquals(intArrayOf(1, 3, 3, 0), inner.genotypeRows())
            val union = engine.get(
                engine.run(RunRequest("join-genotype-phenotype", input, mapOf("type" to "union"))),
            ) as GenotypePhenotype
            assertContentEquals(intArrayOf(1, 3, 3, -1, 0, 2), union.genotypeRows())
        }
    }

    @Test
    fun severalPhenotypesAreJoinedFirstWithTheSameType() = runTest {
        val extra = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(listOf("a", "d")),
                NumericAttribute("height", AttributeType.data, floatArrayOf(7f, 8f)),
            ),
        )
        val engine = LocalEngine(OperationRegistry(joinGenotypePhenotypeOperation))
        val input = engine.put(DataSet.of(phenotype, genotypeTable, extra))
        val inner = engine.get(engine.run(RunRequest("join-genotype-phenotype", input))) as GenotypePhenotype
        assertEquals(listOf("a"), taxaOf(inner.phenotype)) // phenotype joined to extra by intersect: only a is in both
        assertEquals(7f, (inner.phenotype.attribute("height") as NumericAttribute).floatValue(0))
        val union = engine.get(
            engine.run(RunRequest("join-genotype-phenotype", input, mapOf("type" to "union"))),
        ) as GenotypePhenotype
        // Phenotype rows a, b, b, x, c, then d from the second table, then no genotype-only taxa remain.
        assertEquals(listOf("a", "b", "b", "x", "c", "d"), taxaOf(union.phenotype))
        assertContentEquals(intArrayOf(1, 3, 3, -1, 0, 2), union.genotypeRows())
    }

    @Test
    fun specRejectsInputsThatAreNotOneGenotypeAndAtLeastOnePhenotype() = runTest {
        val engine = LocalEngine(OperationRegistry(joinGenotypePhenotypeOperation))
        suspend fun message(vararg items: Any) = assertFailsWith<IllegalArgumentException> {
            engine.run(RunRequest("join-genotype-phenotype", engine.put(DataSet(items.toList()))))
        }.message.orEmpty()
        assertContains(message(phenotype), "exactly one genotype table")
        assertContains(message(genotypeTable), "at least one phenotype table")
        assertContains(message(genotypeTable, genotypes("a"), phenotype), "2 genotype tables")
        assertContains(message(genotypeTable, phenotype, "stray"), "nothing else")
    }
}
