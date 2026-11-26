package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.join.GenotypePhenotype
import net.maizegenetics.tassel.join.JoinType
import net.maizegenetics.tassel.join.alignGenotypePhenotype
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * [GenotypeMarkerSource] on hand-built tables: heterozygotes, partly missing calls, multiallelic sites, indels,
 * monomorphic sites, filters, and the one-read-per-site contract. Expected values are worked out by hand from TASSEL 5's
 * rules (`DiscreteSitesFELM`, `AbstractFixedEffectLM.applySiteFilters`, `ModelEffectUtils`).
 */
class GenotypeMarkerSourceTest {
    /** A copy of what [MarkerColumns] held, since its buffers are reused. */
    private data class Encoded(
        val skip: SkipReason?,
        val observations: List<Int>,
        val classes: List<String>,
        val classAlleles: List<String>,
        val classCounts: List<Int>,
        val classIndex: List<Int>,
        val columns: List<List<Double>>,
        val columnLabels: List<String>,
        val columnAlleles: List<Byte>,
        val comparison: Byte,
        val minorObs: Int,
    )

    private fun snapshot(columns: MarkerColumns) = Encoded(
        skip = columns.skip,
        observations = columns.observation.take(columns.size),
        classes = columns.classLabels,
        classAlleles = columns.classAlleleLabels,
        classCounts = columns.classCounts.toList(),
        classIndex = columns.classIndex.take(if (columns.numberOfClasses > 0) columns.size else 0),
        columns = (0 until columns.columnCount).map { columns.column(it).take(columns.size) },
        columnLabels = columns.columnLabels,
        columnAlleles = columns.columnAlleles.toList(),
        comparison = columns.comparisonAllele,
        minorObs = columns.minorObs,
    )

    private fun genotypes(sites: List<String>): InMemoryGenotypeTable {
        val calls = sites.map { site -> site.split(" ").map { Nucleotide.genotype(it) } }
        val taxa = calls.first().size
        require(calls.all { it.size == taxa }) { "Every site needs $taxa calls" }
        return InMemoryGenotypeTable.of(
            TaxaList.of((0 until taxa).map { "t$it" }),
            PositionList(sites.indices.map { Position(Chromosome.of(1), (it + 1) * 100) }),
        ) { taxon, site -> calls[site][taxon] }
    }

    private fun alignedWith(table: GenotypeTable, type: JoinType = JoinType.intersect): GenotypePhenotype {
        val names = table.taxa.names
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(names),
                NumericAttribute("y", AttributeType.data, FloatArray(names.size) { it.toFloat() }),
            ),
        )
        return alignGenotypePhenotype(table, phenotype, type)
    }

    /** Encodes marker 0 of a one-site table for every observation. */
    private fun encode(site: String, options: MarkerOptions = MarkerOptions()): Encoded =
        encodeSites(listOf(site), options)[0]

    private fun encodeSites(
        sites: List<String>,
        options: MarkerOptions = MarkerOptions(),
        observations: (Int) -> IntArray = { IntArray(it) { i -> i } },
        blockSize: Int = MarkerSource.DEFAULT_MARKER_BLOCK_SIZE,
    ): List<Encoded> {
        val data = alignedWith(genotypes(sites))
        val source = GenotypeMarkerSource.of(data)
        val columns = MarkerColumns(data.numberOfObservations)
        val results = ArrayList<Encoded>()
        source.scan(listOf(observations(data.numberOfObservations)), options, blockSize).forEach { block ->
            for (marker in block.markers) {
                val encoded = block.load(marker, 0, columns)
                assertEquals(encoded, columns.skip == null)
                results.add(snapshot(columns))
            }
        }
        return results
    }

    private val a = Allele.A
    private val c = Allele.C
    private val g = Allele.G

    @Test
    fun genotypeClassesMergePhasedHeterozygotesAndDropMissingCalls() {
        // Classes A (3), C (1), and M = A/C, C/A, A/C (3); N/N is missing.
        val result = encode("AA AC CA AA CC NN AC AA")
        assertNull(result.skip)
        assertEquals(listOf("A", "C", "M"), result.classes)
        assertEquals(listOf("A:A", "C:C", "A:C"), result.classAlleles)
        assertEquals(listOf(3, 1, 3), result.classCounts)
        assertEquals(listOf(0, 1, 2, 3, 4, 6, 7), result.observations)
        assertEquals(listOf(0, 2, 2, 0, 1, 2, 0), result.classIndex)
        // minorObs is the second largest class count: 3, 3, 1.
        assertEquals(3, result.minorObs)
    }

    @Test
    fun genotypeEncodingGivesTassel5AddFInputs() {
        // Allele copies over all genotyped taxa: A 9, C 5. One column for A, the first allele; C is the reference.
        val result = encode("AA AC CA AA CC NN AC AA")
        assertEquals(listOf("A"), result.columnLabels)
        assertEquals(listOf(a), result.columnAlleles)
        assertEquals(c, result.comparison)
        assertEquals(listOf(listOf(2.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0)), result.columns)
    }

    @Test
    fun glmClassOrderPutsAcgtFirstAndMlmPutsHomozygotesThenHeterozygotesBySignedByte() {
        // AA, ZZ (0xEE, a negative Java byte), AC.
        val site = "AA ZZ AC AA ZZ AC"
        assertEquals(listOf("A", "M", "Z"), encode(site).classes)
        val mlm = encode(site, MarkerOptions(classOrder = ClassOrder.Mlm))
        assertEquals(listOf("Z", "A", "M"), mlm.classes)
        assertEquals(listOf(1, 0, 2, 1, 0, 2), mlm.classIndex)
        assertEquals(listOf(2, 2, 2), mlm.classCounts)
    }

    @Test
    fun partlyMissingCallsAreMissing() {
        // TASSEL 5 reads A/N as the string A, a homozygote; TASSEL 6 treats it as missing (a recorded departure).
        val result = encode("AA AN AC CC")
        assertEquals(listOf(0, 2, 3), result.observations)
        assertEquals(listOf("A", "C", "M"), result.classes)
        assertEquals(listOf(1, 1, 1), result.classCounts)
        val additive = encode("AA AN AC CC AA", MarkerOptions(encoding = MarkerEncoding.additive))
        assertEquals(listOf(0, 2, 3, 4), additive.observations)
    }

    @Test
    fun monomorphicAndEmptySitesAreSkippedWithReasons() {
        assertEquals(SkipReason.Monomorphic, encode("AA AA AA AA").skip)
        assertEquals(SkipReason.Monomorphic, encode("AA NN AA AN").skip)
        assertEquals(SkipReason.NoObservations, encode("NN NN AN NN").skip)
        // Two classes that are one allele: Z homozygotes score no allele copies.
        assertEquals(
            SkipReason.Monomorphic,
            encode("AA ZZ AA ZZ", MarkerOptions(encoding = MarkerEncoding.additive)).skip,
        )
        assertEquals(listOf("A", "Z"), encode("AA ZZ AA ZZ").classes)
    }

    @Test
    fun indelHeterozygotesAreOneClassAndAlleleCountsIncludeTheGap() {
        val site = "A- AA -- AA"
        val genotype = encode(site)
        assertEquals(listOf("A", "-", "0"), genotype.classes)
        assertEquals(listOf("A:A", "-:-", "A:-"), genotype.classAlleles)
        val additive = encode(site, MarkerOptions(encoding = MarkerEncoding.additive))
        // A has 5 copies and the gap 3, so the gap is the minor (effect) allele.
        assertEquals(listOf("-"), additive.columnLabels)
        assertEquals(listOf(listOf(1.0, 0.0, 2.0, 0.0)), additive.columns)
        assertEquals(a, additive.comparison)
        // A+ and A- are both class 0 at a three-allele site.
        val three = encode("A+ AA ++ A-")
        assertEquals(listOf("A", "+", "0"), three.classes)
        assertEquals(listOf(1, 1, 2), three.classCounts)
        assertEquals(
            SkipReason.Multiallelic,
            encode("A+ AA ++ A-", MarkerOptions(encoding = MarkerEncoding.additive)).skip,
        )
    }

    @Test
    fun multiallelicSitesAreSkippedOrExpandedNeverCollapsed() {
        val site = "AA AC AG GG CC AA"
        // Copies: A 6, C 3, G 3. The omnibus test has five classes; add_F inputs are A and C (G is last).
        val genotype = encode(site)
        assertEquals(listOf("A", "C", "G", "M", "R"), genotype.classes)
        assertEquals(listOf("A", "C"), genotype.columnLabels)
        assertEquals(g, genotype.comparison)
        assertEquals(listOf(2.0, 1.0, 1.0, 0.0, 0.0, 2.0), genotype.columns[0])
        assertEquals(listOf(0.0, 1.0, 0.0, 0.0, 2.0, 0.0), genotype.columns[1])

        val additive = MarkerOptions(encoding = MarkerEncoding.additive)
        assertEquals(SkipReason.Multiallelic, encode(site, additive).skip)
        val expanded = encode(site, additive.copy(multiallelic = MultiallelicAdditive.Expand))
        assertEquals(listOf("C", "G"), expanded.columnLabels)
        assertEquals(a, expanded.comparison)
        assertEquals(listOf(0.0, 1.0, 0.0, 0.0, 2.0, 0.0), expanded.columns[0])
        assertEquals(listOf(0.0, 0.0, 1.0, 2.0, 0.0, 0.0), expanded.columns[1])
        // Additive plus dominance is a biallelic test.
        val both =
            MarkerOptions(encoding = MarkerEncoding.additiveDominance, multiallelic = MultiallelicAdditive.Expand)
        assertEquals(SkipReason.Multiallelic, encode(site, both).skip)
    }

    @Test
    fun additiveUsesTheMinorAlleleUnlessTold() {
        val site = "AA AA AC CC AA AC"
        // A has 8 copies, C has 4.
        val minor = encode(site, MarkerOptions(encoding = MarkerEncoding.additive))
        assertEquals(listOf("C"), minor.columnLabels)
        assertEquals(a, minor.comparison)
        assertEquals(listOf(listOf(0.0, 0.0, 1.0, 2.0, 0.0, 1.0)), minor.columns)
        val major = encode(site, MarkerOptions(encoding = MarkerEncoding.additive, effectAllele = EffectAllele.Major))
        assertEquals(listOf("A"), major.columnLabels)
        assertEquals(c, major.comparison)
        assertEquals(listOf(listOf(2.0, 2.0, 1.0, 0.0, 2.0, 1.0)), major.columns)
        val named =
            encode(site, MarkerOptions(encoding = MarkerEncoding.additive, effectAllele = EffectAllele.Named(a)))
        assertEquals(major.columns, named.columns)
        assertEquals(
            SkipReason.EffectAlleleAbsent,
            encode(site, MarkerOptions(encoding = MarkerEncoding.additive, effectAllele = EffectAllele.Named(g))).skip,
        )
    }

    @Test
    fun additiveFrequenciesComeFromTheAnalysisObservations() {
        // Over all six taxa A has 7 copies and C has 5, so C is minor; over observations 0 to 3 A has 3 and C has 5,
        // so A is minor.
        val sites = listOf("AA AC CC CC AA AC")
        val options = MarkerOptions(encoding = MarkerEncoding.additive)
        assertEquals(listOf("C"), encodeSites(sites, options)[0].columnLabels)
        val subset = encodeSites(sites, options, observations = { intArrayOf(1, 2, 3) })[0]
        assertEquals(listOf("A"), subset.columnLabels)
        assertEquals(listOf(listOf(1.0, 0.0, 0.0)), subset.columns)
        assertEquals(listOf(0, 1, 2), subset.observations)
    }

    @Test
    fun additiveDominanceAddsAHeterozygoteIndicator() {
        val result = encode("AA AC CC AC CA NN", MarkerOptions(encoding = MarkerEncoding.additiveDominance))
        assertEquals(listOf("C", "dom"), result.columnLabels)
        assertEquals(listOf(c, Allele.UNKNOWN), result.columnAlleles)
        assertEquals(listOf(0, 1, 2, 3, 4), result.observations)
        assertEquals(listOf(0.0, 1.0, 2.0, 1.0, 1.0), result.columns[0])
        assertEquals(listOf(0.0, 1.0, 0.0, 1.0, 1.0), result.columns[1])
    }

    @Test
    fun minClassSizeSetsSmallClassesToMissingOrSkipsTheSite() {
        val site = "AA AA AA CC CC AC"
        // Classes A 3, C 2, M 1; the second largest is 2.
        val kept = encode(site, MarkerOptions(minClassSize = 2))
        assertNull(kept.skip)
        assertEquals(listOf("A", "C"), kept.classes)
        assertEquals(listOf(3, 2), kept.classCounts)
        assertEquals(listOf(0, 1, 2, 3, 4), kept.observations)
        assertEquals(2, kept.minorObs)
        // Only one class is big enough.
        val skipped = encode(site, MarkerOptions(minClassSize = 3))
        assertEquals(SkipReason.ClassSize, skipped.skip)
        assertEquals(2, skipped.minorObs)
        // Nothing is small enough to drop.
        assertEquals(3, encode(site, MarkerOptions(minClassSize = 1)).classes.size)
        // The effect allele is chosen among the observations that stay: the heterozygote is gone, so A and C only.
        val additive = encode(site, MarkerOptions(encoding = MarkerEncoding.additive, minClassSize = 2))
        assertEquals(listOf(0, 1, 2, 3, 4), additive.observations)
        assertEquals(listOf("C"), additive.columnLabels)
    }

    @Test
    fun biallelicOnlyKeepsTwoClassesOrThreeWithOneHeterozygote() {
        val options = MarkerOptions(biallelicOnly = true)
        assertNull(encode("AA CC AC AA", options).skip)
        assertNull(encode("AA CC CC", options).skip)
        assertEquals(SkipReason.NotBiallelic, encode("AA CC GG", options).skip)
        assertEquals(SkipReason.NotBiallelic, encode("AA AC AG", options).skip)
        assertEquals(SkipReason.NotBiallelic, encode("AA CC AC AG", options).skip)
        // Phased heterozygotes are one class, so this is still biallelic.
        assertNull(encode("AA CC AC CA", options).skip)
    }

    @Test
    fun minMafAndMinMacUseTheAnalysisObservations() {
        // 15 A copies and one C copy among 16.
        val site = "AA AA AA AA AA AA AA AC"
        assertNull(encode(site).skip)
        assertEquals(SkipReason.MinorAlleleCount, encode(site, MarkerOptions(minMac = 2)).skip)
        assertNull(encode(site, MarkerOptions(minMac = 1)).skip)
        assertEquals(SkipReason.MinorAlleleFrequency, encode(site, MarkerOptions(minMaf = 0.1)).skip)
        assertNull(encode(site, MarkerOptions(minMaf = 0.05)).skip)
        // Missing calls are not counted, so the frequency is 1 in 4 here.
        assertNull(encode("AA AC NN NN", MarkerOptions(minMaf = 0.2)).skip)
        // The same marker passes for one trait and is skipped for another.
        val data = alignedWith(genotypes(listOf("AA AA AA AC")))
        val source = GenotypeMarkerSource.of(data)
        val columns = MarkerColumns(4)
        val block = source.scan(listOf(intArrayOf(0, 1, 2, 3), intArrayOf(0, 1, 2)), MarkerOptions(minMac = 1)).first()
        assertTrue(block.load(0, 0, columns))
        assertFalse(block.load(0, 1, columns))
        assertEquals(SkipReason.MinorAlleleCount, columns.skip)
    }

    @Test
    fun markersSpanBlocksInOrder() {
        val sites = listOf("AA AC CC", "AA AA CC", "CC AC AC", "AA NN NN", "AA CC AC")
        val whole = encodeSites(sites)
        val blocked = encodeSites(sites, blockSize = 2)
        assertEquals(whole, blocked)
        assertEquals(listOf(null, null, null, SkipReason.Monomorphic, null), whole.map { it.skip })
    }

    @Test
    fun observationSetsReplicatesAndTaxaWithoutGenotypesAreHandled() {
        // Taxa t0..t3 are genotyped; the phenotype has t1 twice, t3, and x with no genotype (a union join also adds t0
        // and t2 as genotype-only observations).
        val table = genotypes(listOf("AA AC CC AC", "AA AA AC CC"))
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(listOf("t1", "t1", "t3", "x")),
                NumericAttribute("y", AttributeType.data, floatArrayOf(1f, 2f, 3f, 4f)),
            ),
        )
        val data = alignGenotypePhenotype(table, phenotype, JoinType.union)
        assertEquals(6, data.numberOfObservations)
        val source = GenotypeMarkerSource.of(data)
        val columns = MarkerColumns(6)
        // Trait 1 uses every observation, trait 2 only observations 1, 2, 4 (t1 replicate, t3, t0).
        val block = source.scan(listOf(IntArray(6) { it }, intArrayOf(1, 2, 4))).first()
        assertTrue(block.load(0, 0, columns))
        val first = snapshot(columns)
        // x has no genotype and is missing; t1 appears twice with the same class.
        assertEquals(listOf(0, 1, 2, 4, 5), first.observations)
        assertEquals(listOf("A", "C", "M"), first.classes)
        assertEquals(listOf(1, 1, 3), first.classCounts)
        assertTrue(block.load(0, 1, columns))
        val second = snapshot(columns)
        assertEquals(listOf(0, 1, 2), second.observations)
        assertEquals(listOf("A", "M"), second.classes)
        assertEquals(listOf(1, 2), second.classCounts)
        assertEquals(listOf(1, 1, 0), second.classIndex)
    }

    @Test
    fun eachSiteIsReadOnceNoMatterHowManyTraitsShareTheScan() {
        val inner = genotypes(listOf("AA AC CC AC CA", "AA AA AC CC AA", "CC AC AC AA AA", "AA CC AC CA AA"))
        val counting = Counting(inner)
        val data = alignedWith(counting)
        val source = GenotypeMarkerSource.of(data)
        val traits = listOf(intArrayOf(0, 1, 2, 3, 4), intArrayOf(0, 2, 4), intArrayOf(1, 3))
        val columns = MarkerColumns(5)
        var loads = 0
        source.scan(traits, blockSize = 3).forEach { block ->
            for (marker in block.markers) {
                // Two passes over each trait: the second must not read again.
                repeat(2) { for (set in traits.indices) block.load(marker, set, columns).also { loads++ } }
            }
        }
        assertEquals(4 * 3 * 2, loads)
        assertEquals(4 * 5, counting.cellsRead, "each of 4 sites x 5 taxa is read exactly once")
        assertEquals(0, counting.wholeSiteOrTaxonReads)
    }

    @Test
    fun onlyGenotypedTaxaAreRead() {
        val inner = genotypes(listOf("AA AC CC AC CA", "AA AA AC CC AA"))
        val counting = Counting(inner)
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(listOf("t4", "t1", "t4")),
                NumericAttribute("y", AttributeType.data, floatArrayOf(1f, 2f, 3f)),
            ),
        )
        val data = alignGenotypePhenotype(counting, phenotype)
        val source = GenotypeMarkerSource.of(data)
        val columns = MarkerColumns(3)
        source.scan(listOf(intArrayOf(0, 1, 2))).forEach { block ->
            for (marker in block.markers) block.load(marker, 0, columns)
        }
        assertEquals(2 * 2, counting.cellsRead, "two sites x two distinct genotyped taxa")
    }

    @Test
    fun dosageComesFromReferenceProbabilityAndAMissingCallIsMissing() = runTest {
        // Site 0: major allele A (3 copies, tied with C, A sorts first). Site 1: major allele A.
        // GLM's dosage (NumericalGenotypePlugin without PCA's imputation) leaves a missing call missing, so t3 and t4
        // are dropped at site 0 and t4 at site 1; the golden glm-refprob drops observations the same way.
        val table = genotypes(listOf("AA AC CC NN NN", "AA AA AC CC NN"))
        val data = alignedWith(table)
        val source = GenotypeMarkerSource.withDosage(data)
        val columns = MarkerColumns(5)
        val results = ArrayList<Encoded>()
        source.scan(listOf(IntArray(5) { it }), MarkerOptions(encoding = MarkerEncoding.dosage)).forEach { block ->
            for (marker in block.markers) {
                assertTrue(block.load(marker, 0, columns))
                results.add(snapshot(columns))
            }
        }
        assertEquals(listOf(0, 1, 2), results[0].observations)
        assertEquals(listOf(listOf(1.0, 0.5, 0.0)), results[0].columns)
        assertEquals(listOf("dosage"), results[0].columnLabels)
        assertEquals(listOf(0, 1, 2, 3), results[1].observations)
        assertEquals(listOf(listOf(1.0, 1.0, 0.5, 0.0)), results[1].columns)
        assertEquals(0, results[0].minorObs)
    }

    @Test
    fun dosageNeedsASourceMadeWithDosageAndCallsNeedNoSuchSource() {
        val source = GenotypeMarkerSource.of(alignedWith(genotypes(listOf("AA AC"))))
        assertFailsWith<IllegalArgumentException> {
            source.scan(listOf(intArrayOf(0, 1)), MarkerOptions(encoding = MarkerEncoding.dosage))
        }
    }

    @Test
    fun theSourceDescribesItsMarkersAndChecksArguments() {
        val data = alignedWith(genotypes(listOf("AA AC", "CC AC")))
        val source = GenotypeMarkerSource.of(data)
        assertEquals(2, source.numberOfMarkers)
        assertEquals("S1_200", source.position(1).name)
        assertFailsWith<IllegalArgumentException> { source.scan(listOf(intArrayOf(5))).toList() }
        assertFailsWith<IllegalArgumentException> { source.scan(emptyList(), blockSize = 0) }
        val block = source.scan(listOf(intArrayOf(0, 1))).first()
        assertFailsWith<IllegalArgumentException> { block.load(0, 0, MarkerColumns(1)) }
        assertFailsWith<IllegalArgumentException> { block.load(7, 0, MarkerColumns(2)) }
        assertFailsWith<IllegalArgumentException> { block.load(0, 3, MarkerColumns(2)) }
        assertFailsWith<IllegalArgumentException> { MarkerOptions(minMaf = 0.9) }
    }

    @Test
    fun explicitAlleleStringsLabelClassesAndColumns() {
        val result = encode("GG GT TT GG", MarkerOptions(classOrder = ClassOrder.Glm))
        assertEquals(listOf("G", "T", "K"), result.classes)
        assertEquals(listOf("G:G", "T:T", "G:T"), result.classAlleles)
    }

    /** A table that counts how often genotypes are read, and does not override `siteBlocks`. */
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
}
