package net.maizegenetics.tassel.association

import kotlinx.coroutines.CoroutineStart
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.coroutineScope
import kotlinx.coroutines.launch
import kotlinx.coroutines.test.runTest
import kotlinx.coroutines.withContext
import kotlinx.coroutines.yield
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class AssociationSinkTest {
    private val provenance = AnalysisProvenance(operation = "glm", version = "test")

    private fun record(
        marker: Int,
        p: Double,
        trait: String = "t",
        test: AssociationTestKind = AssociationTestKind.Genotype,
        status: AssociationStatus = AssociationStatus.Ok,
        reason: SkipReason? = null,
    ) = AssociationRecord(
        trait = trait, markerIndex = marker, marker = "m$marker", chromosome = "1", position = marker,
        test = test, p = p, status = status, skipReason = reason,
    )

    private fun effect(marker: Int, trait: String = "t") = AssociationEffect(
        trait = trait,
        markerIndex = marker,
        marker = "m$marker",
        chromosome = "1",
        position = marker,
        allele = "A",
        estimate = marker.toDouble(),
    )

    private fun chunks(count: Int, perChunk: Int = 3) = (0 until count).map { c ->
        val markers = (c * perChunk until (c + 1) * perChunk)
        AssociationChunk(
            sequence = c,
            records = markers.map { record(it, p = (it % 10) / 10.0) },
            effects = markers.map { effect(it) },
        )
    }

    private fun AssociationTable.markers() = records.map { it.markerIndex }.toList()

    @Test
    fun orderedSinkReleasesChunksInSequenceWhateverTheArrivalOrder() = runTest {
        val collector = AssociationCollector()
        for (chunk in chunks(20).shuffled(kotlin.random.Random(7))) collector.accept(chunk)
        val table = collector.finish(provenance)
        assertEquals((0 until 60).toList(), table.markers())
        assertEquals((0 until 60).toList(), table.effects.map { it.markerIndex }.toList())
    }

    @Test
    fun concurrentProducersGiveTheSameTableAsOneProducer() = runTest {
        val all = chunks(200)
        val serial = AssociationCollector().also { c -> all.forEach { c.accept(it) } }.finish(provenance)

        val concurrent = AssociationCollector()
        withContext(Dispatchers.Default) {
            coroutineScope {
                for (chunk in all.reversed()) {
                    launch {
                        yield()
                        concurrent.accept(chunk)
                    }
                }
            }
        }
        val table = concurrent.finish(provenance)
        assertEquals(serial.records.toList(), table.records.toList())
        assertEquals(serial.effects.toList(), table.effects.toList())
        assertEquals(serial.summary, table.summary)
    }

    @Test
    fun producersStartedOutOfOrderStillGiveSequenceOrder() = runTest {
        val collector = AssociationCollector()
        coroutineScope {
            for (chunk in chunks(10).reversed()) launch(start = CoroutineStart.UNDISPATCHED) { collector.accept(chunk) }
        }
        assertEquals((0 until 30).toList(), collector.finish(provenance).markers())
    }

    @Test
    fun aMissingChunkIsAnErrorAtFinish() = runTest {
        val collector = AssociationCollector()
        collector.accept(chunks(3)[2])
        assertFailsWith<IllegalStateException> { collector.finish(provenance) }
    }

    @Test
    fun aRepeatedChunkIsRefused() = runTest {
        val collector = AssociationCollector()
        val first = chunks(2)[0]
        collector.accept(first)
        assertFailsWith<IllegalArgumentException> { collector.accept(first) }
    }

    @Test
    fun glmMaxPKeepsRowsAtOrBelowTheCutoffAndDropsNaN() = runTest {
        val rows = listOf(record(0, 0.01), record(1, 0.0100001), record(2, Double.NaN), record(3, 0.0), record(4, 0.5))
        val table = AssociationCollector(0.01, MaxPRule.Glm).also {
            it.accept(AssociationChunk(0, rows))
        }.finish(provenance)
        assertEquals(listOf(0, 3), table.markers())
    }

    @Test
    fun mlmMaxPDropsOnlyRowsAboveTheCutoffAndKeepsNaN() = runTest {
        val rows = listOf(record(0, 0.01), record(1, 0.0100001), record(2, Double.NaN), record(3, 0.0), record(4, 0.5))
        val table = AssociationCollector(0.01, MaxPRule.Mlm).also {
            it.accept(AssociationChunk(0, rows))
        }.finish(provenance)
        assertEquals(listOf(0, 2, 3), table.markers())
    }

    @Test
    fun maxPOfOneKeepsEverythingUnderBothRules() = runTest {
        val rows = listOf(record(0, 1.0), record(1, Double.NaN), record(2, 0.3))
        for (rule in MaxPRule.entries) {
            val table = AssociationCollector(1.0, rule).also { it.accept(AssociationChunk(0, rows)) }.finish(provenance)
            assertEquals(3, table.size, rule.name)
        }
    }

    @Test
    fun maxPRejectsOutOfRangeValues() {
        for (bad in listOf(-0.1, 1.5, Double.NaN)) {
            assertFailsWith<IllegalArgumentException> { MaxPSink(bad, MaxPRule.Glm, TableSink()) }
        }
    }

    @Test
    fun theFilterDecidesByTheOmnibusRecordAndTakesTheMarkersOtherRecordsAndEffectsWithIt() = runTest {
        val rows = listOf(
            record(0, 0.9, test = AssociationTestKind.Genotype),
            record(0, 0.001, test = AssociationTestKind.Additive),
            record(1, 0.001, test = AssociationTestKind.Genotype),
            record(1, 0.9, test = AssociationTestKind.Additive),
            // Another trait's marker 1 has its own decision.
            record(1, 0.9, trait = "u"),
        )
        val chunk = AssociationChunk(0, rows, listOf(effect(0), effect(1), effect(1, "u")))
        val table = AssociationCollector(0.01).also { it.accept(chunk) }.finish(provenance)
        assertEquals(listOf(1, 1), table.records.map { it.markerIndex }.toList())
        assertEquals(
            listOf(AssociationTestKind.Genotype, AssociationTestKind.Additive),
            table.records.map {
                it.test
            }.toList(),
        )
        assertEquals(listOf("t" to 1), table.effects.map { it.trait to it.markerIndex }.toList())
    }

    @Test
    fun theSummaryCountsEverythingEvenWhenTheFilterDropsRows() = runTest {
        val rows = listOf(
            record(0, 0.001),
            record(1, 0.5),
            record(2, Double.NaN, status = AssociationStatus.Skipped, reason = SkipReason.Monomorphic),
            record(3, Double.NaN, status = AssociationStatus.Skipped, reason = SkipReason.Monomorphic),
            record(4, Double.NaN, status = AssociationStatus.Skipped, reason = SkipReason.ClassSize),
            record(5, Double.NaN, status = AssociationStatus.Failed),
        )
        val table = AssociationCollector(0.01).also { it.accept(AssociationChunk(0, rows)) }.finish(provenance)
        assertEquals(1, table.size)
        assertEquals(
            AssociationSummary(
                records = 6,
                tested = 2,
                skipped = mapOf(SkipReason.Monomorphic to 2L, SkipReason.ClassSize to 1L),
                failed = 1,
            ),
            table.summary,
        )
        assertEquals(3L, table.summary.skippedTotal)
    }

    @Test
    fun maxPOnTheUnfilteredGoldenReproducesTheMaxPGolden() = runTest {
        fun text(name: String) = normalizeLineEndings(goldenText("golden/glm-options/$name"))
        val parsed = parseGlm(text("append-add-dom.stats.full.txt"), text("append-add-dom.effects.full.txt"))
        // Deliver by marker group, 40 groups per chunk, in reverse arrival order.
        val groups = parsed.records.chunked(3)
        val effectsByGroup = parsed.effects.groupBy { it.markerIndex }
        val chunks = groups.chunked(40).mapIndexed { sequence, batch ->
            AssociationChunk(
                sequence,
                batch.flatten(),
                batch.flatMap { effectsByGroup[it.first().markerIndex].orEmpty() },
            )
        }
        val collector = AssociationCollector(0.01, MaxPRule.Glm)
        coroutineScope { for (chunk in chunks.reversed()) launch { collector.accept(chunk) } }
        val table = collector.finish(provenance)

        val statistics = buildString {
            table.writeGlmStatistics(this, GlmLayout.Genotype, NumberForm.Full, false, appendAddDom = false)
        }
        val effects = buildString { table.writeGlmEffects(this, GlmLayout.Genotype, NumberForm.Full) }
        assertEquals(text("max-p-0.01.stats.full.txt"), statistics)
        assertEquals(text("max-p-0.01.effects.full.txt"), effects)
        assertTrue(table.summary.records == parsed.records.size.toLong())
    }
}
