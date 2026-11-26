package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.table.doubleToString
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

class AssociationTableTest {
    private val provenance = AnalysisProvenance(
        operation = "glm",
        version = "6.0.0-test",
        parameters = linkedMapOf("maxP" to "0.01", "encoding" to "genotype"),
        inputs = listOf(
            AnalysisInput("genotypes.hmp.txt", taxa = 281, markers = 3093),
            AnalysisInput("traits", observations = 270),
        ),
        model = "EarHT ~ Q1 + Q2 + marker",
        kinship = "kinship.txt, centered IBS",
        seed = 42,
        excludedObservations = linkedMapOf("missing trait" to 11),
        excludedMarkers = linkedMapOf("Monomorphic" to 4),
    )

    private val full = AssociationRecord(
        trait = "EarHT", markerIndex = 7, marker = "PZA1", chromosome = "9", position = 123,
        model = AssociationModel.Mlm, test = AssociationTestKind.Additive, encoding = MarkerEncoding.additive,
        effectAllele = "A", comparisonAllele = "C", observations = 250, taxa = 251, alleleCount = 90,
        alleleFrequency = 0.18, minorAlleleCount = 90, minorObs = 12, effect = -1.5, standardError = 0.5,
        ciLower = -2.5, ciUpper = -0.5, statistic = 9.0, numeratorDf = 1.0, denominatorDf = 248.0, p = 0.003,
        log10P = -2.52, rSquared = 0.03, markerMeanSquare = 1.0, errorMeanSquare = 2.0, modelDf = 3.0,
        modelMeanSquare = 4.0, permutationP = 0.004, geneticVariance = 5.0, residualVariance = 6.0,
        heritability = 0.45, varianceEstimation = VarianceEstimation.P3d,
    )

    @Test
    fun everyFieldOfARecordSurvivesTheColumns() {
        val skipped = full.copy(
            status = AssociationStatus.Skipped,
            skipReason = SkipReason.MinorAlleleCount,
            message = "fewer than 5",
            statistic = Double.NaN,
        )
        val failed = full.copy(status = AssociationStatus.Failed, message = "singular", p = Double.NEGATIVE_INFINITY)
        val table = AssociationTable.of(listOf(full, skipped, failed), provenance = provenance)
        assertEquals(listOf(full, skipped, failed).map { it.toString() }, table.records.map { it.toString() }.toList())
        assertEquals(3, table.size)
    }

    @Test
    fun aRecordNeedsASkipReasonExactlyWhenItIsSkipped() {
        assertFailsWith<IllegalArgumentException> { full.copy(status = AssociationStatus.Skipped) }
        assertFailsWith<IllegalArgumentException> { full.copy(skipReason = SkipReason.Monomorphic) }
    }

    @Test
    fun effectsSurviveTheColumns() {
        val effect = AssociationEffect("EarHT", 7, "PZA1", "9", 123, "A", 54, -8.9, 0.25)
        val table = AssociationTable.of(emptyList(), listOf(effect), provenance)
        assertEquals(effect, table.effect(0))
        assertEquals(1, table.effectCount)
    }

    @Test
    fun readingPastTheEndIsRefused() {
        val table = AssociationTable.of(listOf(full), provenance = provenance)
        assertFailsWith<IllegalArgumentException> { table[1] }
        assertFailsWith<IllegalArgumentException> { table.effect(0) }
    }

    @Test
    fun rowsSpanSegmentsAndReportsListThemTraitByTrait() {
        val rows = (0 until 20_000).map { i ->
            full.copy(trait = if (i % 2 == 0) "b" else "a", markerIndex = i, marker = "m$i", p = i.toDouble())
        }
        val table = AssociationTable.of(rows, provenance = provenance)
        assertEquals(20_000, table.size)
        assertEquals(rows[19_999].toString(), table[19_999].toString())
        assertEquals(listOf("b", "a"), table.traits)
        val order = table.recordOrder()
        // "b" was seen first, so it is listed first; rows keep append order within a trait.
        assertEquals(listOf(0, 2, 4), order.take(3))
        assertEquals(10_000, order.indexOfFirst { rows[it].trait == "a" })
        assertEquals(1, table.recordOrder().let { it[10_000] })
    }

    @Test
    fun theBuilderCannotBeReusedAfterBuild() {
        val builder = AssociationTableBuilder()
        builder.build(provenance)
        assertFailsWith<IllegalStateException> { builder.append(full) }
        assertFailsWith<IllegalStateException> { builder.build(provenance) }
    }

    @Test
    fun provenanceWritesDeterministicHeaderLines() {
        assertEquals(
            listOf(
                "##tassel_version=6.0.0-test",
                "##operation=glm",
                "##parameter.maxP=0.01",
                "##parameter.encoding=genotype",
                "##input=genotypes.hmp.txt taxa=281 markers=3093",
                "##input=traits observations=270",
                "##model=EarHT ~ Q1 + Q2 + marker",
                "##kinship=kinship.txt, centered IBS",
                "##seed=42",
                "##excluded_observations.missing trait=11",
                "##excluded_markers.Monomorphic=4",
            ),
            provenance.headerLines(),
        )
        assertEquals(
            listOf("##tassel_version=v", "##operation=x"),
            AnalysisProvenance(operation = "x", version = "v").headerLines(),
        )
        assertEquals(
            "##model=a b c d",
            AnalysisProvenance(operation = "x", version = "v", model = "a\nb\rc\td").headerLines().last(),
        )
    }

    @Test
    fun modernStatisticsCarryTheSchemaProvenanceAndEveryColumn() {
        val skipped = AssociationRecord(
            trait = "t",
            markerIndex = 1,
            marker = "m1",
            chromosome = "1",
            position = 5,
            status = AssociationStatus.Skipped,
            skipReason = SkipReason.Monomorphic,
            message = "one\tclass",
        )
        val table = AssociationTable.of(listOf(full, skipped), provenance = provenance)
        val lines = buildString { table.writeModernStatistics(this) }.trimEnd('\n').split('\n')
        assertEquals("#schema=association/1", lines[0])
        assertEquals(provenance.headerLines(), lines.subList(1, 1 + provenance.headerLines().size))
        val header = lines[1 + provenance.headerLines().size].split('\t')
        assertEquals(MODERN_COLUMNS, header)
        val rows = lines.drop(2 + provenance.headerLines().size).map { it.split('\t') }
        assertEquals(2, rows.size)
        rows.forEach { assertEquals(header.size, it.size) }
        val byName = { row: List<String>, name: String -> row[header.indexOf(name)] }
        assertEquals("-2.52", byName(rows[0], "log10_p"))
        assertEquals("0.003", byName(rows[0], "p"))
        assertEquals("p3d", byName(rows[0], "variance_estimation"))
        assertEquals("additive", byName(rows[0], "test"))
        assertEquals("skipped", byName(rows[1], "status"))
        assertEquals("Monomorphic", byName(rows[1], "skip_reason"))
        assertEquals("one class", byName(rows[1], "message"))
        assertEquals("", byName(rows[1], "p"))
        assertEquals("", byName(rows[1], "observations"))
    }

    @Test
    fun modernNumbersRoundTripAtFullPrecision() {
        val p = 1.2345678901234567E-300
        val table = AssociationTable.of(listOf(full.copy(p = p, log10P = -299.9)), provenance = provenance)
        val text = buildString { table.writeModernStatistics(this) }
        assertTrue(doubleToString(p) in text)
        assertEquals(p, doubleToString(p).toDouble())
    }

    @Test
    fun modernEffectsCarryTheirOwnSchema() {
        val table = AssociationTable.of(
            emptyList(),
            listOf(
                AssociationEffect("t", 3, "m3", "2", 9, "A", 54, -8.5),
                AssociationEffect("t", 3, "m3", "2", 9, "C", estimate = 0.0),
            ),
            provenance,
        )
        val lines = buildString { table.writeModernEffects(this) }.trimEnd('\n').split('\n')
        assertEquals("#schema=association-effects/1", lines[0])
        val n = provenance.headerLines().size
        assertEquals(
            "trait\tmarker_index\tmarker\tchr\tpos\tallele\tobservations\testimate\tstandard_error",
            lines[
                1 +
                    n,
            ],
        )
        assertEquals("t\t3\tm3\t2\t9\tA\t54\t-8.5\t", lines[2 + n])
        assertEquals("t\t3\tm3\t2\t9\tC\t\t0.0\t", lines[3 + n])
    }

    @Test
    fun theLegacyWritersLeaveOutProvenance() {
        val table = AssociationTable.of(listOf(full.copy(test = AssociationTestKind.Genotype)), provenance = provenance)
        val stats = buildString { table.writeGlmStatistics(this) }
        assertTrue(stats.startsWith("Trait\tMarker\tChr\tPos\tmarker_F\tp\tmarker_Rsq\t"))
        assertTrue("##" !in stats)
        assertNull(stats.lines().firstOrNull { it.startsWith("#") })
    }

    @Test
    fun theResultTableViewIsBounded() {
        val table = AssociationTable.of((0 until 50).map { full.copy(markerIndex = it) }, provenance = provenance)
        assertEquals(10, table.toResultTable(maxRows = 10).rowCount)
        assertEquals(50, table.toResultTable().rowCount)
    }

    @Test
    fun aMillionRecordsStayUnderTheStatedBound() {
        val provenance = AnalysisProvenance(operation = "glm", version = "test")
        val builder = AssociationTableBuilder()
        val rows = 1_000_000
        // Records are made and dropped one at a time, so what is held is only what the table stores.
        for (i in 0 until rows) {
            builder.append(
                AssociationRecord(
                    trait = if (i <
                        rows / 2
                    ) {
                        "t1"
                    } else {
                        "t2"
                    },
                    markerIndex = i % (rows / 2), marker = "m${i % (rows / 2)}",
                    chromosome = (i % 10 + 1).toString(), position = i, statistic = i.toDouble(), p = 1.0 / (i + 1),
                    numeratorDf = 2.0, denominatorDf = 250.0, minorObs = 5,
                ),
            )
        }
        val table = builder.build(provenance)
        assertEquals(rows, table.size)
        // 19 double columns (152 bytes) + 13 int columns (52) + 6 byte columns (6) = 210 bytes per record, plus about
        // 30 MB for the 500,000 distinct marker names, estimated at 48 bytes and two per character.
        val bound = 260_000_000L
        assertTrue(table.estimatedBytes < bound, "${table.estimatedBytes} bytes for $rows records, bound $bound")
        assertTrue(table.estimatedBytes > 200L * rows)
    }
}
