package net.maizegenetics.tassel.association

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.TestData
import kotlin.math.abs
import kotlin.math.max
import kotlin.random.Random
import kotlin.test.Test
import kotlin.time.Duration
import kotlin.time.measureTime

/**
 * Microbenchmarks of the GLM kernels for the tassel-6#192 handoff: the reference kernel (#191), and the fast kernel on
 * one worker and on every core, on the tutorial (281 x 3,093, three traits, Q1 and Q2) and a synthetic 1,000 x 100,000
 * panel (three traits, two sharing an observation set, three covariates, 2% missing calls). Genotypes are in memory and
 * reading them is not timed; records go to a counting sink, so the time is the scan, the fits, and building records.
 * Also reports how far the fast kernel's statistics are from the reference's.
 *
 * Stays out of CI. Run with
 * `TASSEL_BENCHMARK=1 ./gradlew :tassel-core:jvmTest --tests '*GlmBenchmarkTest*' --rerun -i`; `TASSEL_BENCHMARK_RUNS`
 * sets the timed runs per row (default 3, after one warm-up), and the median is reported.
 */
class GlmBenchmarkTest {
    @Test
    fun glmKernels() {
        if (System.getenv("TASSEL_BENCHMARK") != "1") return
        val runs = System.getenv("TASSEL_BENCHMARK_RUNS")?.toInt() ?: 3
        runBlocking(Dispatchers.Default) {
            println(
                "| Dataset | Kernel | Workers | Median | Marker-trait fits/s | Paths (full / indicators / pattern hit / miss) |",
            )
            println("|---|---|---|---|---|---|")
            bench("tutorial 281 x 3,093", tutorial(), runs)
            bench(
                "synthetic 1,000 x 100,000",
                syntheticGlmPanel(taxa = 1_000, markers = 100_000),
                runs,
                referenceRuns = 1,
            )
        }
    }

    private suspend fun bench(name: String, data: AssociationData, runs: Int, referenceRuns: Int = runs) {
        val options = GlmOptions()
        val fits = data.traits.size.toLong() * GenotypeMarkerSource.of(data.data).numberOfMarkers
        val reference = time(referenceRuns) { glmReference(data, GenotypeMarkerSource.of(data.data), options, Discard) }
        println(row(name, "reference", "1", reference, fits, ""))
        val cores = associationWorkerCount()
        for (workers in listOf(1, cores).distinct()) {
            val stats = GlmFastStats()
            val fast = time(runs) {
                stats.reset()
                glm(data, GenotypeMarkerSource.of(data.data), options, Discard, workers = workers, stats = stats)
            }
            val paths = "${stats.fullSet} / ${stats.indicators} / ${stats.patternHits} / ${stats.patternMisses}"
            println(row(name, "fast", "$workers", fast, fits, paths) + "  (${speedup(reference, fast)} x reference)")
        }
        val expected = Capture().also { glmReference(data, GenotypeMarkerSource.of(data.data), options, it) }
        val actual = Capture().also { glm(data, GenotypeMarkerSource.of(data.data), options, it) }
        println("|   agreement: ${expected.agreement(actual)} | | | | | |")
    }

    private fun row(name: String, kernel: String, workers: String, time: Duration, fits: Long, paths: String): String {
        val seconds = time.inWholeNanoseconds / 1e9
        return "| $name | $kernel | $workers | ${"%.3f".format(
            seconds,
        )} s | ${"%,.0f".format(fits / seconds)} | $paths |"
    }

    private fun speedup(reference: Duration, fast: Duration): String =
        "%.1f".format(reference.inWholeNanoseconds.toDouble() / fast.inWholeNanoseconds)

    private suspend fun time(runs: Int, block: suspend () -> Unit): Duration {
        block()
        val times = (0 until runs).map { measureTime { block() } }.sorted()
        return times[times.size / 2]
    }

    private fun GlmFastStats.reset() {
        fullSet = 0
        indicators = 0
        patternHits = 0
        patternMisses = 0
        reference = 0
        skipped = 0
    }

    private object Discard : AssociationSink {
        override suspend fun accept(chunk: AssociationChunk) = Unit
    }

    /** Keeps what [agreement] compares from each record, and the effect estimates, in order. */
    private class Capture : AssociationSink {
        val tests = ArrayList<AssociationTestKind>()
        val statistics = ArrayList<Double>()
        val errorDf = ArrayList<Double>()
        val p = ArrayList<Double>()
        val log10P = ArrayList<Double>()
        val estimates = ArrayList<Double>()

        override suspend fun accept(chunk: AssociationChunk) {
            for (record in chunk.records) {
                tests += record.test
                statistics += record.statistic
                errorDf += record.denominatorDf
                p += record.p
                log10P += record.log10P
            }
            for (effect in chunk.effects) estimates += effect.estimate
        }

        /**
         * How far [other] is from these results, against `GlmFastTest`'s rule: the largest difference of F as a fraction
         * of its tolerance (`rel 1e-8`, or twice D4's eps multiple of the error sum of squares, twice that again for `dom_F`; at most 1
         * passes), the largest relative difference of F where `F >= 1`, of `log10 p` where F agrees to `rel 1e-8`,
         * `p >= 1e-300`, and `|log10 p| >= 1e-6` (the largest absolute difference where it is smaller, near p = 1), and
         * of class estimates of at least `1e-3`.
         */
        fun agreement(other: Capture): String {
            check(statistics.size == other.statistics.size && estimates.size == other.estimates.size)
            var tolerance = 0.0
            var strong = 0.0
            var log = 0.0
            var logNearZero = 0.0
            var estimate = 0.0
            for (i in statistics.indices) {
                val a = statistics[i]
                val b = other.statistics[i]
                if (a.isNaN() && b.isNaN()) continue
                val noise =
                    2 * GlmDepartures.cancellationEpsMultiple * EPSILON * errorDf[i] *
                        if (tests[i] == AssociationTestKind.Dominance) 2 else 1
                tolerance = max(tolerance, abs(a - b) / max(1e-8 * max(abs(a), abs(b)), noise))
                if (abs(a) >= 1.0) strong = max(strong, relative(a, b))
                val fAgrees = abs(a - b) <= 1e-8 * max(abs(a), abs(b))
                if (fAgrees && p[i] >= 1e-300) {
                    if (abs(log10P[i]) >= 1e-6) {
                        log = max(log, relative(log10P[i], other.log10P[i]))
                    } else {
                        logNearZero = max(logNearZero, abs(log10P[i] - other.log10P[i]))
                    }
                }
            }
            for (i in estimates.indices) {
                if (abs(estimates[i]) >= 1e-3) estimate = max(estimate, relative(estimates[i], other.estimates[i]))
            }
            return (
                "F differences at most %.2f of the tolerance; max relative difference: F >= 1 %.1e, log10 p %.1e " +
                    "(%.1e absolute near 0), class estimates %.1e"
                ).format(tolerance, strong, log, logNearZero, estimate)
        }

        private fun relative(a: Double, b: Double): Double {
            if (a.isNaN() && b.isNaN() || a == b) return 0.0
            return abs(a - b) / max(abs(a), abs(b))
        }

        private companion object {
            const val EPSILON = 2.220446049250313E-16
        }
    }

    private suspend fun tutorial(): AssociationData {
        val genotypes = HapMap.read(Buffer().also { it.write(TestData.readBytes("tutorial/mdp_genotype.hmp.txt")) })
        suspend fun phenotype(path: String) =
            PhenotypeFiles.read(Buffer().also { it.write(TestData.readBytes(path)) }).table
        val input = DataSet.of(
            genotypes,
            phenotype("tutorial/mdp_traits.txt"),
            phenotype("tutorial/mdp_population_structure.txt"),
        )
        return AssociationData.of(input, AssociationSpec(covariates = listOf("Q1", "Q2")))
    }
}

/** Biallelic sites, minor allele frequency 5% to 50%, mostly homozygous, 2% missing; traits from 20 causal sites. */
internal suspend fun syntheticGlmPanel(taxa: Int, markers: Int): AssociationData {
    val random = Random(192)
    val names = (0 until taxa).map { "s$it" }
    val alleles = listOf("AA", "CC", "AC", "NN").map { Nucleotide.genotype(it) }
    val genotypes = ByteArray(taxa * markers)
    val copies = Array(20) { DoubleArray(taxa) }
    for (site in 0 until markers) {
        val frequency = 0.05 + 0.45 * random.nextDouble()
        for (taxon in 0 until taxa) {
            val u = random.nextDouble()
            val call = when {
                u < 0.02 -> 3
                u < 0.07 -> 2
                random.nextDouble() < frequency -> 1
                else -> 0
            }
            genotypes[site * taxa + taxon] = alleles[call]
            if (site % (markers / 20) == 0 && site / (markers / 20) < 20) {
                copies[site / (markers / 20)][taxon] = if (call == 1) {
                    2.0
                } else if (call == 2) {
                    1.0
                } else {
                    0.0
                }
            }
        }
    }
    val table = InMemoryGenotypeTable.wrap(
        TaxaList.of(names),
        PositionList((0 until markers).map { Position(Chromosome.of(1 + it / 10_000), it + 1, "S$it") }),
        genotypes,
    )
    val pcs = Array(3) { FloatArray(taxa) { random.nextGaussian().toFloat() } }
    fun trait(name: String, missing: Double) = NumericAttribute(
        name,
        AttributeType.data,
        FloatArray(taxa) { t ->
            if (random.nextDouble() < missing) {
                Float.NaN
            } else {
                (50 + copies.sumOf { it[t] } * 0.3 + pcs[0][t] + 0.5 * pcs[1][t] + random.nextGaussian()).toFloat()
            }
        },
    )
    val phenotype = PhenotypeTable(
        listOf(
            PhenotypeTable.taxa(names),
            trait("y1", 0.0),
            trait("y2", 0.0),
            trait("y3", 0.1),
            NumericAttribute("PC1", AttributeType.covariate, pcs[0]),
            NumericAttribute("PC2", AttributeType.covariate, pcs[1]),
            NumericAttribute("PC3", AttributeType.covariate, pcs[2]),
        ),
    )
    return AssociationData.of(DataSet.of(table, phenotype))
}
