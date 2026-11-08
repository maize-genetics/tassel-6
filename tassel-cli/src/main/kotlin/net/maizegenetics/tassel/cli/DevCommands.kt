package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.Context
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import com.github.ajalt.clikt.parameters.types.long
import com.github.ajalt.clikt.parameters.types.restrictTo
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import net.maizegenetics.tassel.distance.centeredIbs
import net.maizegenetics.tassel.distance.centeredIbsPanelColumns
import net.maizegenetics.tassel.distance.centeredIbsPanelSites
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.io.BgzfReader
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.bgzfWorkerCount
import net.maizegenetics.tassel.io.openRandomAccessSource
import net.maizegenetics.tassel.io.parseVcfRegion
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.random.Random
import kotlin.time.TimeSource

/** Developer tools. These are not analyses, so they are hand-written rather than generated. */
class Dev : CliktCommand(name = "dev") {
    override fun help(context: Context): String = "Developer tools: test data and diagnostics"

    override fun run() = Unit
}

/** Base pairs between simulated sites, so a 100,000-site table spans 100 Mb of chromosome 1. */
private const val SITE_SPACING = 1_000

/** Simulated taxa are `taxon_0001`, with enough digits for the count. */
private const val TAXON_PREFIX = "taxon_"

/**
 * Writes a synthetic HapMap file. The same options always produce the same bytes, on any machine:
 * every draw comes from `kotlin.random.Random(seed)`, used in one fixed order.
 */
class Simulate : CliktCommand(name = "simulate") {
    override fun help(context: Context): String = "Generate a synthetic HapMap genotype file"

    private val taxa by option("--taxa", help = "Number of taxa", metavar = "INT")
        .int().restrictTo(min = 1).default(100)
    private val sites by option("--sites", help = "Number of sites", metavar = "INT")
        .int().restrictTo(min = 1).default(1_000)
    private val minMaf by option("--min-maf", help = "Smallest minor allele frequency", metavar = "NUMBER")
        .double().restrictTo(0.0..0.5).default(0.05)
    private val maxMaf by option("--max-maf", help = "Largest minor allele frequency", metavar = "NUMBER")
        .double().restrictTo(0.0..0.5).default(0.5)
    private val missingRate by option(
        "--missing-rate",
        help = "Proportion of genotypes left unknown",
        metavar = "NUMBER",
    ).double().restrictTo(0.0..1.0).default(0.0)
    private val seed by option("--seed", help = "Random seed; the same seed gives the same file", metavar = "INT")
        .long().default(0L)
    private val output by option(
        "--output",
        "-o",
        help = "Where to write the HapMap file; omit to write to standard output",
        metavar = "PATH",
    )

    override fun run() = cliRun {
        require(minMaf <= maxMaf) { "--min-maf ($minMaf) is above --max-maf ($maxMaf)" }
        CliLog.info("simulating $taxa taxa x $sites sites, seed $seed")
        val table = simulateGenotypes(taxa, sites, minMaf, maxMaf, missingRate, seed)
        runBlocking { writeResult(table, output, stderrProgress()) }
    }
}

/**
 * Times [net.maizegenetics.tassel.distance.centeredIbs] after the file has been read, so the number
 * is the kernel and not JVM startup or the matrix write. Not an analysis: hand-written, like [Simulate].
 */
class KinshipBench : CliktCommand(name = "kinship-bench") {
    override fun help(context: Context): String = "Time the Centered IBS kernel on one genotype file"

    private val input by option("--input", help = "HapMap or VCF to read once, then time", metavar = "PATH").required()
    private val runs by option("--runs", help = "Timed runs after one discarded warm-up", metavar = "INT")
        .int().restrictTo(min = 1).default(3)

    override fun run() = cliRun {
        runBlocking {
            val table = readGenotypes(input)
            val panelColumns = centeredIbsPanelColumns(table.numberOfTaxa)
            val panelSites = centeredIbsPanelSites(table.numberOfTaxa, table.numberOfSites)
            val seconds = DoubleArray(runs)
            repeat(runs + 1) { iteration ->
                val started = TimeSource.Monotonic.markNow()
                // The kernel inherits its caller's dispatcher. The CLI and the desktop app run on
                // Dispatchers.Default, which is the parallelism this number is meant to include.
                withContext(Dispatchers.Default) { centeredIbs(table) }
                val elapsed = started.elapsedNow().inWholeNanoseconds / 1e9
                if (iteration > 0) seconds[iteration - 1] = elapsed
            }
            seconds.sort()
            // Same rank as scripts/bench/run-benchmarks.sh: 1-based int((n + 1) / 2).
            val medianAt = (runs + 1) / 2 - 1
            echo(
                "taxa=${table.numberOfTaxa} sites=${table.numberOfSites} " +
                    "panelColumns=$panelColumns panelSites=$panelSites " +
                    "min=${formatSeconds(seconds.first())} median=${formatSeconds(seconds[medianAt])}",
            )
        }
    }
}

/**
 * Times a full BGZF inflate of one file, without parsing records.
 *
 * Not an analysis: hand-written, like [KinshipBench]. The printed median is the scan after one discarded
 * warm-up. Omit `--workers` to use [bgzfWorkerCount].
 */
class BgzfBench : CliktCommand(name = "bgzf-bench") {
    override fun help(context: Context): String = "Time BGZF inflation of one file, at 1 worker or at the host count"

    private val input by option("--input", help = "bgzipped file to inflate", metavar = "PATH").required()
    private val runs by option("--runs", help = "Timed runs after one discarded warm-up", metavar = "INT")
        .int().restrictTo(min = 1).default(3)
    private val workers by option(
        "--workers",
        help = "Inflate workers. Omit for one worker per processor.",
        metavar = "INT",
    ).int().restrictTo(min = 1)

    override fun run() = cliRun {
        val count = workers ?: bgzfWorkerCount()
        runBlocking {
            openRandomAccessSource(input).use { source ->
                val seconds = DoubleArray(runs)
                var bytes = 0L
                repeat(runs + 1) { iteration ->
                    val started = TimeSource.Monotonic.markNow()
                    bytes = scanBgzf(source, count)
                    val elapsed = started.elapsedNow().inWholeNanoseconds / 1e9
                    if (iteration > 0) seconds[iteration - 1] = elapsed
                }
                seconds.sort()
                val medianAt = (runs + 1) / 2 - 1
                echo(
                    "bytes=$bytes workers=$count min=${formatSeconds(seconds.first())} " +
                        "median=${formatSeconds(seconds[medianAt])}",
                )
            }
        }
    }
}

/**
 * Times one indexed region read.
 *
 * Not an analysis: hand-written, like [BgzfBench]. The printed median is the read after one discarded
 * warm-up, so it does not include JVM startup. `--region` uses the same syntax as `bcftools view -r`.
 */
class VcfRegionBench : CliktCommand(name = "vcf-region-bench") {
    override fun help(context: Context): String = "Time one indexed VCF region, after the header and index are open"

    private val input by option("--input", help = "bgzipped VCF", metavar = "PATH").required()
    private val indexPath by option("--index", help = "CSI or TBI index", metavar = "PATH").required()
    private val region by option(
        "--region",
        help = "One region, as bcftools view -r reads it",
        metavar = "REGION",
    ).required()
    private val runs by option("--runs", help = "Timed runs after one discarded warm-up", metavar = "INT")
        .int().restrictTo(min = 1).default(3)

    override fun run() = cliRun {
        val parsed = parseVcfRegion(region)
        runBlocking {
            openRandomAccessSource(input).use { data ->
                openRandomAccessSource(indexPath).use { index ->
                    IndexedVcf.open(data, index).use { vcf ->
                        val seconds = DoubleArray(runs)
                        var sites = 0
                        repeat(runs + 1) { iteration ->
                            val started = TimeSource.Monotonic.markNow()
                            sites = vcf.read(parsed).table.numberOfSites
                            val elapsed = started.elapsedNow().inWholeNanoseconds / 1e9
                            if (iteration > 0) seconds[iteration - 1] = elapsed
                        }
                        seconds.sort()
                        val medianAt = (runs + 1) / 2 - 1
                        echo(
                            "sites=$sites region=$region min=${formatSeconds(seconds.first())} " +
                                "median=${formatSeconds(seconds[medianAt])}",
                        )
                    }
                }
            }
        }
    }
}

private suspend fun scanBgzf(source: RandomAccessSource, workers: Int): Long {
    var total = 0L
    BgzfReader(source, workers = workers).use { reader ->
        val chunk = ByteArray(64 * 1024)
        while (true) {
            val count = reader.read(chunk, 0, chunk.size)
            if (count == 0) break
            total += count
        }
    }
    return total
}

/** Six digits after the decimal, so a tutorial-sized kernel does not round to 0.00. */
private fun formatSeconds(seconds: Double): String {
    val micros = (seconds * 1_000_000.0).toLong()
    val whole = micros / 1_000_000
    val fraction = (micros % 1_000_000).toString().padStart(6, '0')
    return "$whole.$fraction"
}

/**
 * A genotype table drawn from [seed]. Each site gets a minor allele frequency drawn uniformly from
 * `[minMaf, maxMaf]` and two distinct nucleotides; each genotype is then two independent alleles
 * (Hardy-Weinberg), or unknown with probability [missingRate].
 *
 * The draws happen site by site and, within a site, taxon by taxon, so the output depends only on
 * the arguments.
 */
fun simulateGenotypes(
    taxa: Int,
    sites: Int,
    minMaf: Double = 0.05,
    maxMaf: Double = 0.5,
    missingRate: Double = 0.0,
    seed: Long = 0L,
): GenotypeTable {
    require(taxa > 0 && sites > 0) { "taxa and sites must be positive: $taxa x $sites" }
    require(minMaf in 0.0..0.5 && maxMaf in 0.0..0.5) { "allele frequencies must be in [0, 0.5]" }
    require(minMaf <= maxMaf) { "minMaf ($minMaf) is above maxMaf ($maxMaf)" }
    require(missingRate in 0.0..1.0) { "missingRate must be in [0, 1]: $missingRate" }
    require(sites.toLong() * taxa <= Int.MAX_VALUE) { "$taxa taxa x $sites sites is too large to hold in memory" }

    val random = Random(seed)
    val chromosome = Chromosome.of(1)
    val nucleotides = byteArrayOf(Allele.A, Allele.C, Allele.G, Allele.T)
    val positions = ArrayList<Position>(sites)
    val genotypes = ByteArray(sites * taxa)
    for (site in 0 until sites) {
        val maf = if (minMaf == maxMaf) minMaf else minMaf + random.nextDouble() * (maxMaf - minMaf)
        val major = nucleotides[random.nextInt(nucleotides.size)]
        var minor = nucleotides[random.nextInt(nucleotides.size)]
        while (minor == major) minor = nucleotides[random.nextInt(nucleotides.size)]
        positions += Position(chromosome, (site + 1) * SITE_SPACING)
        val offset = site * taxa
        for (taxon in 0 until taxa) {
            genotypes[offset + taxon] = if (random.nextDouble() < missingRate) {
                Genotype.UNKNOWN
            } else {
                val first = if (random.nextDouble() < maf) minor else major
                val second = if (random.nextDouble() < maf) minor else major
                Genotype.of(first, second)
            }
        }
    }
    val names = List(taxa) { TAXON_PREFIX + (it + 1).toString().padStart(taxa.toString().length, '0') }
    return InMemoryGenotypeTable.of(TaxaList.of(names), PositionList(positions), genotypes)
}
