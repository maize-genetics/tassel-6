package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlinx.coroutines.runBlocking
import net.maizegenetics.tassel.filter.SiteFilter
import net.maizegenetics.tassel.filter.filterSites
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.openRandomAccessSource
import net.maizegenetics.tassel.io.parseVcfRegion
import java.nio.file.Files
import java.nio.file.Path
import kotlin.io.path.exists
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * `--region` and `--index` on a generated genotype command.
 *
 * `filter-sites` with its defaults keeps every site, so the HapMap it writes is `HapMap.write` of
 * [IndexedVcf.read]. The expected text is that write, not a second oracle file.
 */
class IndexedVcfCliTest {
    private val tutorialGz = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
    private val tutorialCsi = "$tutorialGz.csi"
    private val tutorialTbi = "$tutorialGz.tbi"
    private val wheatGz = testDataPath("htslib/wheat/wheat.vcf.gz")
    private val wheatCsi = "$wheatGz.csi"
    private val largeGz = testDataPath("htslib/large-pos/large_pos.vcf.gz")

    @Test
    fun filterSitesWritesTheIndexedRead() {
        assertRegionHapMap(tutorialGz, "1:157104-10068726", tutorialCsi)
        assertRegionHapMap(wheatGz, "chr3B:599000000-601000000", wheatCsi)
        assertRegionHapMap(wheatGz, "chr3B", wheatCsi)
    }

    @Test
    fun tutorialCsiAndTbiWriteTheSameHapMap() {
        val region = "1:157104-10068726"
        val fromCsi = filterSitesHapMap(tutorialGz, region, index = tutorialCsi)
        val fromTbi = filterSitesHapMap(tutorialGz, region, index = tutorialTbi)
        assertEquals(hapMapOf(regionTable(tutorialGz, tutorialCsi, region)), fromCsi)
        assertEquals(fromCsi, fromTbi)
    }

    @Test
    fun anIndexWithAnUnrelatedNameWorks() {
        val directory = tempDirectory("tassel-region-alias")
        val alias = directory.resolve("notes.txt")
        Files.copy(Path.of(tutorialCsi), alias)
        val region = "1:157104-157104"
        val written = filterSitesHapMap(tutorialGz, region, index = alias.toString())
        assertEquals(hapMapOf(regionTable(tutorialGz, tutorialCsi, region)), written)
    }

    @Test
    fun lookupUsesCsiBeforeTbi() {
        val region = "1:157104-157104"
        val foreign = planted("tassel-region-csi-first", csi = wheatCsi, tbi = tutorialTbi)
        val ignoredTbi = tasselCommand().test(
            "--quiet filter-sites --input ${foreign.vcf} --region $region --output ${foreign.output}",
        )
        assertEquals(EXIT_USER_ERROR, ignoredTbi.statusCode, ignoredTbi.output)
        assertContains(ignoredTbi.output, "not in the index")
        assertContains(ignoredTbi.output, "chr1A")

        val onlyTbi = planted("tassel-region-tbi-only", csi = null, tbi = tutorialTbi)
        val fromTbi = filterSitesHapMap(onlyTbi.vcf.toString(), region)
        assertEquals(hapMapOf(regionTable(tutorialGz, tutorialTbi, region)), fromTbi)

        val garbageTbi = planted("tassel-region-garbage-tbi", csi = tutorialCsi, tbiText = "not an index")
        val fromCsi = filterSitesHapMap(garbageTbi.vcf.toString(), region)
        assertEquals(hapMapOf(regionTable(tutorialGz, tutorialCsi, region)), fromCsi)
    }

    @Test
    fun aCorruptCsiIsNotReplacedByTheTbi() {
        val planted = planted("tassel-region-corrupt-csi", csiText = "not an index", tbi = tutorialTbi)
        val result = tasselCommand().test(
            "--quiet filter-sites --input ${planted.vcf} --region 1:157104-157104 --output ${planted.output}",
        )
        assertUserError(result.statusCode, result.output, "cannot read index", planted.vcf.toString() + ".csi")
    }

    @Test
    fun regionCombinesWithMinMaf() {
        val region = "1:157104-10068726"
        val table = regionTable(tutorialGz, tutorialCsi, region)
        val filtered = runBlocking { filterSites(table, SiteFilter(minMaf = 0.05)) }
        assertTrue(filtered.numberOfSites < table.numberOfSites, "min-maf kept every site")
        val directory = tempDirectory("tassel-region-maf")
        val base = directory.resolve("kept")
        val result = tasselCommand().test(
            "--quiet filter-sites --input $tutorialGz --region $region --min-maf 0.05 --output $base",
        )
        assertEquals(0, result.statusCode, result.output)
        assertEquals(hapMapOf(filtered), outputText(directory.resolve("kept.hmp.txt")))
    }

    @Test
    fun anEmptyRangeAndARegionPastTheIndexLimitSucceed() {
        assertRegionHapMap(wheatGz, "chr3B:800000000-851934019", wheatCsi)
        val directory = tempDirectory("tassel-region-past-limit")
        val base = directory.resolve("past")
        val result = tasselCommand().test(
            "--quiet filter-sites --input $tutorialGz --region 1:536870913-536871000 --index $tutorialTbi --output $base",
        )
        assertEquals(0, result.statusCode, result.output)
        assertEquals(
            hapMapOf(regionTable(tutorialGz, tutorialTbi, "1:536870913-536871000")),
            outputText(directory.resolve("past.hmp.txt")),
        )
        assertTrue(directory.resolve("past.hmp.txt").exists())
    }

    @Test
    fun regionOnAFileThatIsNotBgzippedVcfIsAUserError() {
        assertUserError(
            "filter-sites --input ${testDataPath("tutorial/mdp_genotype.hmp.txt")} --region 1:1-1000",
            testDataPath("tutorial/mdp_genotype.hmp.txt"),
            ".vcf.gz",
        )
        assertUserError(
            "filter-sites --input ${testDataPath("tutorial/mdp_genotype_w_depth.vcf")} --region 1:1-1000",
            testDataPath("tutorial/mdp_genotype_w_depth.vcf"),
            ".vcf.gz",
        )
        val plainGzip = testDataPath("htslib/variants/plain.vcf.gz")
        assertUserError(
            "filter-sites --input $plainGzip --region 1:1-100",
            plainGzip,
            "not bgzip compressed",
        )
    }

    @Test
    fun indexWithoutRegionIsAUserError() {
        val result = tasselCommand().test(
            "filter-sites --input $tutorialGz --index /no/such/index.csi",
        )
        assertUserError(result.statusCode, result.output, "--index requires --region")
        assertTrue(!result.output.contains("Cannot open"), result.output)
    }

    @Test
    fun aMissingIndexIsAUserError() {
        val planted = planted("tassel-region-no-index")
        val result = tasselCommand().test(
            "filter-sites --input ${planted.vcf} --region 1:157104-157104",
        )
        assertUserError(
            result.statusCode,
            result.output,
            "no index found for ${planted.vcf}",
            "bcftools index -c ${planted.vcf}",
            "tabix -p vcf",
            "2^29",
        )
    }

    @Test
    fun anUnreadableOrCorruptIndexIsAUserError() {
        val missing = "/no/such/index.csi"
        val unreadable = tasselCommand().test(
            "filter-sites --input $tutorialGz --region 1:157104-157104 --index $missing",
        )
        assertUserError(unreadable.statusCode, unreadable.output, "Cannot open $missing")

        val directory = tempDirectory("tassel-region-corrupt")
        val bad = directory.resolve("bad.csi")
        Files.writeString(bad, "not an index")
        val corrupt = tasselCommand().test(
            "filter-sites --input $tutorialGz --region 1:157104-157104 --index $bad",
        )
        assertUserError(corrupt.statusCode, corrupt.output, "cannot read index $bad")
    }

    @Test
    fun anIndexThatDoesNotMatchTheHeaderIsAUserError() {
        val result = tasselCommand().test(
            "filter-sites --input $wheatGz --region chr3B --index $tutorialCsi",
        )
        assertUserError(result.statusCode, result.output, "does not match", wheatGz, tutorialCsi, "10")
    }

    @Test
    fun regionErrorsAreUserErrors() {
        assertUserError(
            "filter-sites --input $tutorialGz --region chr3B:",
            "chr3B:",
            "malformed",
        )
        assertUserError(
            "filter-sites --input $tutorialGz --region 1:10068726-157104",
            "1:10068726-157104",
            "inverted",
        )
        assertUserError(
            "filter-sites --input $tutorialGz --region chr3B:1-2,chr3B:3-4",
            "chr3B:1-2,chr3B:3-4",
            "comma",
        )
        val huge = "chr3B:" + "9".repeat(40)
        assertUserError(
            "filter-sites --input $tutorialGz --region $huge",
            huge,
            "does not fit in a Long",
        )
        assertUserError(
            "filter-sites --input $tutorialGz --region 99:1-1000",
            "99",
            "not in the index",
        )
        assertUserError(
            "filter-sites --input $largeGz --region chr1:3000000000-3000000000",
            "2147483647",
            "3000000000",
            "chr1",
        )
    }

    private fun assertRegionHapMap(vcf: String, region: String, index: String) {
        val written = filterSitesHapMap(vcf, region)
        assertEquals(hapMapOf(regionTable(vcf, index, region)), written, region)
    }

    private fun filterSitesHapMap(vcf: String, region: String, index: String? = null): String {
        val directory = tempDirectory("tassel-region")
        val base = directory.resolve("out")
        val indexArg = index?.let { "--index $it" }.orEmpty()
        val result = tasselCommand().test(
            "--quiet filter-sites --input $vcf --region $region $indexArg --output $base",
        )
        assertEquals(0, result.statusCode, result.output)
        val written = directory.resolve("out.hmp.txt")
        assertTrue(written.exists(), "filter-sites wrote no out.hmp.txt for $region")
        return outputText(written)
    }

    private fun regionTable(vcf: String, index: String, region: String): GenotypeTable = runBlocking {
        openRandomAccessSource(vcf).use { data ->
            openRandomAccessSource(index).use { indexSource ->
                IndexedVcf.open(data, indexSource).use { indexed ->
                    indexed.read(parseVcfRegion(region)).table
                }
            }
        }
    }

    private fun hapMapOf(table: GenotypeTable): String {
        val path = tempDirectory("tassel-region-expected").resolve("expected.hmp.txt")
        writeHapMap(table, path)
        return outputText(path)
    }

    private fun planted(
        name: String,
        csi: String? = null,
        tbi: String? = null,
        csiText: String? = null,
        tbiText: String? = null,
    ): Planted {
        val directory = tempDirectory(name)
        val vcf = directory.resolve("sample.vcf.gz")
        Files.copy(Path.of(tutorialGz), vcf)
        if (csi != null) Files.copy(Path.of(csi), directory.resolve("sample.vcf.gz.csi"))
        if (csiText != null) Files.writeString(directory.resolve("sample.vcf.gz.csi"), csiText)
        if (tbi != null) Files.copy(Path.of(tbi), directory.resolve("sample.vcf.gz.tbi"))
        if (tbiText != null) Files.writeString(directory.resolve("sample.vcf.gz.tbi"), tbiText)
        return Planted(vcf, directory.resolve("out"))
    }

    private fun assertUserError(args: String, vararg needles: String) {
        val result = tasselCommand().test(args)
        assertUserError(result.statusCode, result.output, *needles)
    }

    private fun assertUserError(status: Int, output: String, vararg needles: String) {
        assertEquals(EXIT_USER_ERROR, status, output)
        assertTrue(status != EXIT_INTERNAL_ERROR, output)
        for (needle in needles) assertContains(output, needle)
    }

    private class Planted(val vcf: Path, val output: Path)
}
