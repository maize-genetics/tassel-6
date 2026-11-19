package net.maizegenetics.tassel.app.export

import kotlin.test.Test
import kotlin.test.assertEquals

/** The file name each genotype format produces. The writer is chosen separately, in [ExportTest]. */
class GenotypeExportTest {
    @Test
    fun theFileNameFollowsTheChosenFormat() {
        assertEquals("mdp_genotype.hmp.txt", genotypeFileName("mdp_genotype.hmp.txt"))
        assertEquals("filtered.hmp.txt", genotypeFileName("filtered"))
        assertEquals(
            "tutorial.hmp.txt + filter-sites.hmp.txt",
            genotypeFileName("tutorial.hmp.txt + filter-sites"),
        )
        assertEquals("mdp_genotype.vcf", genotypeFileName("mdp_genotype.hmp.txt", GenotypeExportFormat.Vcf))
        assertEquals("wheat.bcf", genotypeFileName("wheat.VCF.GZ", GenotypeExportFormat.Bcf))
        assertEquals("mdp_genotype.vcf", genotypeFileName("mdp_genotype.vcf.gz", GenotypeExportFormat.Vcf))
    }
}
