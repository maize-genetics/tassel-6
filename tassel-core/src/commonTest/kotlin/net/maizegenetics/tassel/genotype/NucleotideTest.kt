package net.maizegenetics.tassel.genotype

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

/**
 * Encoding tables compared with TASSEL 5 `v5.2.98`. The expected values were printed by a small Java program calling
 * `net.maizegenetics.dna.snp.NucleotideAlignmentConstants` in `sTASSEL.jar`.
 */
class NucleotideTest {
    @Test
    fun genotypeForEveryCharacterMatchesTassel5() {
        val expected = hexBytes(DIPLOID_BY_CHAR)
        for (code in 0..255) assertEquals(expected[code], Nucleotide.genotype(code.toChar()), "char $code")
        assertFailsWith<IllegalArgumentException> { Nucleotide.genotype('\u0100') }
    }

    @Test
    fun alleleForEveryCharacterMatchesTassel5() {
        val expected = hexBytes(ALLELE_BY_CHAR)
        for (code in 0..255) assertEquals(expected[code], Nucleotide.allele(code.toChar()), "char $code")
    }

    @Test
    fun iupacForEveryGenotypeMatchesTassel5() {
        for (genotype in 0..255) {
            assertEquals(IUPAC_BY_GENOTYPE[genotype], Nucleotide.iupac(genotype.toByte()), "genotype $genotype")
        }
    }

    @Test
    fun complementsMatchTassel5() {
        val expected = hexBytes(COMPLEMENT_BY_GENOTYPE)
        for (genotype in 0..255) {
            assertEquals(expected[genotype], Nucleotide.complementGenotype(genotype.toByte()), "genotype $genotype")
        }
        for (code in 0..300) {
            val expectedCode = IUPAC_COMPLEMENTS[code.toChar()] ?: 'X'
            assertEquals(expectedCode, Nucleotide.complementIupac(code.toChar()), "char $code")
        }
    }

    @Test
    fun alleleStringsMatchTassel5() {
        assertEquals("ACGT+-XXXXXXXXZN", Nucleotide.alleleStrings.joinToString(""))
        assertEquals("A:G", Nucleotide.genotypeString(0x02))
        assertEquals("N:N", Nucleotide.genotypeString(Genotype.UNKNOWN))
        assertEquals("+:Z", Nucleotide.genotypeString(0x4E))
    }

    @Test
    fun twoCharacterCodesMatchTassel5() {
        for (entry in PAIRS.split(' ')) {
            val (code, hex) = entry.split('=')
            assertEquals(hex.toInt(16).toByte(), Nucleotide.genotype(code), code)
        }
    }

    @Test
    fun everyPrintablePairMatchesTassel5() {
        // FNV-1a over Nucleotide.genotype for every pair of printable ASCII characters, which TASSEL 5 maps through a
        // hashed table, so most invalid pairs return a valid-looking genotype.
        var hash = 0x811c9dc5.toInt()
        for (first in 0x20 until 0x7F) {
            for (second in 0x20 until 0x7F) {
                val genotype = Nucleotide.genotype("${first.toChar()}${second.toChar()}")
                hash = (hash xor (genotype.toInt() and 0xFF)) * 0x01000193
            }
        }
        assertEquals("14796a47", hash.toUInt().toString(16))
    }

    @Test
    fun iupacCodesRoundTrip() {
        for (code in "ACGT+-RYSWKMNZ0") assertEquals(code, Nucleotide.iupac(Nucleotide.genotype(code)), "$code")
        assertEquals(0xFF.toByte(), Nucleotide.genotype('X'))
        assertEquals(0x21.toByte(), Nucleotide.genotype('S'))
    }

    @Test
    fun allelePairsRoundTrip() {
        val alleles = "ACGT+-"
        for (first in alleles) {
            for (second in alleles) {
                val genotype = Nucleotide.genotype("$first$second")
                assertEquals(Genotype.of(Nucleotide.allele(first), Nucleotide.allele(second)), genotype)
                assertEquals("$first:$second", Nucleotide.genotypeString(genotype))
            }
        }
    }

    @Test
    fun rejectsLongCodes() {
        assertFailsWith<IllegalArgumentException> { Nucleotide.genotype("") }
        assertFailsWith<IllegalArgumentException> { Nucleotide.genotype("A/C") }
    }

    @Test
    fun genotypeBitOperations() {
        assertEquals(0x23.toByte(), Genotype.of(Allele.G, Allele.T))
        assertEquals(Allele.UNKNOWN, Genotype.firstAllele(Genotype.UNKNOWN))
        assertEquals(Allele.RARE, Genotype.secondAllele(Genotype.RARE))
        assertEquals(true, Genotype.isHeterozygous(0x0F))
        assertEquals(false, Genotype.isHeterozygous(Genotype.UNKNOWN))
        assertEquals(false, Genotype.isHomozygous(Genotype.UNKNOWN))
        assertEquals(true, Genotype.isHomozygous(Genotype.RARE))
    }

    private fun hexBytes(hex: String): ByteArray = ByteArray(hex.length / 2) {
        hex.substring(it * 2, it * 2 + 2).toInt(16).toByte()
    }

    private companion object {
        const val DIPLOID_BY_CHAR =
            "5466666666664466666666666666666666666666666666665566666666664566" +
                "666666f466f466e466342444045514665435256605661566666666f566f566e5" +
                "6600661166666622666666236601ff666666022133666603ff13ee6666666666" +
                "6666666666666666666666666666666666666666666666666666666666666666" +
                "663020f100f110e1663121f001f011e066664ff25ff266e26632226602661266" +
                "5166666666664166506666666666406652ff66ef666642ff0f661f66663f2f66" +
                "66664e665e66666666664f665f666666036613666633236666f366e3666666f3" +
                "0ffe1fee663f2ffe0eff1eef663e2eff66666666666666666666436653666666"

        const val ALLELE_BY_CHAR =
            "0606060606060606060606060606060606060606060606060606060606060606" +
                "0606060606060606060606040605060606060606060606060606060606060606" +
                "06000601060606020606060606060f0606060606030606060606060606060606" +
                "06000601060606020606060606060f0606060606030606060606060606060606" +
                "0606060606060606060606060606060606060606060606060606060606060606" +
                "0606060606060606060606060606060606060606060606060606060606060606" +
                "0606060606060606060606060606060606060606060606060606060606060606" +
                "0606060606060606060606060606060606060606060606060606060606060606"

        const val IUPAC_BY_GENOTYPE =
            "AMRW00XXXXXXXXAAMCSY00XXXXXXXXCCRSGK00XXXXXXXXGGWYKT00XXXXXXXXTT" +
                "0000+0XXXXXXXX++00000-XXXXXXXX--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" +
                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" +
                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXACGT+-XXXXXXXXZNACGT+-XXXXXXXXNN"

        const val COMPLEMENT_BY_GENOTYPE =
            "333231303435363738393a3b3c3d3e3f232221202425262728292a2b2c2d2e2f" +
                "131211101415161718191a1b1c1d1e1f030201000405060708090a0b0c0d0e0f" +
                "434241404445464748494a4b4c4d4e4f535251505455565758595a5b5c5d5e5f" +
                "636261606465666768696a6b6c6d6e6f737271707475767778797a7b7c7d7e7f" +
                "838281808485868788898a8b8c8d8e8f939291909495969798999a9b9c9d9e9f" +
                "a3a2a1a0a4a5a6a7a8a9aaabacadaeafb3b2b1b0b4b5b6b7b8b9babbbcbdbebf" +
                "c3c2c1c0c4c5c6c7c8c9cacbcccdcecfd3d2d1d0d4d5d6d7d8d9dadbdcdddedf" +
                "e3e2e1e0e4e5e6e7e8e9eaebecedeeeff3f2f1f0f4f5f6f7f8f9fafbfcfdfeff"

        /** Every character with an IUPAC complement other than `X`. */
        val IUPAC_COMPLEMENTS = mapOf(
            '+' to '+', '-' to '-', '0' to '0', 'A' to 'T', 'C' to 'G', 'G' to 'C', 'K' to 'M', 'M' to 'K',
            'N' to 'N', 'R' to 'Y', 'S' to 'S', 'T' to 'A', 'W' to 'W', 'Y' to 'R', 'Z' to 'Z',
        )

        const val PAIRS =
            "AA=00 AC=01 AG=02 AT=03 A+=04 A-=05 AN=0f AX=0f AZ=0e CA=10 CC=11 CG=12 CT=13 C+=14 C-=15 CN=1f CX=1f " +
                "CZ=1e GA=20 GC=21 GG=22 GT=23 G+=24 G-=25 GN=2f GX=2f GZ=2e TA=30 TC=31 TG=32 TT=33 T+=34 T-=35 " +
                "TN=3f TX=3f TZ=3e +A=40 +C=41 +G=42 +T=43 ++=44 +-=45 +N=4f +X=4f +Z=4e -A=50 -C=51 -G=52 -T=53 " +
                "-+=54 --=55 -N=5f -X=5f -Z=5e NA=f0 NC=f1 NG=f2 NT=f3 N+=f4 N-=f5 NN=ff NX=ff NZ=fe XA=f0 XC=f1 " +
                "XG=f2 XT=f3 X+=f4 X-=f5 XN=ff XX=ff XZ=fe ZA=e0 ZC=e1 ZG=e2 ZT=e3 Z+=e4 Z-=e5 ZN=ef ZX=ef ZZ=ee"
    }
}
