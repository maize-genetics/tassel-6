package net.maizegenetics.tassel.taxa

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class TaxaListTest {
    @Test
    fun looksUpTaxaByName() {
        val taxa = TaxaList.of("33-16", "38-11", "B73")
        assertEquals(3, taxa.size)
        assertEquals(2, taxa.indexOf("B73"))
        assertEquals(1, taxa.indexOf(Taxon("38-11")))
        assertEquals(-1, taxa.indexOf("Mo17"))
        assertTrue(Taxon("33-16") in taxa)
        assertFalse(Taxon("Mo17") in taxa)
        assertEquals(listOf("33-16", "38-11", "B73"), taxa.names)
        assertEquals(TaxaList.of("33-16", "38-11", "B73"), taxa)
    }

    @Test
    fun rejectsDuplicateNames() {
        assertFailsWith<IllegalArgumentException> { TaxaList.of("B73", "Mo17", "B73") }
    }
}
