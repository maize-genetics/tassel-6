package net.maizegenetics.tassel.io.deflate

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class HuffmanTableTest {
    @Test
    fun completeTwoSymbolTreeFillsThePrimaryTable() {
        val table = table()
        table.build(intArrayOf(1, 1), 0, 2, allowEmpty = false)
        assertEquals(0, symbol(table, 0))
        assertEquals(1, symbol(table, 1))
        assertEquals(0, symbol(table, 2))
        assertEquals(1, symbol(table, 3))
    }

    @Test
    fun singleCodeIsAccepted() {
        val table = table()
        table.build(intArrayOf(0, 1), 0, 2, allowEmpty = false)
        assertEquals(1, symbol(table, 0))
        val failure = assertFailsWith<IllegalArgumentException> { symbol(table, 1) }
        assertTrue(failure.message!!.contains("invalid length symbol"), failure.message)
    }

    @Test
    fun emptyDistanceAlphabetIsAcceptedOnlyWhenAsked() {
        val table = table()
        table.build(intArrayOf(0, 0), 0, 2, allowEmpty = true)
        assertEquals(true, table.empty)
        val rejected = assertFailsWith<IllegalArgumentException> {
            table.build(intArrayOf(0, 0), 0, 2, allowEmpty = false)
        }
        assertTrue(rejected.message!!.contains("incomplete code lengths"), rejected.message)
    }

    @Test
    fun overSubscribedAndIncompleteLengthsFail() {
        val over = assertFailsWith<IllegalArgumentException> {
            table().build(intArrayOf(1, 1, 1), 0, 3, allowEmpty = false)
        }
        assertTrue(over.message!!.contains("over-subscribed code lengths"), over.message)
        val incomplete = assertFailsWith<IllegalArgumentException> {
            table().build(intArrayOf(2, 2), 0, 2, allowEmpty = false)
        }
        assertTrue(incomplete.message!!.contains("incomplete code lengths"), incomplete.message)
    }

    private fun table() = HuffmanTable(128, 7, "invalid length symbol")

    private fun symbol(table: HuffmanTable, peek: Int): Int {
        val entry = table.lookup(peek)
        if (entry and HuffmanTable.FLAG_INVALID != 0) {
            throw IllegalArgumentException(table.invalidMessage)
        }
        return entry and HuffmanTable.PAYLOAD_MASK
    }
}
