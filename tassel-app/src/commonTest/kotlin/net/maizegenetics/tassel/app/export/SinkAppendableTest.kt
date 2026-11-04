package net.maizegenetics.tassel.app.export

import kotlinx.io.Buffer
import kotlinx.io.readString
import net.maizegenetics.tassel.table.resultTable
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/** The bridge from the core's `Appendable` writers to a `kotlinx.io.Sink`. */
class SinkAppendableTest {
    @Test
    fun everyAppendOverloadReachesTheSink() {
        val buffer = Buffer()
        val out = SinkAppendable(buffer)

        out.append('T')
        out.append("ASSEL")
        out.append(" 6 and more", 0, 2)
        out.append(null as CharSequence?)
        out.flush()

        assertEquals("TASSEL 6null", buffer.readString())
    }

    @Test
    fun charactersAreBatchedRatherThanWrittenOneAtATime() {
        val buffer = Buffer()
        val out = SinkAppendable(buffer, chunkSize = 8)

        repeat(16) { out.append('x') }
        out.flush()

        // Two full chunks and nothing left over: 16 characters cost two writes, not sixteen.
        assertEquals(2, out.writes)
        assertEquals(16, buffer.readString().length)
    }

    @Test
    fun aLongStringCostsOneWriteRatherThanOnePerCharacter() {
        val buffer = Buffer()
        val out = SinkAppendable(buffer, chunkSize = 8)

        out.append("a".repeat(1000))
        out.flush()

        assertEquals(1, out.writes)
        assertEquals(1000, buffer.readString().length)
    }

    @Test
    fun nothingBufferedMeansNothingWritten() {
        val buffer = Buffer()
        val out = SinkAppendable(buffer)

        out.flush()
        out.flush()

        assertEquals(0, out.writes)
        assertTrue(buffer.exhausted())
    }

    @Test
    fun whatReachesTheSinkIsWhatTheCoreWriterWouldHaveProduced() {
        val table = resultTable("Overall Summary", listOf("Stat Type", "Value")) {
            row("Number of Taxa", 281.0)
            row("Number of Sites", 3093.0)
            row("Missing", null)
        }
        val buffer = Buffer()

        buffer.appendText { table.writeDelimited(it) }

        assertEquals(table.toDelimitedString(), buffer.readString())
    }

    @Test
    fun aChunkSmallerThanARowStillWritesTheSameBytes() {
        val table = resultTable("Site Summary", listOf("Site", "MAF")) {
            repeat(50) { row("S$it", it / 100.0) }
        }
        val buffer = Buffer()
        val out = SinkAppendable(buffer, chunkSize = 4)

        table.writeDelimited(out)
        out.flush()

        assertEquals(table.toDelimitedString(), buffer.readString())
    }
}
