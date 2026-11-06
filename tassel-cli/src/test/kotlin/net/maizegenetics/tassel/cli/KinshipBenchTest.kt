package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class KinshipBenchTest {
    @Test
    fun printsTheKernelShapeAndAMedian() {
        val directory = tempDirectory("tassel-kinship-bench")
        val input = directory.resolve("sim.hmp.txt")
        val simulated = tasselCommand().test(
            "--quiet dev simulate --taxa 8 --sites 20 --seed 0 --output $input",
        )
        assertEquals(0, simulated.statusCode, simulated.output)

        val result = tasselCommand().test("--quiet dev kinship-bench --input $input --runs 1")
        assertEquals(0, result.statusCode, result.output)
        assertTrue(result.output.contains("taxa=8"), result.output)
        assertTrue(result.output.contains("sites=20"), result.output)
        assertTrue(result.output.contains("panelColumns="), result.output)
        assertTrue(result.output.contains("median="), result.output)
    }
}
