package net.maizegenetics.tassel.io

import kotlin.test.Test
import kotlin.test.assertEquals

class BgzfWorkersWasmTest {
    @Test
    fun wasmUsesOneWorker() {
        assertEquals(1, bgzfWorkerCount())
    }
}
