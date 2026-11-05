package net.maizegenetics.tassel.ops

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertSame

class OperationRegistryTest {
    private val registry = OperationRegistry(scaleValues, readLines, spin)

    @Test
    fun looksUpById() {
        assertSame(scaleValues, registry["scale-values"])
        assertNull(registry["missing"])
        assertEquals(listOf("scale-values", "read-lines", "spin"), registry.descriptors.map { it.id })
    }

    @Test
    fun looksUpByLegacyPluginAndParameterName() {
        assertSame(scaleValues, registry.byLegacyPlugin("ScaleValuesPlugin"))
        assertNull(registry.byLegacyPlugin("scalevaluesplugin"))
        assertEquals("factor", registry.parameterByLegacyName("ScaleValuesPlugin", "scaleFactor")?.name)
        assertEquals("path", registry.parameterByLegacyName("ReadLinesPlugin", "inputFile")?.name)
        assertNull(registry.parameterByLegacyName("ScaleValuesPlugin", "scalefactor"))
        assertNull(registry.parameterByLegacyName("ReadLinesPlugin", "scaleFactor"))
        assertEquals("rounding", scaleValues.parameterByLegacyName("roundingMode")?.name)
    }

    @Test
    fun rejectsDuplicates() {
        assertFailsWith<IllegalArgumentException> { OperationRegistry(listOf(spin, spin)) }
        val sameLegacy = operation<Unit, Unit>(id = "other", legacyPlugin = "ReadLinesPlugin") { run { _, _ -> } }
        assertFailsWith<IllegalArgumentException> { OperationRegistry(readLines, sameLegacy) }
    }
}
