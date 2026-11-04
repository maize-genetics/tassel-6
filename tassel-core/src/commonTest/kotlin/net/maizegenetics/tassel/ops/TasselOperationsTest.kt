package net.maizegenetics.tassel.ops

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

class TasselOperationsTest {
    @Test
    fun theRegistryIsBuilt() {
        val registry = TasselOperations.registry
        assertEquals(registry.operations.size, registry.descriptors.size)
        assertNull(registry["not-an-operation"])
        assertNull(registry.byLegacyPlugin("NoSuchPlugin"))
    }

    /** Guards the one line each later operation PR adds to [TasselOperations.registry]. */
    @Test
    fun everyOperationIsReachableByItsIdAndItsTassel5PluginName() {
        val registry = TasselOperations.registry
        for (operation in registry.operations) {
            assertEquals(operation, registry[operation.id], "operation '${operation.id}' is not registered by id")
            val plugin = operation.legacyPlugin
            if (plugin != null) {
                assertEquals(operation, registry.byLegacyPlugin(plugin), "'$plugin' maps to another operation")
            }
        }
        val ids = registry.operations.map { it.id }
        assertEquals(ids.distinct(), ids, "duplicate operation id")
    }

    @Test
    fun registeringOperationsKeepsTheLookups() {
        val registry = OperationRegistry(scaleValues)
        assertEquals(listOf("scale-values"), registry.operations.map { it.id })
        assertEquals(scaleValues, registry["scale-values"])
        assertEquals(scaleValues, registry.byLegacyPlugin("ScaleValuesPlugin"))
        assertTrue(registry.parameterByLegacyName("ScaleValuesPlugin", "scaleFactor") != null)
    }
}
