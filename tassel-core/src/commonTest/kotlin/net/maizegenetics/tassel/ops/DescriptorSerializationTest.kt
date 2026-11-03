package net.maizegenetics.tassel.ops

import kotlinx.serialization.json.Json
import kotlinx.serialization.json.jsonArray
import kotlinx.serialization.json.jsonObject
import kotlinx.serialization.json.jsonPrimitive
import kotlin.test.Test
import kotlin.test.assertEquals

class DescriptorSerializationTest {
    @Test
    fun descriptorDescribesEveryParameter() {
        val d = scaleValues.descriptor
        assertEquals("scale-values", d.id)
        assertEquals("ScaleValuesPlugin", d.legacyPlugin)
        assertEquals("List", d.inputType)
        assertEquals(
            listOf(
                DoubleParameterDescriptor("factor", "Multiplier", "scaleFactor", default = 1.0, min = 0.0, max = 10.0),
                IntParameterDescriptor(
                    "repeats",
                    "Times to apply the factor",
                    "numRepeats",
                    default = 1,
                    min = 1,
                    max = 5,
                ),
                IntParameterDescriptor("limit", "Keep at most this many values", "maxValues", min = 0, max = 1000),
                BooleanParameterDescriptor("negate", "Negate the result", "negateValues", default = false),
                EnumParameterDescriptor(
                    "rounding",
                    "Rounding",
                    "roundingMode",
                    default = "NONE",
                    choices = listOf("NONE", "FLOOR", "CEILING"),
                ),
                FileParameterDescriptor("report", "Optional report", "reportFile", kind = FileKind.Output),
            ),
            d.parameters,
        )
        assertEquals(
            FileParameterDescriptor("path", "File to read", "inputFile", required = true),
            readLines.descriptor.parameters.single(),
        )
    }

    @Test
    fun descriptorJsonRoundTrips() {
        for (descriptor in OperationRegistry(scaleValues, readLines, spin).descriptors) {
            val json = Json.encodeToString(OperationDescriptor.serializer(), descriptor)
            assertEquals(descriptor, Json.decodeFromString(OperationDescriptor.serializer(), json), json)
        }
    }

    @Test
    fun jsonUsesTypeDiscriminatorAndKebabFileKinds() {
        val json = Json.parseToJsonElement(
            Json.encodeToString(OperationDescriptor.serializer(), scaleValues.descriptor),
        )
        val parameters = json.jsonObject.getValue("parameters").jsonArray.map { it.jsonObject }
        assertEquals(
            listOf("double", "int", "int", "boolean", "enum", "file"),
            parameters.map {
                it.getValue("type").jsonPrimitive.content
            },
        )
        assertEquals("output", parameters.last().getValue("kind").jsonPrimitive.content)
        assertEquals(
            "siteMinAlleleFreq",
            parseLegacy("""{"type":"double","name":"minMaf","legacy":"siteMinAlleleFreq"}"""),
        )
    }

    private fun parseLegacy(json: String): String? =
        Json.decodeFromString(ParameterDescriptor.serializer(), json).legacy
}
