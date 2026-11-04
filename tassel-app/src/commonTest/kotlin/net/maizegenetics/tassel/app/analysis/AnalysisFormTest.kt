package net.maizegenetics.tassel.app.analysis

import net.maizegenetics.tassel.ops.BooleanParameterDescriptor
import net.maizegenetics.tassel.ops.DoubleParameterDescriptor
import net.maizegenetics.tassel.ops.EnumParameterDescriptor
import net.maizegenetics.tassel.ops.FileKind
import net.maizegenetics.tassel.ops.FileParameterDescriptor
import net.maizegenetics.tassel.ops.IntParameterDescriptor
import net.maizegenetics.tassel.ops.InvalidParametersException
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.ParameterDescriptor
import net.maizegenetics.tassel.ops.StringListParameterDescriptor
import net.maizegenetics.tassel.ops.StringParameterDescriptor
import net.maizegenetics.tassel.ops.TasselOperations
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertTrue

/**
 * The form generator: one control per parameter, derived from the spec alone.
 *
 * [everyDescriptorType] lists all seven subtypes of the sealed `ParameterDescriptor` hierarchy. It is the test
 * that fails first if an eighth is added without teaching [formField] about it; [formField]'s own `when` has no
 * `else`, so the compiler fails before this does.
 */
class AnalysisFormTest {
    private val everyDescriptorType: List<ParameterDescriptor> = listOf(
        IntParameterDescriptor("count", help = "How many", default = 7, min = 1, max = 10),
        DoubleParameterDescriptor("minMaf", help = "Minimum frequency", default = 0.05, min = 0.0, max = 0.5),
        BooleanParameterDescriptor("covariance", help = "Use covariance", default = true),
        EnumParameterDescriptor("method", help = "Method", default = "Centered_IBS", choices = listOf("Centered_IBS")),
        FileParameterDescriptor("pedigree", help = "Pedigree file", kind = FileKind.Input),
        StringParameterDescriptor("startChr", help = "First chromosome"),
        StringListParameterDescriptor("taxaList", help = "Taxa to keep", default = listOf("B73", "Mo17")),
    )

    @Test
    fun everySealedDescriptorTypeGetsAControl() {
        val kinds = everyDescriptorType.map { formField(it).kind }

        assertEquals(
            listOf<FieldKind>(
                FieldKind.Number(min = "1", max = "10", whole = true),
                FieldKind.Number(min = "0.0", max = "0.5", whole = false),
                FieldKind.Toggle,
                FieldKind.Choice(listOf("Centered_IBS")),
                FieldKind.Path(FileKind.Input),
                FieldKind.Text,
                FieldKind.TextList,
            ),
            kinds,
        )
    }

    @Test
    fun everyDefaultIsTheTextTheEngineParses() {
        val defaults = everyDescriptorType.associate { it.name to formField(it).default }

        assertEquals("7", defaults["count"])
        assertEquals("0.05", defaults["minMaf"])
        assertEquals("true", defaults["covariance"])
        assertEquals("Centered_IBS", defaults["method"])
        // A parameter with no default is simply not supplied, which is what the empty string means.
        assertEquals("", defaults["pedigree"])
        assertEquals("", defaults["startChr"])
        // TASSEL 5 writes a list as one comma-separated value: -taxaList B73,Mo17.
        assertEquals("B73,Mo17", defaults["taxaList"])
    }

    @Test
    fun aFormIsGeneratedFromADescriptorAndNothingElse() {
        val form = analysisForm(
            OperationDescriptor(
                id = "demo",
                help = "A demonstration",
                inputType = "GenotypeTable",
                outputType = "ResultTable",
                parameters = everyDescriptorType,
            ),
        )

        assertEquals("demo", form.operationId)
        assertEquals("A demonstration", form.title)
        assertEquals(everyDescriptorType.map { it.name }, form.fields.map { it.name })
        assertEquals(everyDescriptorType.map { it.help }, form.fields.map { it.help })
    }

    @Test
    fun blankFieldsAreNotSupplied() {
        val form = analysisForm(descriptorOf(everyDescriptorType))

        assertEquals(setOf("count", "minMaf", "covariance", "method", "taxaList"), form.parameters.keys)
        assertEquals("7", form.parameters["count"])
    }

    @Test
    fun typingInAFieldReplacesItsValueAndClearsItsLastComplaint() {
        val form = analysisForm(descriptorOf(everyDescriptorType))
            .copy(problems = mapOf("count" to "too big", "minMaf" to "too small"))

        val edited = form.with("count", "3")

        assertEquals("3", edited.values["count"])
        assertEquals(mapOf("minMaf" to "too small"), edited.problems)
        // The original is untouched; the state is immutable (.cursor/rules/compose-ui.mdc).
        assertEquals("7", form.values["count"])
    }

    @Test
    fun theEnginesOwnComplaintsLandOnTheFieldsTheyBelongTo() {
        val kinship = TasselOperations.registry["kinship"]!!
        val form = analysisForm(kinship.descriptor).with("maxAlleles", "99")

        val invalid = runCatching { kinship.bind(form.parameters) }.exceptionOrNull()
        assertIs<InvalidParametersException>(invalid)
        val problems = form.withProblems(invalid).problems

        assertEquals(setOf("maxAlleles"), problems.keys)
        assertTrue(problems.getValue("maxAlleles").isNotEmpty())
    }

    @Test
    fun everyRegisteredOperationGetsAFormThatTheOperationItselfAccepts() {
        // No operation is named here: whatever the registry holds gets a form, which is the point (ADR-0002).
        for (operation in TasselOperations.registry.operations) {
            val form = analysisForm(operation.descriptor)

            assertEquals(operation.parameters.map { it.name }, form.fields.map { it.name }, operation.id)
            // The defaults the form starts with are values the operation binds without complaint.
            operation.bind(form.parameters)
        }
    }

    private fun descriptorOf(parameters: List<ParameterDescriptor>) = OperationDescriptor(
        id = "demo",
        inputType = "GenotypeTable",
        outputType = "ResultTable",
        parameters = parameters,
    )
}
