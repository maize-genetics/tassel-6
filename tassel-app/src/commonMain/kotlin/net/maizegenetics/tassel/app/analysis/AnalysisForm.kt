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

/**
 * What kind of control one parameter needs. This is the whole of the app's knowledge about parameter types: one
 * sealed hierarchy mapped onto another by [formField], so a composable switches on a handful of control shapes
 * instead of on `ParameterDescriptor`, and the mapping itself stays a pure function a test can walk.
 */
sealed interface FieldKind {
    /** A number typed as text, with the spec's inclusive bounds spelled out for the help line. */
    data class Number(val min: String?, val max: String?, val whole: Boolean) : FieldKind

    /** A switch. */
    data object Toggle : FieldKind

    /** One of [choices], which are the enum constant names the engine parses. */
    data class Choice(val choices: List<String>) : FieldKind

    /** Free text. */
    data object Text : FieldKind

    /** Free text holding one comma-separated list, which is how TASSEL 5 writes `-taxaList W22,B73`. */
    data object TextList : FieldKind

    /** A path. The app has no file system in the browser, so this is typed rather than browsed. */
    data class Path(val kind: FileKind) : FieldKind
}

/**
 * One row of a generated form.
 *
 * @property name parameter name, which is also the key [AnalysisForm.values] holds it under.
 * @property help the spec's one-line help.
 * @property required whether the engine refuses to run without it.
 * @property default the spec's default as the text the engine parses, or `""` when it has none.
 * @property kind the control to draw.
 */
data class FormField(
    val name: String,
    val help: String,
    val required: Boolean,
    val default: String,
    val kind: FieldKind,
)

/**
 * A form for one operation and the text the user has typed into it.
 *
 * Values are text because [net.maizegenetics.tassel.engine.RunRequest] carries text: it is the one form every front
 * end can send over any transport (ADR-0004), and it means the form never parses anything itself. `""` means "not
 * supplied", which is what [Operation.bind] already takes it for.
 *
 * @property operationId id of the operation this form runs.
 * @property title the operation's help line, or its id when it has none.
 * @property fields the parameters, in the spec's declaration order.
 * @property values current text by parameter name.
 * @property problems what the engine rejected last time, by parameter name.
 */
data class AnalysisForm(
    val operationId: String,
    val title: String,
    val fields: List<FormField>,
    val values: Map<String, String>,
    val problems: Map<String, String> = emptyMap(),
) {
    /** This form with [name] set to [value], and that field's last complaint cleared. */
    fun with(name: String, value: String): AnalysisForm =
        copy(values = values + (name to value), problems = problems - name)

    /** The parameters to send the engine: everything the user left blank is simply not supplied. */
    val parameters: Map<String, String> get() = values.filterValues { it.isNotEmpty() }
}

/**
 * The form for [descriptor], filled in with the spec's defaults.
 *
 * Nothing here names an operation or a parameter, so an operation added to `TasselOperations.registry` gets its form
 * with no change to the app (ADR-0002).
 */
fun analysisForm(descriptor: OperationDescriptor): AnalysisForm {
    val fields = descriptor.parameters.map { formField(it) }
    return AnalysisForm(
        operationId = descriptor.id,
        title = descriptor.help.ifEmpty { descriptor.id },
        fields = fields,
        values = fields.associate { it.name to it.default },
    )
}

/**
 * The control for one parameter.
 *
 * The `when` is exhaustive over the sealed [ParameterDescriptor] hierarchy and deliberately has no `else`: an eighth
 * descriptor type must be a compile error here rather than a field that silently fails to render. `tassel-cli`'s
 * `OperationCommand.registerParameter` is the same pattern for options.
 */
fun formField(parameter: ParameterDescriptor): FormField {
    val (default, kind) = when (parameter) {
        is IntParameterDescriptor ->
            parameter.default?.toString().orEmpty() to
                FieldKind.Number(parameter.min?.toString(), parameter.max?.toString(), whole = true)
        is DoubleParameterDescriptor ->
            // The spec's own text for the default, so the form and the engine agree on it to the last digit.
            parameter.default?.toString().orEmpty() to
                FieldKind.Number(parameter.min?.toString(), parameter.max?.toString(), whole = false)
        is BooleanParameterDescriptor -> parameter.default?.toString().orEmpty() to FieldKind.Toggle
        is EnumParameterDescriptor -> parameter.default.orEmpty() to FieldKind.Choice(parameter.choices)
        is FileParameterDescriptor -> parameter.default.orEmpty() to FieldKind.Path(parameter.kind)
        is StringParameterDescriptor -> parameter.default.orEmpty() to FieldKind.Text
        is StringListParameterDescriptor -> parameter.default?.joinToString(",").orEmpty() to FieldKind.TextList
    }
    return FormField(parameter.name, parameter.help, parameter.required, default, kind)
}

/**
 * This form carrying what the engine complained about.
 *
 * The app never parses or range-checks a value itself: it sends the text, and the engine's own `bind` decides. So
 * the form rejects exactly what the engine rejects, shows every problem at once, and needs no change when a spec's
 * range does. A remote engine would reach the same code by the same route.
 */
fun AnalysisForm.withProblems(invalid: InvalidParametersException): AnalysisForm =
    copy(problems = invalid.problems.associate { it.parameter to it.message })
