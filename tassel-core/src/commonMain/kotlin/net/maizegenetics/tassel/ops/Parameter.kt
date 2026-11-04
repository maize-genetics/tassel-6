package net.maizegenetics.tassel.ops

import kotlinx.serialization.SerialName
import kotlinx.serialization.Serializable
import kotlin.math.floor
import kotlin.reflect.KProperty

/**
 * One typed parameter of an [Operation], created by declaring a delegated property in an
 * [OperationScope]. The property name becomes [name].
 *
 * A parameter is either required ([required] is true, no [default]), optional with a [default], or
 * optional without a default, in which case its value is `null` when not supplied. TASSEL 5 forbids a
 * parameter that is both required and has a default; so does this type.
 *
 * Inside an operation's `run` block, read the bound value with `parameter()` (see [RunScope]).
 *
 * @param T the Kotlin value type: `Int`, `Double`, `Boolean`, an enum, or `String` (file paths),
 *   nullable for optional parameters without a default.
 * @property name parameter name, the Kotlin property name (lower camel case).
 * @property help one-line help text shown by the CLI, UI forms, and MCP tool schemas.
 * @property legacy TASSEL 5 command-line name (`PluginParameter.cmdLineName()`), or `null` if the
 *   parameter has no TASSEL 5 equivalent.
 * @property required whether a value must be supplied.
 * @property default value used when none is supplied; `null` for required parameters and for
 *   optional parameters without a default.
 * @see net.maizegenetics.plugindef.PluginParameter
 */
class Parameter<T> internal constructor(
    internal val operationId: String,
    val name: String,
    val help: String,
    val legacy: String?,
    val required: Boolean,
    val default: T?,
    internal val type: ValueType,
) {
    /** Serializable description of this parameter. */
    val descriptor: ParameterDescriptor = type.describe(this)

    /** Lets `val x by int(...)` expose the [Parameter] itself as the property value. */
    operator fun getValue(thisRef: Any?, property: KProperty<*>): Parameter<T> = this

    /**
     * Converts a supplied value to this parameter's type and checks its range. `null` and the empty
     * string mean "not supplied", as in TASSEL 5 `AbstractPlugin.convert`.
     */
    internal fun resolve(raw: Any?): Any? {
        if (raw == null || raw == "") {
            if (required) throw InvalidValue("required parameter is missing")
            return default
        }
        val value = (if (raw is String) type.parse(raw) else type.cast(raw))
            ?: throw InvalidValue("expected ${type.label}, got ${describeRaw(raw)}")
        type.rangeProblem(value)?.let { throw InvalidValue(it) }
        return value
    }

    override fun toString(): String = "Parameter($operationId.$name)"

    private fun describeRaw(raw: Any): String = if (raw is String) "\"$raw\"" else "$raw (${raw::class.simpleName})"
}

/**
 * Kind of path a file parameter holds. Front ends use it to pick an open or save dialog and to decide
 * whether the path must already exist.
 *
 * @see net.maizegenetics.plugindef.PluginParameter.PARAMETER_TYPE
 */
@Serializable
enum class FileKind {
    /** An existing file to read (TASSEL 5 `IN_FILE`). */
    @SerialName("input")
    Input,

    /** A file to write (TASSEL 5 `OUT_FILE`). */
    @SerialName("output")
    Output,

    /** An existing directory to read (TASSEL 5 `IN_DIR`). */
    @SerialName("input-directory")
    InputDirectory,

    /** A directory to write into (TASSEL 5 `OUT_DIR`). */
    @SerialName("output-directory")
    OutputDirectory,
}

/**
 * A parameter declared in an [OperationScope] but not yet bound to a property name. Delegating a
 * property to it (`val minMaf by double(...)`) creates and registers the [Parameter].
 */
class ParameterDeclaration<T> internal constructor(
    internal val scope: OperationScope<*, *>,
    internal val type: ValueType,
    internal val default: T?,
    internal val legacy: String?,
    internal val help: String,
    internal val isRequired: Boolean = false,
) {
    init {
        default?.let { value -> type.rangeProblem(value)?.let { throw IllegalArgumentException("default: $it") } }
    }

    /** Registers the parameter under the delegated property's name. */
    operator fun provideDelegate(thisRef: Any?, property: KProperty<*>): Parameter<T> = scope.register(
        Parameter(scope.id, property.name, help, legacy, isRequired, default, type),
    )
}

/**
 * Marks an optional parameter without a default as required: binding fails unless a value is
 * supplied, and the value read in `run` is non-null.
 *
 * ```kotlin
 * val input by file(kind = FileKind.Input, legacy = "inputFile", help = "Genotype file").required()
 * ```
 *
 * @see net.maizegenetics.plugindef.PluginParameter.Builder.required
 */
fun <T : Any> ParameterDeclaration<T?>.required(): ParameterDeclaration<T> {
    require(default == null) { "a parameter with a default cannot be required" }
    return ParameterDeclaration(scope, type, null, legacy, help, isRequired = true)
}

/**
 * Validated parameter values for one run of the operation [operationId], produced by
 * [Operation.bind]. Every parameter of the operation has an entry: the supplied value, the default,
 * or `null`.
 */
class ParameterValues internal constructor(val operationId: String, private val values: Map<String, Any?>) {
    /** Value of [parameter], typed. Throws [IllegalArgumentException] for another operation's parameter. */
    @Suppress("UNCHECKED_CAST")
    operator fun <T> get(parameter: Parameter<T>): T {
        require(parameter.operationId == operationId) { "$parameter does not belong to '$operationId'" }
        return values[parameter.name] as T
    }

    /** Value of the parameter called [name]. Throws [IllegalArgumentException] for an unknown name. */
    operator fun get(name: String): Any? {
        require(name in values) { "'$operationId' has no parameter '$name'" }
        return values[name]
    }

    /** All values by parameter name, in declaration order. */
    fun toMap(): Map<String, Any?> = values

    override fun equals(other: Any?): Boolean =
        other is ParameterValues && other.operationId == operationId && other.values == values

    override fun hashCode(): Int = 31 * operationId.hashCode() + values.hashCode()

    override fun toString(): String = "ParameterValues($operationId, $values)"
}

/** One reason a supplied parameter value was rejected. */
@Serializable
data class ParameterProblem(
    /** Name of the parameter as supplied (may be unknown to the operation). */
    val parameter: String,
    /** Human-readable reason. */
    val message: String,
) {
    override fun toString(): String = "$parameter: $message"
}

/**
 * Thrown by [Operation.bind] when supplied values are invalid. Lists every problem found, not just the
 * first, so a form or an agent can fix them all at once.
 */
class InvalidParametersException(val operationId: String, val problems: List<ParameterProblem>) :
    IllegalArgumentException("Invalid parameters for '$operationId': ${problems.joinToString("; ")}")

internal class InvalidValue(message: String) : Exception(message)

/** How one parameter type parses, type-checks, range-checks, and describes values. */
internal abstract class ValueType(val label: String) {
    /** Parses command-line text; `null` if it is not a valid value. */
    abstract fun parse(text: String): Any?

    /** Accepts an already-typed value; `null` if it has the wrong type. */
    abstract fun cast(value: Any): Any?

    /** A message if [value] is outside the allowed range, else `null`. */
    open fun rangeProblem(value: Any): String? = null

    abstract fun describe(parameter: Parameter<*>): ParameterDescriptor
}

private val INTEGER = Regex("""[+-]?\d+""")
private val DECIMAL = Regex("""[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?""")
private val INT_RANGE = Int.MIN_VALUE.toDouble()..Int.MAX_VALUE.toDouble()

internal class IntType(private val range: IntRange?) : ValueType("an integer") {
    init {
        require(range == null || !range.isEmpty()) { "empty range $range" }
    }

    // TASSEL 5 parses through BigDecimal.intValueExact, so "5.0" and "1e2" are integers but "5.5" is not.
    override fun parse(text: String): Any? = when {
        INTEGER.matches(text) -> text.toIntOrNull()
        DECIMAL.matches(text) -> text.toDouble().takeIf { it == floor(it) && it in INT_RANGE }?.toInt()
        else -> null
    }

    override fun cast(value: Any): Any? = value as? Int

    override fun rangeProblem(value: Any): String? =
        if (range == null || value as Int in range) null else "$value is outside range [${range.first}..${range.last}]"

    override fun describe(parameter: Parameter<*>) = with(parameter) {
        IntParameterDescriptor(name, help, legacy, required, default as Int?, range?.first, range?.last)
    }
}

internal class DoubleType(private val range: ClosedFloatingPointRange<Double>?) : ValueType("a number") {
    init {
        require(range == null || !range.isEmpty()) { "empty range $range" }
    }

    override fun parse(text: String): Any? = if (DECIMAL.matches(text)) text.toDouble() else null

    override fun cast(value: Any): Any? = when (value) {
        is Double -> value
        is Int -> value.toDouble()
        else -> null
    }

    override fun rangeProblem(value: Any): String? = if (range == null || value as Double in range) {
        null
    } else {
        "$value is outside range [${range.start}..${range.endInclusive}]"
    }

    override fun describe(parameter: Parameter<*>) = with(parameter) {
        DoubleParameterDescriptor(name, help, legacy, required, default as Double?, range?.start, range?.endInclusive)
    }
}

internal object BooleanType : ValueType("true or false") {
    override fun parse(text: String): Any? = text.lowercase().toBooleanStrictOrNull()

    override fun cast(value: Any): Any? = value as? Boolean

    override fun describe(parameter: Parameter<*>) = with(parameter) {
        BooleanParameterDescriptor(name, help, legacy, required, default as Boolean?)
    }
}

internal class EnumType(private val entries: List<Enum<*>>) :
    ValueType("one of ${entries.joinToString(", ") { it.name }}") {
    override fun parse(text: String): Any? = entries.firstOrNull { it.name == text }

    override fun cast(value: Any): Any? = value.takeIf { it in entries }

    override fun describe(parameter: Parameter<*>) = with(parameter) {
        EnumParameterDescriptor(name, help, legacy, required, (default as Enum<*>?)?.name, entries.map { it.name })
    }
}

internal class FileType(private val kind: FileKind) : ValueType("a path") {
    override fun parse(text: String): Any = text

    override fun cast(value: Any): Any? = value as? String

    override fun describe(parameter: Parameter<*>) = with(parameter) {
        FileParameterDescriptor(name, help, legacy, required, default as String?, kind)
    }
}
