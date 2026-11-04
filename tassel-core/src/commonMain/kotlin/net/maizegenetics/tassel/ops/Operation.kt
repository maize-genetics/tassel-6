package net.maizegenetics.tassel.ops

import kotlin.reflect.KClass

/** Restricts implicit receivers in the operation DSL, so parameters cannot be declared inside `run`. */
@DslMarker
annotation class OperationDsl

/**
 * Declares an operation: its id, TASSEL 5 plugin name, typed parameters, and the code that runs it.
 * This spec is the single source for the CLI, UI forms, MCP tools, and bindings (ADR-0002).
 *
 * ```kotlin
 * val FilterSites = operation<GenotypeTable, GenotypeTable>(
 *     id = "filter-sites",
 *     legacyPlugin = "FilterSiteBuilderPlugin",
 * ) {
 *     val minMaf by double(default = 0.0, range = 0.0..0.5, legacy = "siteMinAlleleFreq",
 *         help = "Minimum minor allele frequency")
 *     val minCount by int(default = 0, legacy = "siteMinCount",
 *         help = "Minimum number of non-missing taxa")
 *     run { input, progress -> filterSites(input, minMaf(), minCount(), progress) }
 * }
 * ```
 *
 * Operations that take no input use `Unit` for [I].
 *
 * @param I input data type.
 * @param O output data type.
 * @param id stable lower kebab-case id (for example `filter-sites`); also the CLI subcommand name.
 * @param legacyPlugin TASSEL 5 plugin class simple name this operation replaces, or `null`.
 * @param help one-line description of the operation.
 * @param block declares parameters and calls `run` exactly once.
 * @throws IllegalArgumentException if the id, a parameter name, or a legacy name is invalid or
 *   duplicated, a default is outside its range, or `run` is not called exactly once.
 * @see net.maizegenetics.plugindef.AbstractPlugin
 * @see net.maizegenetics.plugindef.PluginParameter
 */
inline fun <reified I : Any, reified O : Any> operation(
    id: String,
    legacyPlugin: String? = null,
    help: String = "",
    noinline block: OperationScope<I, O>.() -> Unit,
): Operation<I, O> = buildOperation(id, legacyPlugin, help, I::class, O::class, block)

@PublishedApi
internal fun <I : Any, O : Any> buildOperation(
    id: String,
    legacyPlugin: String?,
    help: String,
    inputType: KClass<I>,
    outputType: KClass<O>,
    block: OperationScope<I, O>.() -> Unit,
): Operation<I, O> {
    require(ID.matches(id)) { "operation id must be lower kebab-case: '$id'" }
    require(legacyPlugin == null || legacyPlugin.isNotBlank()) { "$id: blank legacyPlugin" }
    val scope = OperationScope<I, O>(id).apply(block)
    val body = requireNotNull(scope.body) { "$id: run { } was not called" }
    return Operation(id, legacyPlugin, help, inputType, outputType, scope.parameters.toList(), body)
}

private val ID = Regex("[a-z][a-z0-9]*(-[a-z0-9]+)*")
private val PARAMETER_NAME = Regex("[a-z][A-Za-z0-9]*")

/**
 * Receiver of the [operation] block. Each factory returns a [ParameterDeclaration] to delegate a
 * property to; the property name becomes the parameter name.
 *
 * Factories that take a `default` declare an optional parameter with that default. Factories without
 * one declare an optional parameter whose value is `null` when not supplied (TASSEL 5 parameters
 * built with a `null` default); add [required] to make it required. Ranges are inclusive at both
 * ends, like the closed Guava ranges TASSEL 5 uses.
 *
 * @see net.maizegenetics.plugindef.PluginParameter.Builder
 */
@OperationDsl
class OperationScope<I : Any, O : Any> internal constructor(internal val id: String) {
    internal val parameters = mutableListOf<Parameter<*>>()
    internal var body: (suspend RunScope.(I, ProgressReporter) -> O)? = null

    /** Integer parameter with a default. */
    fun int(
        default: Int,
        range: IntRange? = null,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<Int> = ParameterDeclaration(this, IntType(range), default, legacy, help)

    /** Integer parameter with no default; `null` unless supplied. */
    fun int(range: IntRange? = null, legacy: String? = null, help: String = ""): ParameterDeclaration<Int?> =
        ParameterDeclaration(this, IntType(range), null, legacy, help)

    /** Double parameter with a default. */
    fun double(
        default: Double,
        range: ClosedFloatingPointRange<Double>? = null,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<Double> = ParameterDeclaration(this, DoubleType(range), default, legacy, help)

    /** Double parameter with no default; `null` unless supplied. */
    fun double(
        range: ClosedFloatingPointRange<Double>? = null,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<Double?> = ParameterDeclaration(this, DoubleType(range), null, legacy, help)

    /** Boolean parameter. Text values are `true` or `false`, in any letter case. */
    fun boolean(default: Boolean = false, legacy: String? = null, help: String = ""): ParameterDeclaration<Boolean> =
        ParameterDeclaration(this, BooleanType, default, legacy, help)

    /** Enum parameter with a default. Text values are constant names, matched exactly as in TASSEL 5. */
    inline fun <reified E : Enum<E>> enum(
        default: E,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<E> {
        @Suppress("UNCHECKED_CAST")
        return enumDeclaration(enumValues<E>().asList(), default, legacy, help) as ParameterDeclaration<E>
    }

    /** Enum parameter with no default; `null` unless supplied. */
    inline fun <reified E : Enum<E>> enum(legacy: String? = null, help: String = ""): ParameterDeclaration<E?> =
        enumDeclaration(enumValues<E>().asList(), null, legacy, help)

    @PublishedApi
    internal fun <E : Enum<E>> enumDeclaration(
        entries: List<E>,
        default: E?,
        legacy: String?,
        help: String,
    ): ParameterDeclaration<E?> = ParameterDeclaration(this, EnumType(entries), default, legacy, help)

    /** File or directory path parameter with a default. Paths are not checked for existence here. */
    fun file(
        default: String,
        kind: FileKind = FileKind.Input,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<String> = ParameterDeclaration(this, FileType(kind), default, legacy, help)

    /** File or directory path parameter with no default; `null` unless supplied. */
    fun file(
        kind: FileKind = FileKind.Input,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<String?> = ParameterDeclaration(this, FileType(kind), null, legacy, help)

    /** Text parameter with a default, for a value with no narrower type, such as a chromosome name. */
    fun string(default: String, legacy: String? = null, help: String = ""): ParameterDeclaration<String> =
        ParameterDeclaration(this, StringType, default, legacy, help)

    /** Text parameter with no default; `null` unless supplied. */
    fun string(legacy: String? = null, help: String = ""): ParameterDeclaration<String?> =
        ParameterDeclaration(this, StringType, null, legacy, help)

    /**
     * Text-list parameter with a default. As text it is one comma-separated list, as TASSEL 5 takes
     * `-taxaList W22,B73`; entries are not trimmed and none of them may be blank.
     */
    fun stringList(
        default: List<String>,
        legacy: String? = null,
        help: String = "",
    ): ParameterDeclaration<List<String>> = ParameterDeclaration(this, StringListType, default, legacy, help)

    /** Text-list parameter with no default; `null` unless supplied. */
    fun stringList(legacy: String? = null, help: String = ""): ParameterDeclaration<List<String>?> =
        ParameterDeclaration(this, StringListType, null, legacy, help)

    /**
     * Sets the code that runs the operation. Read parameter values in [block] by calling the
     * parameter (`minMaf()`). Report progress through `progress`; when run by an engine, each report
     * is also a cancellation point.
     */
    fun run(block: suspend RunScope.(input: I, progress: ProgressReporter) -> O) {
        require(body == null) { "$id: run { } called more than once" }
        body = block
    }

    internal fun <T> register(parameter: Parameter<T>): Parameter<T> {
        val name = parameter.name
        require(PARAMETER_NAME.matches(name)) { "$id: parameter name must be lower camel case: '$name'" }
        require(parameters.none { it.name == name }) { "$id: duplicate parameter '$name'" }
        parameter.legacy?.let { legacy ->
            require(legacy.isNotBlank()) { "$id.$name: blank legacy name" }
            require(parameters.none { it.legacy == legacy }) { "$id: duplicate legacy name '$legacy'" }
        }
        parameters += parameter
        return parameter
    }
}

/**
 * Receiver of an operation's `run` block. Gives typed access to the bound parameter values.
 *
 * @property values the values this run was bound with.
 */
@OperationDsl
class RunScope internal constructor(val values: ParameterValues) {
    /** The value bound to this parameter for the current run. */
    operator fun <T> Parameter<T>.invoke(): T = values[this]
}

/**
 * A declared operation, created by [operation]. Immutable and safe to share.
 *
 * @property id stable lower kebab-case id.
 * @property legacyPlugin TASSEL 5 plugin class simple name, or `null`.
 * @property help one-line description.
 * @property inputType class of the input data (`Unit` when the operation takes none).
 * @property outputType class of the output data.
 * @property parameters parameters in declaration order.
 * @see net.maizegenetics.plugindef.AbstractPlugin
 */
class Operation<I : Any, O : Any> internal constructor(
    val id: String,
    val legacyPlugin: String?,
    val help: String,
    val inputType: KClass<I>,
    val outputType: KClass<O>,
    val parameters: List<Parameter<*>>,
    private val body: suspend RunScope.(I, ProgressReporter) -> O,
) {
    /** Serializable description of this operation and its parameters. */
    val descriptor: OperationDescriptor = OperationDescriptor(
        id = id,
        legacyPlugin = legacyPlugin,
        help = help,
        inputType = inputType.simpleName ?: "?",
        outputType = outputType.simpleName ?: "?",
        parameters = parameters.map { it.descriptor },
    )

    /** The parameter called [name], or `null`. */
    fun parameter(name: String): Parameter<*>? = parameters.firstOrNull { it.name == name }

    /**
     * The parameter whose TASSEL 5 command-line name is [legacyName], or `null`. Matching is
     * case-sensitive, as in TASSEL 5 `AbstractPlugin.setParameter`.
     */
    fun parameterByLegacyName(legacyName: String): Parameter<*>? = parameters.firstOrNull { it.legacy == legacyName }

    /**
     * Validates [values] and fills in defaults. Keys are parameter names. Each value is either the
     * parameter's Kotlin type (an `Int` is also accepted for a double) or command-line text, which is
     * parsed. `null` or `""` means "not supplied".
     *
     * @throws InvalidParametersException listing unknown names, unparsable or mistyped values,
     *   out-of-range values, and missing required parameters.
     * @see net.maizegenetics.plugindef.AbstractPlugin.setParameter
     */
    fun bind(values: Map<String, Any?> = emptyMap()): ParameterValues {
        val problems = mutableListOf<ParameterProblem>()
        for (key in values.keys) {
            if (parameter(key) != null) continue
            val hint = parameterByLegacyName(key)?.let { "; it is the TASSEL 5 name of '${it.name}'" } ?: ""
            problems += ParameterProblem(key, "unknown parameter$hint")
        }
        val bound = LinkedHashMap<String, Any?>()
        for (parameter in parameters) {
            try {
                bound[parameter.name] = parameter.resolve(values[parameter.name])
            } catch (e: InvalidValue) {
                problems += ParameterProblem(parameter.name, e.message ?: "invalid value")
            }
        }
        if (problems.isNotEmpty()) throw InvalidParametersException(id, problems)
        return ParameterValues(id, bound)
    }

    /** Runs the operation on [input] with already-bound [values]. */
    suspend fun run(input: I, values: ParameterValues, progress: ProgressReporter = ProgressReporter.None): O {
        require(values.operationId == id) { "values were bound for '${values.operationId}', not '$id'" }
        return RunScope(values).body(input, progress)
    }

    /** Binds [parameters] (see [bind]) and runs the operation on [input]. */
    suspend fun run(
        input: I,
        parameters: Map<String, Any?> = emptyMap(),
        progress: ProgressReporter = ProgressReporter.None,
    ): O = run(input, bind(parameters), progress)

    @Suppress("UNCHECKED_CAST")
    internal suspend fun runUnchecked(input: Any, values: ParameterValues, progress: ProgressReporter): O {
        require(inputType.isInstance(input)) {
            "'$id' expects ${inputType.simpleName} input, got ${input::class.simpleName}"
        }
        return run(input as I, values, progress)
    }

    override fun toString(): String = "Operation($id)"
}
