package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.Context
import com.github.ajalt.clikt.parameters.options.OptionWithValues
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.choice
import kotlinx.coroutines.runBlocking
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.ops.BooleanParameterDescriptor
import net.maizegenetics.tassel.ops.DoubleParameterDescriptor
import net.maizegenetics.tassel.ops.EnumParameterDescriptor
import net.maizegenetics.tassel.ops.FileKind
import net.maizegenetics.tassel.ops.FileParameterDescriptor
import net.maizegenetics.tassel.ops.IntParameterDescriptor
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.ops.ParameterDescriptor
import net.maizegenetics.tassel.ops.TasselOperations

/** `--input` and `--output` are the CLI's own; no spec may declare a parameter with either name. */
private val RESERVED_NAMES = setOf("input", "output")

/**
 * A subcommand for every registered operation, in the registry's presentation order (ADR-0002 and
 * ADR-0006). Nothing here knows the name of an operation or of a parameter: a new spec gets its
 * subcommand, its options, and its `--help` with no change to this module.
 */
fun operationCommands(registry: OperationRegistry = TasselOperations.registry): List<CliktCommand> =
    registry.operations.map { OperationCommand(it.descriptor, registry) }

/**
 * One generated subcommand. Options are the spec's parameter names in kebab-case; the spec's
 * `legacy` TASSEL 5 names are not registered, because translating those is `tassel translate`'s job
 * (ADR-0006, tassel-6#26).
 */
class OperationCommand(
    private val descriptor: OperationDescriptor,
    private val registry: OperationRegistry = TasselOperations.registry,
) : CliktCommand(name = descriptor.id) {
    private val input = if (descriptor.inputType != "Unit") {
        option("--input", "-i", help = "Genotype file to analyze (HapMap)", metavar = "PATH")
            .required()
            .also { registerOption(it) }
    } else {
        null
    }

    private val output = option(
        "--output",
        "-o",
        help = "Base name for the result files; omit to write to standard output",
        metavar = "PATH",
    ).also { registerOption(it) }

    /** Reads each option as the command-line text the engine takes, or `null` when not supplied. */
    private val parameters: List<Pair<String, () -> String?>> = descriptor.parameters.map { parameter ->
        require(parameter.name !in RESERVED_NAMES) {
            "'${descriptor.id}' declares a parameter called '${parameter.name}', which the CLI reserves"
        }
        parameter.name to registerParameter(parameter)
    }

    override fun help(context: Context): String = descriptor.help

    override fun run() = cliRun {
        if (input != null && descriptor.inputType != "GenotypeTable") {
            // The only reader the CLI has is HapMap, so an operation over anything else needs a
            // reader here before its subcommand can run.
            unsupported("'${descriptor.id}' takes a ${descriptor.inputType}, which the CLI cannot read yet")
        }
        val supplied = parameters.mapNotNull { (name, read) -> read()?.let { name to it } }.toMap()
        // Bind before reading, so bad options are reported without first parsing a large file. The
        // engine binds again from the same text; that is the contract RunRequest carries.
        registry[descriptor.id]?.bind(supplied)
        val genotypes = input?.let { readGenotypes(it.value, stderrProgress("Reading")) }
        val engine = LocalEngine(registry)
        runBlocking {
            val handle = genotypes?.let { engine.put(it) }
            val request = RunRequest(descriptor.id, handle, supplied)
            CliLog.info("running ${descriptor.id}")
            val result = engine.get(engine.run(request, stderrProgress(descriptor.id)))
            writeResult(result, output.value, stderrProgress("Writing"))
        }
    }

    private fun registerParameter(parameter: ParameterDescriptor): () -> String? {
        val name = "--" + kebabCase(parameter.name)
        return when (parameter) {
            is BooleanParameterDescriptor -> {
                val fallback = parameter.default
                if (fallback == null) {
                    // No spec declares a boolean without a default yet; a plain value keeps "not
                    // supplied" distinguishable from "supplied false".
                    text(option(name, help = parameter.help, metavar = "true|false"))
                } else {
                    val flag = option(name, help = parameter.help).flag(
                        "--no-" + kebabCase(parameter.name),
                        default = fallback,
                        defaultForHelp = fallback.toString(),
                    )
                    registerOption(flag)
                    val read: () -> String = { flag.value.toString() }
                    read
                }
            }
            is EnumParameterDescriptor -> {
                val choices = enumChoices(parameter.choices)
                text(
                    option(name, help = parameter.help).choice(choices, metavar = choiceMetavar(choices.keys)),
                    default = parameter.default,
                    isRequired = parameter.required,
                    defaultForHelp = parameter.default?.let { kebabCase(it) },
                )
            }
            is IntParameterDescriptor -> text(
                option(name, help = withRange(parameter.help, parameter.min, parameter.max), metavar = "INT"),
                parameter.default?.toString(),
                parameter.required,
            )
            is DoubleParameterDescriptor -> text(
                option(name, help = withRange(parameter.help, parameter.min, parameter.max), metavar = "NUMBER"),
                // The spec's own text for the default, so --help and the engine agree.
                parameter.default?.toString(),
                parameter.required,
            )
            is FileParameterDescriptor -> text(
                option(name, help = parameter.help, metavar = fileMetavar(parameter.kind)),
                parameter.default,
                parameter.required,
            )
        }
    }

    /**
     * Registers a string-valued option and returns its value, applying the spec's default if any.
     * [defaultForHelp] is how `--help` spells that default when it is not how the user types it.
     */
    private fun text(
        option: OptionWithValues<String?, String, String>,
        default: String? = null,
        isRequired: Boolean = false,
        defaultForHelp: String? = null,
    ): () -> String? {
        if (default != null) {
            val withDefault = option.default(default, defaultForHelp ?: default)
            registerOption(withDefault)
            return { withDefault.value }
        }
        if (isRequired) {
            val required = option.required()
            registerOption(required)
            return { required.value }
        }
        registerOption(option)
        return { option.value }
    }
}

/**
 * The only spelling the CLI accepts for an enum value: the spec's constant name in kebab-case
 * (`--method centered-ibs`), which is the house style in `.cursor/rules/cli.mdc`. TASSEL 5's own
 * spelling (`-method Centered_IBS`) belongs to `tassel translate` under ADR-0006, not here. The
 * engine is sent the constant name, which is the only form the spec parses.
 *
 * @return accepted value to constant name.
 */
private fun enumChoices(choices: List<String>): Map<String, String> {
    val accepted = choices.associateBy { kebabCase(it) }
    require(accepted.size == choices.size) { "enum constants share a kebab-case spelling: $choices" }
    return accepted
}

/**
 * What `--help` shows for an enum. The accepted values are kebab-case and so already lowercase,
 * which is what makes this safe: the help formatter lowercases a metavar that offers no
 * alternatives, and lowercasing a kebab-case value cannot change it. So the help always spells a
 * value exactly as it has to be typed, whether there is one choice or several.
 */
private fun choiceMetavar(accepted: Collection<String>): String =
    if (accepted.size == 1) accepted.first() else accepted.joinToString("|", "(", ")")

/** Spells out the spec's inclusive range, which `--help` would otherwise not show. */
private fun withRange(help: String, min: Any?, max: Any?): String = when {
    min != null && max != null -> "$help (from $min to $max)"
    min != null -> "$help (at least $min)"
    max != null -> "$help (at most $max)"
    else -> help
}

private fun fileMetavar(kind: FileKind): String = when (kind) {
    FileKind.Input, FileKind.Output -> "PATH"
    FileKind.InputDirectory, FileKind.OutputDirectory -> "DIR"
}

/** `maxAlleles` becomes `max-alleles` and `Centered_IBS` becomes `centered-ibs`. */
internal fun kebabCase(name: String): String = buildString {
    for ((index, character) in name.withIndex()) {
        when {
            character == '_' || character == ' ' -> append('-')
            character.isUpperCase() -> {
                if (index > 0 && name[index - 1] != '_' && !name[index - 1].isUpperCase()) append('-')
                append(character.lowercaseChar())
            }
            else -> append(character)
        }
    }
}
