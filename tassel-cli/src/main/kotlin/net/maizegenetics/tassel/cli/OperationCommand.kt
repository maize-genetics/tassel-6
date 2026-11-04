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
                    option(name, help = parameter.help)
                        .choice(choices, metavar = choiceMetavar(parameter.choices, choices.keys), ignoreCase = true),
                    parameter.default,
                    parameter.required,
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

    /** Registers a string-valued option and returns its value, applying the spec's default if any. */
    private fun text(
        option: OptionWithValues<String?, String, String>,
        default: String? = null,
        isRequired: Boolean = false,
    ): () -> String? {
        if (default != null) {
            val withDefault = option.default(default)
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
 * Enum values the CLI accepts: the spec's constant names, which TASSEL 5 also accepts
 * (`-method Centered_IBS`), plus a kebab-case spelling where that differs
 * (`--method centered-ibs`, as `.cursor/rules/cli.mdc` writes it). Matching ignores letter case.
 * Either way the engine is sent the constant name, which is the only form the spec parses.
 */
private fun enumChoices(choices: List<String>): Map<String, String> = buildMap {
    for (choice in choices) {
        put(choice, choice)
        val kebab = kebabCase(choice)
        if (!kebab.equals(choice, ignoreCase = true)) put(kebab, choice)
    }
}

/**
 * What `--help` shows for an enum: the spec's constant names, which are also what TASSEL 5 accepts.
 * A one-constant enum lists its kebab-case alias too, because the help formatter lowercases a
 * metavar that offers no alternatives, and `(centered_ibs)` would be the wrong spelling to show.
 */
private fun choiceMetavar(choices: List<String>, accepted: Set<String>): String =
    (if (choices.size > 1) choices else accepted.toList()).joinToString("|", "(", ")")

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
