package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.Context
import com.github.ajalt.clikt.parameters.options.OptionWithValues
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.multiple
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
import net.maizegenetics.tassel.ops.StringListParameterDescriptor
import net.maizegenetics.tassel.ops.StringParameterDescriptor
import net.maizegenetics.tassel.ops.TasselOperations

/**
 * `--input`, `--output`, `--region`, and `--index` are the CLI's own.
 * No spec may declare a parameter with one of these names.
 */
internal val RESERVED_NAMES = setOf("input", "output", "region", "index")

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
 * (ADR-0006, tassel-6#26). `--region` and `--index` are registered here, like `--input` and
 * `--output`, so every genotype command has them without a spec change (ADR-0002).
 */
class OperationCommand(
    private val descriptor: OperationDescriptor,
    private val registry: OperationRegistry = TasselOperations.registry,
) : CliktCommand(name = descriptor.id) {
    /** The `--input` paths, in command-line order; empty for an operation that takes no input. */
    private val inputs: () -> List<String> = when (descriptor.inputType) {
        "Unit" -> {
            { emptyList() }
        }
        "DataSet" -> {
            val option = option(
                "--input",
                "-i",
                help = "File to include. Repeat it for each input, in order. A phenotype file " +
                    "(<Phenotype>, <Trait>, <Covariate>, <Factor>, or <Header>) is read as a phenotype; " +
                    "anything else is a genotype (HapMap, VCF, or BCF). Detected from the content, not the name.",
                metavar = "PATH",
            ).multiple(required = true).also { registerOption(it) }
            ({ option.value })
        }
        "PhenotypeTable" -> {
            val option = option(
                "--input",
                "-i",
                help = "Phenotype file to analyze (<Phenotype>, or TASSEL 5's <Trait>, <Covariate>, <Factor>, " +
                    "and <Header> layouts)",
                metavar = "PATH",
            ).required().also { registerOption(it) }
            ({ listOf(option.value) })
        }
        else -> {
            val option = option(
                "--input",
                "-i",
                help = "Genotype file to analyze (HapMap, VCF, or BCF)",
                metavar = "PATH",
            )
                .required()
                .also { registerOption(it) }
            ({ listOf(option.value) })
        }
    }

    private val takesGenotypes = descriptor.inputType == "GenotypeTable" || descriptor.inputType == "DataSet"

    private val region = if (takesGenotypes) {
        option(
            "--region",
            help = "Read one region of a bgzipped, indexed VCF or BCF (.bcf), for example " +
                "chr3B:600000000-650000000. A bare contig name reads that contig. Coordinates are " +
                "1-based and inclusive. A BCF region needs a .csi.",
            metavar = "TEXT",
        ).also { registerOption(it) }
    } else {
        null
    }

    private val index = if (takesGenotypes) {
        option(
            "--index",
            help = "CSI or TBI index for --region. When omitted, look beside the VCF for .csi, then .tbi. " +
                "A BCF uses .csi only.",
            metavar = "PATH",
        ).also { registerOption(it) }
    } else {
        null
    }

    private val output = option(
        "--output",
        "-o",
        help = "Where to write the result. For a genotype, the suffix selects the format, ignoring case: " +
            "filtered.hmp.txt, filtered.vcf, filtered.vcf.gz, or filtered.bcf. " +
            "Any other name stays HapMap, and so does standard output. A phenotype is written as <Phenotype> " +
            "text, and a result of several tables as one file per table. " +
            "Omit to write to standard output.",
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
        if (descriptor.inputType !in setOf("Unit", "GenotypeTable", "PhenotypeTable", "DataSet")) {
            // An operation over anything else needs a reader here before its subcommand can run.
            unsupported("'${descriptor.id}' takes a ${descriptor.inputType}, which the CLI cannot read yet")
        }
        val supplied = parameters.mapNotNull { (name, read) -> read()?.let { name to it } }.toMap()
        // Bind before reading, so bad options are reported without first parsing a large file. The
        // engine binds again from the same text; that is the contract RunRequest carries.
        registry[descriptor.id]?.bind(supplied)
        val regionText = region?.value
        val indexText = index?.value
        if (indexText != null && regionText == null) {
            throw IllegalArgumentException("--index requires --region")
        }
        val engine = LocalEngine(registry)
        runBlocking {
            val data = readInput(regionText, indexText)
            val handle = data?.let { engine.put(it) }
            val request = RunRequest(descriptor.id, handle, supplied)
            CliLog.info("running ${descriptor.id}")
            val result = engine.get(engine.run(request, stderrProgress(descriptor.id)))
            writeResult(result, output.value, stderrProgress("Writing"))
        }
    }

    /** Reads `--input` as the type the operation takes, or `null` when it takes none. */
    private suspend fun readInput(region: String?, index: String?): Any? {
        val paths = inputs()
        val progress = stderrProgress("Reading")
        return when (descriptor.inputType) {
            "GenotypeTable" -> {
                val path = paths.single()
                if (sniffInput(path) == InputKind.Phenotype) {
                    throw IllegalArgumentException(
                        "$path is a phenotype file, but '${descriptor.id}' takes a genotype (HapMap, VCF, or BCF)",
                    )
                }
                readGenotypes(path, progress, region, index)
            }
            "PhenotypeTable" -> {
                val path = paths.single()
                // A missing file falls through to readPhenotype, which names it.
                if (isRegularFile(path) && sniffInput(path) != InputKind.Phenotype) {
                    throw IllegalArgumentException(
                        "$path is not a phenotype file, but '${descriptor.id}' takes one; " +
                            "it must start with <Phenotype>, <Trait>, <Covariate>, <Factor>, or <Header>",
                    )
                }
                readPhenotype(path, progress)
            }
            "DataSet" -> readDataSet(paths, progress, region, index)
            else -> null
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
            is StringParameterDescriptor -> text(
                option(name, help = parameter.help, metavar = "TEXT"),
                parameter.default,
                parameter.required,
            )
            is StringListParameterDescriptor -> text(
                // One comma-separated value, so the engine gets the same text a spec file would.
                option(name, help = parameter.help, metavar = "TEXT,TEXT"),
                parameter.default?.joinToString(","),
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
