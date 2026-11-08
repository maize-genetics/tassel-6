package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.Context
import com.github.ajalt.clikt.core.context
import com.github.ajalt.clikt.core.main
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.output.MordantHelpFormatter
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.tassel.tasselVersion

class Tassel : CliktCommand(name = "tassel") {
    private val quiet by option("--quiet", "-q", help = "Only print warnings and errors").flag()

    init {
        context {
            helpFormatter = { MordantHelpFormatter(it, showRequiredTag = true, showDefaultValues = true) }
        }
        versionOption(tasselVersion())
    }

    override fun help(context: Context): String =
        "TASSEL 6: genotype analysis. Every analysis subcommand is generated from its operation spec."

    override fun run() {
        CliLog.quiet = quiet
    }
}

/** The root command with every subcommand: one per registered operation (ADR-0002), then `dev`. */
fun tasselCommand(): CliktCommand =
    Tassel().subcommands(operationCommands() + Dev().subcommands(Simulate(), KinshipBench(), BgzfBench()))

fun main(args: Array<String>) = tasselCommand().main(args)
