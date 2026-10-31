package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.main
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.tassel.tasselVersion

class Tassel : CliktCommand(name = "tassel") {
    init {
        versionOption(tasselVersion())
    }

    override fun run() = Unit
}

fun main(args: Array<String>) = Tassel().main(args)
