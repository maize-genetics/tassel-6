package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.tassel.ops.BooleanParameterDescriptor
import net.maizegenetics.tassel.ops.DoubleParameterDescriptor
import net.maizegenetics.tassel.ops.EnumParameterDescriptor
import net.maizegenetics.tassel.ops.FileParameterDescriptor
import net.maizegenetics.tassel.ops.IntParameterDescriptor
import net.maizegenetics.tassel.ops.ParameterDescriptor
import net.maizegenetics.tassel.ops.TasselOperations
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNotNull
import kotlin.test.assertTrue

/** A long flag such as `--max-alleles`, or a one-letter short flag. */
private val KEBAB_FLAG = Regex("--[a-z0-9]+(-[a-z0-9]+)*|-[a-z]")

/** Help text with every run of whitespace collapsed, so assertions do not depend on line wrapping. */
private fun help(command: String): String {
    val result = tasselCommand().test("$command --help")
    assertEquals(0, result.statusCode, "`tassel $command --help` failed:\n${result.output}")
    return result.output.replace(Regex("\\s+"), " ")
}

private fun flat(text: String) = text.replace(Regex("\\s+"), " ")

class OperationCommandTest {
    private val registry = TasselOperations.registry

    @Test
    fun everyRegisteredOperationBecomesASubcommand() {
        val generated = operationCommands(registry).map { it.commandName }
        assertEquals(registry.operations.map { it.id }, generated)
        val root = tasselCommand().registeredSubcommandNames()
        assertEquals(registry.operations.map { it.id } + "dev", root)
    }

    @Test
    fun rootHelpListsEveryOperationWithItsSpecHelp() {
        val text = help("")
        for (operation in registry.operations) {
            assertContains(text, operation.id)
            if (operation.help.isNotEmpty()) assertContains(text, flat(operation.help))
        }
    }

    @Test
    fun helpShowsEachParametersSpecHelpAndDefault() {
        for (operation in registry.operations) {
            val text = help(operation.id)
            assertContains(text, flat(operation.help), message = "${operation.id} help text")
            for (parameter in operation.descriptor.parameters) {
                val flag = "--" + kebabCase(parameter.name)
                assertContains(text, flag, message = "${operation.id} is missing $flag")
                if (parameter.help.isNotEmpty()) {
                    assertContains(text, flat(parameter.help).take(40), message = "${operation.id} $flag help")
                }
                defaultText(parameter)?.let {
                    assertContains(text, "default: $it", message = "${operation.id} $flag default")
                }
            }
        }
    }

    @Test
    fun legacyTassel5NamesAreNotOptions() {
        for (operation in registry.operations) {
            val options = tasselCommand().registeredSubcommands()
                .single { it.commandName == operation.id }
                .registeredOptions()
                .flatMap { it.names }
            for (parameter in operation.parameters) {
                val legacy = parameter.legacy ?: continue
                // A legacy name that is already the kebab-case parameter name (summary's
                // `overview`) is not a second spelling; TASSEL 5's own `-overview` is the same word.
                if (legacy == kebabCase(parameter.name)) continue
                assertFalse(
                    "--$legacy" in options || "-$legacy" in options,
                    "${operation.id} exposes the TASSEL 5 name '$legacy' as an option",
                )
            }
            for (option in options) {
                assertTrue(KEBAB_FLAG.matches(option), "${operation.id} option '$option' is not kebab-case")
            }
        }
    }

    @Test
    fun noSpecClaimsTheNamesTheCliReserves() {
        for (operation in registry.operations) {
            for (parameter in operation.parameters) {
                assertTrue(
                    parameter.name !in setOf("input", "output"),
                    "'${operation.id}' declares '${parameter.name}', which the CLI uses for its own option",
                )
            }
        }
    }

    @Test
    fun theFlagsBlockedByTheMissingStringParameterTypeAreNotHandWritten() {
        // tassel-6#45 adds string and string-list parameters to the DSL; until then filter-sites
        // cannot declare startChr/endChr and filter-taxa cannot declare taxaList, so the generator
        // must not show those flags. When #45 declares them, these assertions flip by themselves.
        val siteFlags = help("filter-sites")
        assertFalse("--start-chr" in siteFlags || "--end-chr" in siteFlags)
        assertFalse("--taxa-list" in help("filter-taxa"))
    }

    @Test
    fun aValueOutsideTheSpecsRangeIsAUserError() {
        val result = tasselCommand().test(
            "kinship --input ${testDataPath("tutorial/mdp_genotype.hmp.txt")} --max-alleles 9",
        )
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "--max-alleles")
        assertContains(result.output, "outside range [2..6]")
    }

    @Test
    fun anUnknownOptionIsAUserError() {
        val result = tasselCommand().test("kinship --nonesuch 1")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "no such option")
    }

    @Test
    fun aMissingInputFileIsAUserError() {
        val result = tasselCommand().test("summary --input no/such/file.hmp.txt")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "no such file")
    }

    @Test
    fun helpSpellsEveryEnumValueExactlyAsItHasToBeTyped() {
        // The help formatter lowercases a metavar with no alternatives, which would corrupt a
        // constant name such as Centered_IBS. Kebab-case values survive it. A Clikt upgrade that
        // changed how a metavar is rendered would break this test rather than the user's copy-paste.
        for (operation in registry.operations) {
            val text = help(operation.id)
            val enums = operation.descriptor.parameters.filterIsInstance<EnumParameterDescriptor>()
            for (parameter in enums) {
                val flag = "--" + kebabCase(parameter.name)
                val shown = Regex("\\Q$flag\\E=(\\S+)").find(text)?.groupValues?.get(1)
                assertNotNull(shown, "${operation.id} $flag has no metavar in:\n$text")
                val offered = shown.trim('<', '>', '(', ')').split("|")
                assertEquals(parameter.choices.map { kebabCase(it) }, offered, "${operation.id} $flag choices")
                parameter.default?.let { assertContains(text, "default: ${kebabCase(it)}") }
            }
        }
    }

    @Test
    fun anEnumValueCopiedOutOfTheHelpTextIsAccepted() {
        val directory = tempDirectory("tassel-enum")
        val input = directory.resolve("in.hmp.txt")
        writeHapMap(simulateGenotypes(taxa = 4, sites = 8, seed = 3), input)
        val shown = Regex("--method=(\\S+)").find(help("kinship"))!!.groupValues[1]
        val offered = shown.trim('<', '>', '(', ')').split("|")
        assertEquals(listOf("centered-ibs"), offered)

        for (value in offered) {
            val result = tasselCommand().test(
                "--quiet kinship --input $input --method $value --output ${directory.resolve(value)}",
            )
            assertEquals(0, result.statusCode, result.output)
        }
    }

    @Test
    fun theTassel5SpellingOfAnEnumValueIsRejectedWithTheKebabCaseOne() {
        // -method Centered_IBS is TASSEL 5's; it belongs to `tassel translate` (ADR-0006), so here
        // it has to fail during parsing with an error naming the spelling this CLI takes.
        val input = testDataPath("tutorial/mdp_genotype.hmp.txt")
        val result = tasselCommand().test("kinship --input $input --method Centered_IBS")

        assertEquals(EXIT_USER_ERROR, result.statusCode)
        val message = result.output.replace(Regex("\\s+"), " ")
        assertContains(message, "invalid choice: Centered_IBS")
        assertContains(message, "centered-ibs")
    }

    /** How `--help` should spell a parameter's default: the way the user would type it. */
    private fun defaultText(parameter: ParameterDescriptor): String? = when (parameter) {
        is IntParameterDescriptor -> parameter.default?.toString()
        is DoubleParameterDescriptor -> parameter.default?.toString()
        is BooleanParameterDescriptor -> parameter.default?.toString()
        is EnumParameterDescriptor -> parameter.default?.let { kebabCase(it) }
        is FileParameterDescriptor -> parameter.default
    }
}
