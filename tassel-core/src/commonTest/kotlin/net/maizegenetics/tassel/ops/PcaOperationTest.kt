package net.maizegenetics.tassel.ops

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.analysis.DEFAULT_COMPONENTS
import net.maizegenetics.tassel.analysis.PrincipalComponentsResult
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotNull
import kotlin.test.assertNull

/**
 * The `pca` spec: TASSEL 5's parameter names and defaults, and that running it is the same as
 * calling the kernel.
 *
 * @see net.maizegenetics.analysis.data.PrincipalComponentsPlugin
 */
class PcaOperationTest {
    @Test
    fun isRegisteredWithTassel5LegacyNames() {
        assertEquals(pcaOperation, TasselOperations.registry.get("pca"))
        assertEquals(pcaOperation, TasselOperations.registry.byLegacyPlugin("PrincipalComponentsPlugin"))
        assertEquals(GenotypeTable::class, pcaOperation.inputType)
        assertEquals(PrincipalComponentsResult::class, pcaOperation.outputType)
    }

    @Test
    fun carriesEveryTassel5ParameterWithItsDefault() {
        val expected = mapOf(
            "covariance" to ("covariance" to true),
            "limitBy" to ("limitBy" to PcaLimit.number_of_components),
            "numberOfComponents" to ("ncomponents" to DEFAULT_COMPONENTS),
            "minimumEigenvalue" to ("minEigenval" to 0.0),
            "totalVariance" to ("totalVar" to 0.5),
            "reportEigenvalues" to ("reportEigenvalues" to true),
            "reportEigenvectors" to ("reportEigenvectors" to true),
        )
        assertEquals(expected.keys.toList(), pcaOperation.parameters.map { it.name })
        expected.forEach { (name, legacyAndDefault) ->
            val (legacy, default) = legacyAndDefault
            assertEquals(default, pcaOperation.parameter(name)?.default, name)
            assertEquals(name, pcaOperation.parameterByLegacyName(legacy)?.name, legacy)
        }
    }

    @Test
    fun runsWithDefaultsAndWithBothTablesTurnedOff() = runTest {
        val defaults = pcaOperation.run(fourSites)
        // Four taxa, so the five default components are capped at four.
        val expected = minOf(DEFAULT_COMPONENTS, fourSites.numberOfTaxa) + 1
        assertEquals(expected, defaults.principalComponents.numberOfAttributes)
        assertNotNull(defaults.eigenvalues)
        assertNotNull(defaults.eigenvectors)

        val quiet = pcaOperation.run(
            fourSites,
            mapOf("reportEigenvalues" to false, "reportEigenvectors" to false, "numberOfComponents" to 2),
        )
        assertNull(quiet.eigenvalues)
        assertNull(quiet.eigenvectors)
        assertEquals(3, quiet.principalComponents.numberOfAttributes)
    }

    @Test
    fun acceptsTassel5CommandLineTextForEveryParameter() = runTest {
        val values = pcaOperation.bind(
            mapOf(
                "covariance" to "false",
                "limitBy" to "number_of_components",
                "numberOfComponents" to "3",
                "minimumEigenvalue" to "0.25",
                "totalVariance" to "0.9",
            ),
        )
        val result = pcaOperation.run(fourSites, values)
        assertEquals(4, result.principalComponents.numberOfAttributes)
    }

    /**
     * TASSEL 5 chooses the count for these two with `Arrays.binarySearch` over an array its SVD
     * leaves unsorted, so there is no defined behavior to port; #19 scoped them out. The parameters
     * stay declared so the legacy translator and rTASSEL keep resolving their names.
     */
    @Test
    fun refusesTheTwoLimitModesThatWereNotPorted() = runTest {
        listOf(PcaLimit.min_eigenvalue, PcaLimit.total_variance).forEach { limit ->
            val failure = assertFailsWith<IllegalArgumentException> {
                pcaOperation.run(fourSites, mapOf("limitBy" to limit.name))
            }
            assertEquals(true, failure.message?.contains(limit.name), failure.message)
        }
    }

    @Test
    fun rejectsATotalVarianceOutsideTassel5sRange() {
        assertFailsWith<InvalidParametersException> { pcaOperation.bind(mapOf("totalVariance" to 1.5)) }
    }

    private companion object {
        /** Four taxa over four sites; the spec tests only need the kernel to run. */
        val fourSites: GenotypeTable by lazy {
            val names = listOf("t0", "t1", "t2", "t3")
            val genotypes = listOf("AACC", "ACCC", "CCAC", "ACAA")
            val header = (HapMap.FIXED_COLUMNS + names).joinToString("\t")
            val lines = (0 until 4).map { site ->
                val fixed = listOf("s$site", "A/C", "1", "${site + 1}", "+", "NA", "NA", "NA", "NA", "NA", "NA")
                (fixed + genotypes.map { "${it[site]}${it[site]}" }).joinToString("\t")
            }
            val text = (listOf(header) + lines).joinToString("\n", postfix = "\n")
            HapMap.read(Buffer().also { it.write(text.encodeToByteArray()) })
        }
    }
}
