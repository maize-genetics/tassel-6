package net.maizegenetics.tassel

import kotlin.test.Test
import kotlin.test.assertTrue

class VersionTest {
    @Test
    fun versionIsTassel6SemVer() {
        val version = tasselVersion()
        val semVer = Regex("""6\.\d+\.\d+(-[0-9A-Za-z.-]+)?""")
        assertTrue(semVer.matches(version), "not a TASSEL 6 SemVer string: '$version'")
    }
}
