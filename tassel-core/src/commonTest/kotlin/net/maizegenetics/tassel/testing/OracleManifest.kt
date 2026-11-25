package net.maizegenetics.tassel.testing

/**
 * `testdata/oracle/MANIFEST`, written by `scripts/oracle/generate-oracles.sh` (ADR-0029, kind 2).
 *
 * The syntax is [GoldenManifest]'s: `[section]` headers and `key = value` lines, where a key that
 * starts with `oracle/` is a fixture file (relative to `testdata/`) and its value is its
 * [Tolerance]. Tool versions and image digests are plain properties of the `[tools]` section.
 */
object OracleManifest {
    const val PATH = "oracle/MANIFEST"
    const val OUTPUT_PREFIX = "oracle/"

    fun load(): GoldenManifest = parse(TestData.readText(PATH))

    fun parse(text: String): GoldenManifest = GoldenManifest.parse(text, OUTPUT_PREFIX)
}
