package net.maizegenetics.tassel.spike.gram

/**
 * Checksums of [generateDosages] at [PARK_MILLER_SEED], decimal `u64`.
 * Locked by [dosageChecksum] on every target. Issue #113 copies these strings.
 */
val PUBLISHED_CHECKSUMS: Map<GramShape, String> = mapOf(
    GramShape(281, 3_093) to "4544405675908146118",
    GramShape(1_000, 10_000) to "14324033804557704143",
    GramShape(1_000, 100_000) to "3888745768927589105",
    GramShape(5_000, 20_000) to "6644050518704627923",
)

/** First 16 site-major codes of every shape with at least 16 taxa. The first site does not depend on later taxa. */
val TUTORIAL_PREVIEW: List<Int> = listOf(2, 2, 3, 1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 0, 1, 1)
