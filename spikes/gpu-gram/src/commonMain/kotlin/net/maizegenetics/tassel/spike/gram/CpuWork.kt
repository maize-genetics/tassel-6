package net.maizegenetics.tassel.spike.gram

/** Times the raw-dosage f64 Gram. One discarded warm-up, then [runs] samples. */
fun timeDosageShapes(
    runs: Int,
    shapes: List<GramShape> = GRAM_SHAPES,
    log: (String) -> Unit = ::println,
): List<CpuSizeJson> = shapes.map { shape ->
    log("dosage ${shape.label}")
    val codes = generateDosages(shape.taxa, shape.sites)
    val checksum = dosageChecksum(codes).toString()
    val preview = codes.take(16).map { it.toInt() }
    val dosage = timeRuns(runs) { dosageGram(codes, shape.taxa, shape.sites) }
    log("dosage ${shape.label} median ${dosage.medianSeconds}")
    CpuSizeJson(
        taxa = shape.taxa,
        sites = shape.sites,
        checksum = checksum,
        preview = preview,
        dosage = dosage.toJson(),
    )
}
