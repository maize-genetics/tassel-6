plugins {
    alias(libs.plugins.kotlin.jvm)
    alias(libs.plugins.ktlint)
    alias(libs.plugins.kover)
    application
}

kotlin {
    jvmToolchain(libs.versions.jdk.get().toInt())
}

dependencies {
    implementation(project(":tassel-core"))
    implementation(libs.clikt)
    // Engine.run is suspending (ADR-0004), so the CLI needs runBlocking. tassel-core already ships
    // coroutines at run time; this puts them on the compile classpath too.
    implementation(libs.kotlinx.coroutines.core)
    testImplementation(libs.kotlin.test)
}

// The end-to-end tests read testdata/tutorial and testdata/golden, as tassel-core's tests do.
val testDataDir = rootProject.layout.projectDirectory.dir("testdata")
tasks.withType<Test>().configureEach {
    inputs.dir(testDataDir).withPathSensitivity(PathSensitivity.RELATIVE).withPropertyName("testData")
    environment("TASSEL_TESTDATA_DIR", testDataDir.asFile.absolutePath)
}

application {
    applicationName = "tassel"
    mainClass = "net.maizegenetics.tassel.cli.MainKt"
    // Mordant loads JNA for terminal detection; JDK 25 warns unless native access is granted.
    applicationDefaultJvmArgs = listOf("--enable-native-access=ALL-UNNAMED")
}
