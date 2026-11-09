import org.jetbrains.kotlin.gradle.targets.js.testing.KotlinJsTest
import org.jetbrains.kotlin.gradle.targets.native.tasks.KotlinNativeTest

plugins {
    alias(libs.plugins.kotlin.multiplatform)
    alias(libs.plugins.kotlin.serialization)
    alias(libs.plugins.ktlint)
    alias(libs.plugins.kover)
}

ktlint {
    filter {
        exclude { it.file.path.contains("/build/generated/") }
    }
}

abstract class GenerateBuildInfo : DefaultTask() {
    @get:Input
    abstract val version: Property<String>

    @get:OutputDirectory
    abstract val outputDir: DirectoryProperty

    @TaskAction
    fun generate() {
        val file = outputDir.file("net/maizegenetics/tassel/BuildInfo.kt").get().asFile
        file.parentFile.mkdirs()
        file.writeText(
            """
            |package net.maizegenetics.tassel
            |
            |internal object BuildInfo {
            |    const val VERSION: String = "${version.get()}"
            |}
            |
            """.trimMargin(),
        )
    }
}

val generateBuildInfo = tasks.register<GenerateBuildInfo>("generateBuildInfo") {
    version = project.version.toString()
    outputDir = layout.buildDirectory.dir("generated/buildInfo/commonMain/kotlin")
}

// commonTest reads testdata/ through TestData, which finds it through this variable on every target.
val testDataDir = rootProject.layout.projectDirectory.dir("testdata")
val testDataEnv = "TASSEL_TESTDATA_DIR"

tasks.withType<AbstractTestTask>().configureEach {
    inputs.dir(testDataDir).withPathSensitivity(PathSensitivity.RELATIVE).withPropertyName("testData")
}
tasks.withType<Test>().configureEach { environment(testDataEnv, testDataDir.asFile.absolutePath) }
tasks.withType<KotlinJsTest>().configureEach { environment(testDataEnv, testDataDir.asFile.absolutePath) }
tasks.withType<KotlinNativeTest>().configureEach { environment(testDataEnv, testDataDir.asFile.absolutePath) }

kotlin {
    jvmToolchain(libs.versions.jdk.get().toInt())

    jvm()
    wasmJs {
        browser {
            // Core tests run on Node; the browser runtime is exercised by tassel-app.
            testTask { enabled = false }
        }
        nodejs()
    }
    macosArm64()
    linuxX64()

    sourceSets {
        // Node fs must not live in wasmJsMain: that source set is linked into the browser (ADR-0021).
        // The default hierarchy forbids a source set depending on wasmJsMain, so this directory is compiled
        // only into wasmJsTest. The production wasm klib, and the browser bundle, do not include it.
        named("wasmJsTest") {
            kotlin.srcDir("src/wasmJsNodeMain/kotlin")
        }
        commonMain {
            kotlin.srcDir(generateBuildInfo.map { it.outputDir })
            dependencies {
                implementation(libs.kotlinx.coroutines.core)
                // Readers and writers take kotlinx.io Source and Sink, so this is part of the public API (ADR-0017).
                api(libs.kotlinx.io.core)
                api(libs.kotlinx.serialization.json)
            }
        }
        commonTest.dependencies {
            implementation(libs.kotlin.test)
            implementation(libs.kotlinx.coroutines.test)
        }
    }
}
