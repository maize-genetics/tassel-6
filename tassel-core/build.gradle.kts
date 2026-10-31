plugins {
    alias(libs.plugins.kotlin.multiplatform)
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
        commonMain {
            kotlin.srcDir(generateBuildInfo.map { it.outputDir })
        }
        commonTest.dependencies {
            implementation(libs.kotlin.test)
        }
    }
}
