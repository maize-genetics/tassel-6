@file:OptIn(ExperimentalWasmDsl::class)

import org.gradle.api.tasks.JavaExec
import org.gradle.jvm.toolchain.JavaLanguageVersion
import org.jetbrains.kotlin.gradle.ExperimentalWasmDsl
import org.jetbrains.kotlin.gradle.targets.js.nodejs.NodeJsExec
import org.jetbrains.kotlin.gradle.targets.jvm.KotlinJvmTarget

plugins {
    alias(libs.plugins.kotlin.multiplatform)
    alias(libs.plugins.kotlin.serialization)
    alias(libs.plugins.ktlint)
}

kotlin {
    jvmToolchain(libs.versions.jdk.get().toInt())

    jvm()
    wasmJs {
        browser {
            commonWebpackConfig {
                outputFileName = "gpuGram.js"
            }
            // GPU runs are a manual browser session. check must not launch one.
            testTask { enabled = false }
        }
        nodejs()
        binaries.executable()
    }

    sourceSets {
        commonMain.dependencies {
            implementation(project(":tassel-core"))
            implementation(libs.kotlinx.coroutines.core)
            implementation(libs.kotlinx.serialization.json)
        }
        commonTest.dependencies {
            implementation(libs.kotlin.test)
            implementation(libs.kotlinx.coroutines.test)
        }
    }
}

val jvmCompilation = kotlin.targets.getByName<KotlinJvmTarget>("jvm").compilations.getByName("main")
val launcher = javaToolchains.launcherFor {
    languageVersion.set(JavaLanguageVersion.of(libs.versions.jdk.get().toInt()))
}
val reportDirectory = layout.buildDirectory.dir("reports/gpu-gram")
val distributionDirectory = layout.buildDirectory.dir("dist/wasmJs/productionExecutable")
val shaderFile = layout.projectDirectory.file("wgsl/gram_f32.wgsl")

fun JavaExec.spikeMain(main: String) {
    group = "application"
    classpath(jvmCompilation.output.allOutputs, jvmCompilation.runtimeDependencyFiles)
    mainClass.set(main)
    javaLauncher.set(launcher)
    dependsOn(jvmCompilation.compileTaskProvider)
}

tasks.register<JavaExec>("cpuBench") {
    spikeMain("net.maizegenetics.tassel.spike.gram.CpuBenchKt")
    description = "Time the f64 Gram on the JVM and write build/reports/gpu-gram/jvm-cpu.json."
    val output = reportDirectory.map { it.file("jvm-cpu.json").asFile.absolutePath }
    argumentProviders.add { listOf(output.get()) }
    outputs.file(reportDirectory.map { it.file("jvm-cpu.json") })
    outputs.upToDateWhen { false }
}

tasks.named<NodeJsExec>("wasmJsNodeProductionRun") {
    val report = layout.buildDirectory.file("reports/gpu-gram/node-cpu.json")
    // The Node process starts in the npm package directory, so a relative report path would miss
    // the module's build/reports directory. Set the variable the wasm main already reads.
    environment("TASSEL_GPU_GRAM_OUT", report.get().asFile.absolutePath)
    outputs.file(report)
    outputs.upToDateWhen { false }
}

tasks.register<JavaExec>("serveGpuGram") {
    spikeMain("net.maizegenetics.tassel.spike.gram.ServeGpuGramKt")
    description = "Serve the browser spike and collect JSON reports under build/reports/gpu-gram/."
    dependsOn(tasks.named("wasmJsBrowserDistribution"))
    val distribution = distributionDirectory.map { it.asFile.absolutePath }
    val reports = reportDirectory.map { it.asFile.absolutePath }
    argumentProviders.add {
        listOf(distribution.get(), shaderFile.asFile.absolutePath, reports.get(), "8765")
    }
    outputs.upToDateWhen { false }
}
