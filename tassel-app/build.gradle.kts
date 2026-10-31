plugins {
    alias(libs.plugins.kotlin.multiplatform)
    alias(libs.plugins.kotlin.compose)
    alias(libs.plugins.compose.multiplatform)
}

kotlin {
    jvmToolchain(libs.versions.jdk.get().toInt())

    jvm("desktop")
    wasmJs {
        browser {
            commonWebpackConfig {
                outputFileName = "tasselApp.js"
            }
        }
        binaries.executable()
    }

    sourceSets {
        commonMain.dependencies {
            implementation(project(":tassel-core"))
            implementation(libs.compose.runtime)
            implementation(libs.compose.foundation)
            implementation(libs.compose.ui)
            implementation(libs.compose.material3)
        }
        named("desktopMain").dependencies {
            implementation(compose.desktop.currentOs)
            implementation(libs.kotlinx.coroutines.swing)
        }
    }
}

compose.desktop {
    application {
        mainClass = "net.maizegenetics.tassel.app.MainKt"
        // Compose runs and packages with the Gradle JVM unless told to use the toolchain.
        javaHome = javaToolchains
            .launcherFor { languageVersion = JavaLanguageVersion.of(libs.versions.jdk.get()) }
            .get().metadata.installationPath.asFile.absolutePath
    }
}
