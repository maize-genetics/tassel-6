plugins {
    alias(libs.plugins.kotlin.multiplatform)
    alias(libs.plugins.kotlin.compose)
    alias(libs.plugins.compose.multiplatform)
    alias(libs.plugins.ktlint)
    alias(libs.plugins.kover)
}

// Both targets ship the tutorial dataset: on desktop it lands on the classpath, in the browser it is served next to
// index.html, which is the only way the web build can offer it at all (there is no file system there).
val bundleTutorialData = tasks.register<Sync>("bundleTutorialData") {
    from(rootProject.layout.projectDirectory.file("testdata/tutorial/mdp_genotype.hmp.txt"))
    into(layout.buildDirectory.dir("generated/tutorialData"))
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
            // The view model holds a StateFlow and launches the load (.cursor/rules/compose-ui.mdc).
            implementation(libs.kotlinx.coroutines.core)
        }
        commonTest.dependencies {
            implementation(libs.kotlin.test)
            implementation(libs.kotlinx.coroutines.test)
        }
        named("desktopMain") {
            resources.srcDir(bundleTutorialData)
            dependencies {
                implementation(compose.desktop.currentOs)
                implementation(libs.kotlinx.coroutines.swing)
            }
        }
        named("wasmJsMain") {
            resources.srcDir(bundleTutorialData)
        }
    }
}

compose.desktop {
    application {
        mainClass = "net.maizegenetics.tassel.app.MainKt"
        // Skiko loads its native renderer with System.load; JDK 25 warns unless native access is granted, and a
        // future JDK will block it outright. Same reason as tassel-cli's applicationDefaultJvmArgs for JNA.
        jvmArgs += "--enable-native-access=ALL-UNNAMED"
        // Compose runs and packages with the Gradle JVM unless told to use the toolchain.
        javaHome = javaToolchains
            .launcherFor { languageVersion = JavaLanguageVersion.of(libs.versions.jdk.get()) }
            .get().metadata.installationPath.asFile.absolutePath
    }
}
