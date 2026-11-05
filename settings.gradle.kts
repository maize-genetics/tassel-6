rootProject.name = "tassel-6"

pluginManagement {
    repositories {
        gradlePluginPortal()
        mavenCentral()
        google()
    }
}

plugins {
    // Provisions the JDK 25 toolchain (ADR-0011) when it is not installed locally.
    id("org.gradle.toolchains.foojay-resolver-convention") version "1.0.0"
}

dependencyResolutionManagement {
    repositories {
        mavenCentral()
        google()
    }
}

include(":tassel-core", ":tassel-cli", ":tassel-app")
