plugins {
    alias(libs.plugins.kotlin.multiplatform) apply false
    alias(libs.plugins.kotlin.jvm) apply false
    alias(libs.plugins.kotlin.compose) apply false
    alias(libs.plugins.compose.multiplatform) apply false
    alias(libs.plugins.ktlint)
    alias(libs.plugins.kover)
}

dependencies {
    kover(project(":tassel-core"))
    kover(project(":tassel-cli"))
    kover(project(":tassel-app"))
}
