plugins {
    alias(libs.plugins.kotlin.jvm)
    application
}

kotlin {
    jvmToolchain(libs.versions.jdk.get().toInt())
}

dependencies {
    implementation(project(":tassel-core"))
    implementation(libs.clikt)
    testImplementation(libs.kotlin.test)
}

application {
    applicationName = "tassel"
    mainClass = "net.maizegenetics.tassel.cli.MainKt"
    // Mordant loads JNA for terminal detection; JDK 25 warns unless native access is granted.
    applicationDefaultJvmArgs = listOf("--enable-native-access=ALL-UNNAMED")
}
