package net.maizegenetics.tassel.spike.gram

import java.nio.file.Files
import java.nio.file.Path
import kotlin.test.Test
import kotlin.test.assertTrue

class ShaderFileTest {
    @Test
    fun gramShaderDeclaresBothEntryPoints() {
        val source = Files.readString(shaderFile())
        assertTrue(source.contains("// @@NAIVE"))
        assertTrue(source.contains("// @@TILED"))
        assertTrue(source.contains("fn gram_naive"))
        assertTrue(source.contains("fn gram_tiled"))
        assertTrue(source.contains("workgroupBarrier"))
    }

    private fun shaderFile(): Path {
        val candidates = listOf(
            Path.of("wgsl/gram_f32.wgsl"),
            Path.of("spikes/gpu-gram/wgsl/gram_f32.wgsl"),
        )
        return candidates.firstOrNull { Files.exists(it) }
            ?: error("gram_f32.wgsl not found from ${Path.of("").toAbsolutePath()}")
    }
}
