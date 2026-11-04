package net.maizegenetics.tassel.app.theme

import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.darkColorScheme
import androidx.compose.material3.lightColorScheme
import androidx.compose.runtime.Composable
import androidx.compose.runtime.CompositionLocalProvider
import androidx.compose.runtime.staticCompositionLocalOf
import androidx.compose.ui.graphics.Color

/**
 * The app theme: Material 3 plus the genotype palette. Every color the UI draws comes from here, either from
 * `MaterialTheme.colorScheme` or from [LocalAlleleColors] (`.cursor/rules/compose-ui.mdc`).
 */
@Composable
fun TasselTheme(dark: Boolean = false, content: @Composable () -> Unit) {
    CompositionLocalProvider(LocalAlleleColors provides if (dark) DarkAlleleColors else LightAlleleColors) {
        MaterialTheme(
            colorScheme = if (dark) darkColorScheme() else lightColorScheme(),
            content = content,
        )
    }
}

/** The genotype palette for the surrounding [TasselTheme]. */
val LocalAlleleColors = staticCompositionLocalOf { LightAlleleColors }

/**
 * Nucleotide colors for the light scheme. The four bases keep the hues TASSEL 5's genotype viewer and most other
 * alignment browsers use (A green, C blue, G orange, T red) so a maize geneticist can read the grid without a legend.
 */
val LightAlleleColors: AlleleColors = AlleleColors(
    a = Color(0xFF2E7D32),
    c = Color(0xFF1565C0),
    g = Color(0xFFEF6C00),
    t = Color(0xFFC62828),
    insert = Color(0xFF6A1B9A),
    gap = Color(0xFF4E342E),
    heterozygous = Color(0xFFFDD835),
    missing = Color(0xFFE0E0E0),
)

/** The same hues lightened for the dark scheme, where they sit on a dark surface. */
val DarkAlleleColors: AlleleColors = AlleleColors(
    a = Color(0xFF81C784),
    c = Color(0xFF64B5F6),
    g = Color(0xFFFFB74D),
    t = Color(0xFFE57373),
    insert = Color(0xFFBA68C8),
    gap = Color(0xFFA1887F),
    heterozygous = Color(0xFFFFF176),
    missing = Color(0xFF424242),
)
