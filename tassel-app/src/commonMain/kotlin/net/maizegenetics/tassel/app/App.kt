package net.maizegenetics.tassel.app

import androidx.compose.runtime.Composable
import androidx.compose.runtime.collectAsState
import androidx.compose.runtime.getValue
import androidx.compose.runtime.remember
import androidx.compose.runtime.rememberCoroutineScope
import net.maizegenetics.tassel.app.data.GenotypeFiles
import net.maizegenetics.tassel.app.state.AppViewModel
import net.maizegenetics.tassel.app.theme.TasselTheme
import net.maizegenetics.tassel.app.ui.AppScreen
import kotlin.coroutines.CoroutineContext
import kotlin.coroutines.EmptyCoroutineContext

/**
 * The whole app. Both entry points call this with their own [GenotypeFiles], which is the only thing that differs
 * between desktop and the browser.
 *
 * @param files where genotype files come from on this target.
 * @param readContext the context HapMap parsing runs on; see [AppViewModel].
 */
@Composable
fun App(files: GenotypeFiles, readContext: CoroutineContext = EmptyCoroutineContext) {
    val scope = rememberCoroutineScope()
    val viewModel = remember(files, readContext) { AppViewModel(files, scope, readContext) }
    val state by viewModel.state.collectAsState()
    TasselTheme {
        AppScreen(state = state, onEvent = viewModel::handle)
    }
}
