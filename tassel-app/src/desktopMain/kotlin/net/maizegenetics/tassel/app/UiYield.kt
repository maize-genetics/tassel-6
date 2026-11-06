package net.maizegenetics.tassel.app

import kotlinx.coroutines.yield

internal actual suspend fun yieldToEventLoop() {
    yield()
}
