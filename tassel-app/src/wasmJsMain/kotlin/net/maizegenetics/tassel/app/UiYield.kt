package net.maizegenetics.tassel.app

import kotlinx.coroutines.delay

/**
 * `yield` on wasm resumes as a microtask, and Chrome does not run `requestAnimationFrame` until a
 * macrotask. `delay` schedules a timer, which is that macrotask (ADR-0019: the reporter decides
 * what a progress report costs). One millisecond is enough to reach the next frame without
 * measuring a long sleep.
 */
internal actual suspend fun yieldToEventLoop() {
    delay(1)
}
