package net.maizegenetics.tassel.app

/**
 * Gives the UI a turn at the event loop after a progress report.
 *
 * Desktop `yield`s. The browser delays on a timer so the page can paint: a coroutine `yield` there
 * stays in the same turn, and kinship's reports never became frames (tassel-6#70).
 */
internal expect suspend fun yieldToEventLoop()
