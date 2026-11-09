package net.maizegenetics.tassel.app

/**
 * One line on the target's log, for a measurement the UI does not draw.
 *
 * Desktop prints it. The browser logs it. A region read uses this for the time and the bytes it actually fetched,
 * which is what the handoff quotes for a file larger than memory.
 */
internal expect fun logAppMessage(line: String)
