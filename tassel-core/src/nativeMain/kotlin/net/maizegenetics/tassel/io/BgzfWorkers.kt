package net.maizegenetics.tassel.io

import kotlinx.cinterop.ExperimentalForeignApi
import platform.posix._SC_NPROCESSORS_ONLN
import platform.posix.sysconf

@OptIn(ExperimentalForeignApi::class)
actual fun bgzfWorkerCount(): Int {
    val processors = sysconf(_SC_NPROCESSORS_ONLN)
    return if (processors > 0) processors.toInt() else 1
}
