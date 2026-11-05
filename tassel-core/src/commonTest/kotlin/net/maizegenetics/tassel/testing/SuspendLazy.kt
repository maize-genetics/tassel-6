package net.maizegenetics.tassel.testing

/**
 * A value produced once by a suspending [produce], then handed out to every later caller.
 *
 * Reading the tutorial file and running an analysis over it both suspend (ADR-0019), so the shared fixtures that
 * used to be `by lazy` cannot be any more: `lazy` takes a plain lambda, and `commonTest` has no `runBlocking` to
 * bridge with on Node or Native. This is the same "parse once, reuse everywhere" the fixtures had, in the one
 * shape that works on every target.
 *
 * Not thread safe, and does not need to be: tests on every target run these one at a time.
 */
class SuspendLazy<T : Any>(private val produce: suspend () -> T) {
    private var value: T? = null

    suspend operator fun invoke(): T = value ?: produce().also { value = it }
}
