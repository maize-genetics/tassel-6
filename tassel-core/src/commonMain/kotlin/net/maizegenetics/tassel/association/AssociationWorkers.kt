package net.maizegenetics.tassel.association

/**
 * How many workers association kernels split markers across. One on wasmJs, where the dispatcher is a single thread;
 * the processor count on the JVM and on Native.
 *
 * A sibling of `kinshipWorkerCount` (ADR-0020): `association` does not call into `distance`. Results do not depend on
 * this count (ADR-0025): each marker is fitted by one worker from inputs no other marker changes, and results are put
 * back in table order. It is not an operation-spec parameter.
 */
internal expect fun associationWorkerCount(): Int
