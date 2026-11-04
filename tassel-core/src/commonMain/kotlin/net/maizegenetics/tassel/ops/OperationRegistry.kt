package net.maizegenetics.tassel.ops

/**
 * An immutable set of operations, looked up by id or by TASSEL 5 names. Front ends (CLI, UI, MCP,
 * legacy translator) are driven from a registry.
 *
 * @param operations the operations; ids and legacy plugin names must be unique.
 * @throws IllegalArgumentException on a duplicate id or legacy plugin name.
 */
class OperationRegistry(operations: Iterable<Operation<*, *>>) {
    /** Registry of [operations]. */
    constructor(vararg operations: Operation<*, *>) : this(operations.asList())

    /** All operations, in registration order. */
    val operations: List<Operation<*, *>> = operations.toList()

    private val byId = this.operations.associateBy { it.id }
    private val byLegacyPlugin = this.operations.filter { it.legacyPlugin != null }.associateBy { it.legacyPlugin!! }

    init {
        require(byId.size == this.operations.size) { "duplicate operation id in ${this.operations.map { it.id }}" }
        require(byLegacyPlugin.size == this.operations.count { it.legacyPlugin != null }) {
            "duplicate legacy plugin in ${this.operations.mapNotNull { it.legacyPlugin }}"
        }
    }

    /** Descriptors of all operations, in registration order. */
    val descriptors: List<OperationDescriptor> get() = operations.map { it.descriptor }

    /** The operation with this [id], or `null`. */
    operator fun get(id: String): Operation<*, *>? = byId[id]

    /** The operation replacing the TASSEL 5 plugin with this class simple name, or `null`. */
    fun byLegacyPlugin(legacyPlugin: String): Operation<*, *>? = byLegacyPlugin[legacyPlugin]

    /**
     * The parameter that TASSEL 5 plugin [legacyPlugin] called [legacyName] (as in
     * `run_pipeline.pl -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05`), or `null`.
     */
    fun parameterByLegacyName(legacyPlugin: String, legacyName: String): Parameter<*>? =
        byLegacyPlugin(legacyPlugin)?.parameterByLegacyName(legacyName)
}
