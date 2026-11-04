package net.maizegenetics.tassel.ops

import net.maizegenetics.tassel.filter.filterSitesOperation
import net.maizegenetics.tassel.filter.filterTaxaOperation

/**
 * Every operation TASSEL 6 ships. The CLI, the app, the MCP tools, and the legacy translator are
 * all built from [registry] (ADR-0002 and ADR-0006); none of them should assemble its own.
 *
 * Adding an operation is one line in [registry], so the parallel porting PRs conflict on that line
 * and nothing else. Keep the list in the order front ends should present the operations.
 */
object TasselOperations {
    /** The operations, in presentation order. */
    val registry: OperationRegistry =
        OperationRegistry(
            // summaryOperation,     // #15 GenotypeSummaryPlugin
            filterSitesOperation, // #16 FilterSiteBuilderPlugin
            filterTaxaOperation, // #16 FilterTaxaBuilderPlugin
            // kinshipOperation,     // #18 KinshipPlugin
            // pcaOperation,         // #19 PrincipalComponentsPlugin
        )
}
