package net.maizegenetics.tassel.app.state

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.ops.OperationDescriptor

/** `OperationDescriptor.inputType` of an operation that takes several selected datasets (ADR-0024). */
internal const val DATASET_INPUT_TYPE: String = "DataSet"

/**
 * The operations that can run on [selection].
 *
 * One dataset offers the operations whose input type is that dataset's type. Two or more offer the `DataSet`
 * operations, which are the joins. Nothing selected offers nothing. The match is on the descriptor's input type
 * alone, so an operation added to the registry shows up with no change here.
 */
fun operationsFor(operations: List<OperationDescriptor>, selection: List<Dataset>): List<OperationDescriptor> =
    when (selection.size) {
        0 -> emptyList()
        1 -> operations.filter { it.inputType == selection[0].inputType }
        else -> operations.filter { it.inputType == DATASET_INPUT_TYPE }
    }

/** The [DataSet] a multi-selection becomes, items in selection order. */
fun dataSetOf(selection: List<Dataset>): DataSet = DataSet(selection.map { it.data.value })
