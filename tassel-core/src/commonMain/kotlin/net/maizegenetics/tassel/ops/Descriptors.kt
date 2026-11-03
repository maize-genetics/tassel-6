package net.maizegenetics.tassel.ops

import kotlinx.serialization.SerialName
import kotlinx.serialization.Serializable

/**
 * Serializable description of an [Operation]: what CLI generation, UI forms, MCP tool schemas, and a
 * remote engine need to know about it, without its code.
 *
 * @property id operation id.
 * @property legacyPlugin TASSEL 5 plugin class simple name, or `null`.
 * @property help one-line description.
 * @property inputType simple class name of the input (`Unit` when there is none).
 * @property outputType simple class name of the output.
 * @property parameters parameters in declaration order.
 */
@Serializable
data class OperationDescriptor(
    val id: String,
    val legacyPlugin: String? = null,
    val help: String = "",
    val inputType: String,
    val outputType: String,
    val parameters: List<ParameterDescriptor> = emptyList(),
)

/**
 * Serializable description of a [Parameter]. The subtype gives the value type; in JSON it is the
 * `type` field (`int`, `double`, `boolean`, `enum`, or `file`).
 *
 * @see net.maizegenetics.plugindef.PluginParameter
 */
@Serializable
sealed interface ParameterDescriptor {
    /** Parameter name. */
    val name: String

    /** One-line help text. */
    val help: String

    /** TASSEL 5 command-line name, or `null`. */
    val legacy: String?

    /** Whether a value must be supplied. */
    val required: Boolean
}

/** An integer parameter. [min] and [max] are inclusive bounds, if any. */
@Serializable
@SerialName("int")
data class IntParameterDescriptor(
    override val name: String,
    override val help: String = "",
    override val legacy: String? = null,
    override val required: Boolean = false,
    val default: Int? = null,
    val min: Int? = null,
    val max: Int? = null,
) : ParameterDescriptor

/** A double parameter. [min] and [max] are inclusive bounds, if any. */
@Serializable
@SerialName("double")
data class DoubleParameterDescriptor(
    override val name: String,
    override val help: String = "",
    override val legacy: String? = null,
    override val required: Boolean = false,
    val default: Double? = null,
    val min: Double? = null,
    val max: Double? = null,
) : ParameterDescriptor

/** A boolean parameter. */
@Serializable
@SerialName("boolean")
data class BooleanParameterDescriptor(
    override val name: String,
    override val help: String = "",
    override val legacy: String? = null,
    override val required: Boolean = false,
    val default: Boolean? = null,
) : ParameterDescriptor

/** An enum parameter. [default] and [choices] are constant names, in declaration order. */
@Serializable
@SerialName("enum")
data class EnumParameterDescriptor(
    override val name: String,
    override val help: String = "",
    override val legacy: String? = null,
    override val required: Boolean = false,
    val default: String? = null,
    val choices: List<String>,
) : ParameterDescriptor

/** A file or directory path parameter. */
@Serializable
@SerialName("file")
data class FileParameterDescriptor(
    override val name: String,
    override val help: String = "",
    override val legacy: String? = null,
    override val required: Boolean = false,
    val default: String? = null,
    val kind: FileKind = FileKind.Input,
) : ParameterDescriptor
