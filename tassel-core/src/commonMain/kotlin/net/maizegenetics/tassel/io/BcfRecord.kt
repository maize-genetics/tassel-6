package net.maizegenetics.tassel.io

import kotlinx.io.Source
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position

/**
 * One BCF 2.2 record from an inflated stream: the site and the calls [VcfRecordParser] would produce from the
 * line htslib prints.
 *
 * Returns null when [source] is exhausted before the 8-byte length prefix. A prefix whose body is short throws
 * [TRUNCATED_BCF_RECORD]. [Bcf.read] owns BGZF framing and the multi-record loop.
 *
 * The genotype block is turned back into the ASCII htslib would print, including haploid `0`, `10/0`, and
 * end-of-vector padding, and then [call] applies the VCF rules. Phase is a `|` in that text. [call] drops it.
 *
 * The layout is VCF 4.3 section 6, as htslib 1.25 writes it. `rlen` and `QUAL` are not stored. A missing `QUAL` is
 * the float bits [BCF_FLOAT_MISSING] and is not decoded as a [Float].
 */
internal fun readBcfRecord(source: Source, header: BcfHeader, calls: ByteArray): Position? {
    if (!source.request(1)) return null
    if (!source.request(BCF_RECORD_PREFIX_BYTES.toLong())) throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
    val lengths = source.readByteArray(BCF_RECORD_PREFIX_BYTES)
    val sharedLength = lengths.uint32Le(0)
    val indivLength = lengths.uint32Le(BCF_RECORD_PREFIX_BYTES / 2)
    val bodyLength = bcfBodyLength(sharedLength, indivLength)
    if (!source.request(bodyLength.toLong())) throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
    return decodeRecord(source.readByteArray(bodyLength), sharedLength.toInt(), header, calls)
}

/**
 * Shared plus individual body length.
 *
 * The same checks [readBcfRecord] applies to the 8-byte prefix, so a whole-file read can pull the body before
 * it calls the decoder. A length that is not a record throws [TRUNCATED_BCF_RECORD].
 */
internal fun bcfBodyLength(shared: Long, indiv: Long): Int {
    if (shared > Int.MAX_VALUE.toLong() || indiv > Int.MAX_VALUE.toLong()) {
        throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
    }
    val body = shared + indiv
    if (shared < SHARED_PREFIX || body > Int.MAX_VALUE.toLong()) {
        throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
    }
    return body.toInt()
}

/** Missing and end-of-vector, promoted from every integer width the way htslib promotes them. */
internal const val BCF_INT_MISSING = Int.MIN_VALUE

/** The value htslib prints nothing for. In a genotype it ends that sample. */
internal const val BCF_INT_VECTOR_END = Int.MIN_VALUE + 1

/** `QUAL` `.`, and the missing value of a float vector. Compared as bits, not as a [Float]. */
internal const val BCF_FLOAT_MISSING = 0x7F800001

/** End-of-vector for a float vector. */
internal const val BCF_FLOAT_VECTOR_END = 0x7F800002

/** A length prefix whose body is not all there. */
internal const val TRUNCATED_BCF_RECORD = "BCF record is truncated"

/**
 * A cursor over one typed-value region.
 *
 * The high nibble of the descriptor is the element count. Fifteen means the next typed integer is the real count.
 * [readInts], [readFloatBits], [readString], and [skip] consume one value. A FORMAT value's count is per sample, so
 * the genotype block passes the sample count and the payload is that many times longer.
 */
internal class BcfCursor(private val bytes: ByteArray, start: Int = 0, private val end: Int = bytes.size) {
    var offset: Int = start
        private set

    /** Integers in this vector, with each width's sentinels promoted to [BCF_INT_MISSING] and [BCF_INT_VECTOR_END]. */
    fun readInts(times: Int = 1): IntArray {
        val typed = descriptor()
        val count = elementCount(typed.count, times)
        if (typed.type == TYPE_NULL) {
            if (count != 0) throw IllegalArgumentException("BCF integer has type ${typed.type}")
            return IntArray(0)
        }
        if (typed.type != TYPE_INT8 && typed.type != TYPE_INT16 && typed.type != TYPE_INT32) {
            throw IllegalArgumentException("BCF integer has type ${typed.type}")
        }
        return IntArray(count) { readInt(typed.type) }
    }

    /** Raw float bits, including [BCF_FLOAT_MISSING] and [BCF_FLOAT_VECTOR_END]. */
    fun readFloatBits(times: Int = 1): IntArray {
        val typed = descriptor()
        if (typed.type != TYPE_FLOAT) throw IllegalArgumentException("BCF float has type ${typed.type}")
        return IntArray(elementCount(typed.count, times)) { int32() }
    }

    /** A character vector. Empty, which is how a missing ID is stored, returns an empty string. */
    fun readString(): String {
        val typed = descriptor()
        if (typed.count == 0) {
            if (typed.type != TYPE_NULL && typed.type != TYPE_CHAR) {
                throw IllegalArgumentException("BCF string has type ${typed.type}")
            }
            return ""
        }
        if (typed.type != TYPE_CHAR) throw IllegalArgumentException("BCF string has type ${typed.type}")
        val start = offset
        advance(typed.count.toLong())
        return bytes.decodeToString(start, offset)
    }

    /** Moves past one typed value without interpreting its payload. [times] is the FORMAT sample count. */
    fun skip(times: Int = 1) {
        val typed = descriptor()
        advance(elementCount(typed.count, times).toLong() * width(typed.type))
    }

    private fun descriptor(): Typed {
        val bits = uint8()
        val type = bits and 0x0f
        if (type == TYPE_NULL) return Typed(TYPE_NULL, 0)
        val count = bits ushr 4
        return Typed(type, if (count == OVERFLOW_COUNT) overflowCount() else count)
    }

    private fun overflowCount(): Int {
        val bits = uint8()
        val type = bits and 0x0f
        val count = bits ushr 4
        if (count == OVERFLOW_COUNT) throw IllegalArgumentException("BCF typed size is nested")
        if (count != 1) throw IllegalArgumentException("BCF typed size must be one integer")
        val value = readInt(type)
        if (value < 0) throw IllegalArgumentException("BCF typed value has a negative size")
        return value
    }

    private fun readInt(type: Int): Int = when (type) {
        TYPE_INT8 -> when (val value = int8()) {
            INT8_MISSING -> BCF_INT_MISSING
            INT8_VECTOR_END -> BCF_INT_VECTOR_END
            else -> value
        }
        TYPE_INT16 -> when (val value = int16()) {
            INT16_MISSING -> BCF_INT_MISSING
            INT16_VECTOR_END -> BCF_INT_VECTOR_END
            else -> value
        }
        TYPE_INT32 -> int32()
        else -> throw IllegalArgumentException("BCF integer has type $type")
    }

    private fun uint8(): Int {
        if (offset >= end) throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
        val value = bytes[offset].toInt() and 0xff
        offset += 1
        return value
    }

    private fun int8(): Int {
        val value = uint8()
        return if (value < 0x80) value else value - 0x100
    }

    private fun int16(): Int {
        val unsigned = uint8() or (uint8() shl 8)
        return (unsigned shl 16) shr 16
    }

    private fun int32(): Int = uint8() or (uint8() shl 8) or (uint8() shl 16) or (uint8() shl 24)

    private fun advance(count: Long) {
        if (count < 0 || offset.toLong() + count > end.toLong()) throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
        offset += count.toInt()
    }

    private class Typed(val type: Int, val count: Int)
}

private fun decodeRecord(body: ByteArray, sharedLength: Int, header: BcfHeader, calls: ByteArray): Position {
    val chromosome = contigName(body, header)
    val pos1 = body.int32Le(POS_OFFSET).toLong() + 1L
    val alleleInfo = body.uint32Le(ALLELE_INFO_OFFSET)
    val nInfo = (alleleInfo and 0xffffL).toInt()
    val nAllele = (alleleInfo ushr 16).toInt()
    val sampleInfo = body.uint32Le(SAMPLE_OFFSET)
    val nSample = (sampleInfo and 0xffffffL).toInt()
    val nFormat = (sampleInfo ushr 24).toInt()

    val shared = BcfCursor(body, SHARED_PREFIX, sharedLength)
    val id = shared.readString()
    val snpId = id.takeUnless { it.isEmpty() || it == "." }
    val site = snpId ?: "$chromosome:$pos1"
    val position = requireVcfPosition(pos1, site)
    if (nAllele <= 0) throw IllegalArgumentException("VCF site $site has no REF or ALT allele: []")
    val alleleText = Array(nAllele) { shared.readString() }
    shared.skip()
    repeat(nInfo) {
        shared.skip()
        shared.skip()
    }
    if (shared.offset != sharedLength) {
        throw IllegalArgumentException("BCF record shared block does not match l_shared")
    }
    if (nSample != calls.size) {
        throw IllegalArgumentException("VCF site $site has $nSample genotypes, expected ${calls.size}")
    }

    val reference = alleleText[0]
    val alternate = if (nAllele == 1) "." else alleleText.copyOfRange(1, nAllele).joinToString(",")
    val variants = splitVariants(variantsOf(reference, alternate), site)
    val alleles = allelesOf(variants, site)
    val indiv = BcfCursor(body, sharedLength, body.size)
    var sawGt = false
    repeat(nFormat) {
        val key = indiv.readInts()
        if (key.size != 1) throw IllegalArgumentException("BCF FORMAT key must be one integer")
        if (key[0] == header.gtId && !sawGt) {
            // The descriptor count is per sample. Values are grouped by sample, not by allele.
            val values = indiv.readInts(nSample)
            if (nSample > 0) {
                val perSample = values.size / nSample
                for (sample in 0 until nSample) {
                    val start = sample * perSample
                    val text = genotypeText(values, start, start + perSample)
                    calls[sample] = call(text, 0, text.size, alleles, site)
                }
                sawGt = true
            }
        } else {
            indiv.skip(nSample)
        }
    }
    if (nSample > 0 && !sawGt) throw IllegalArgumentException("VCF site $site has no GT field in FORMAT")
    if (indiv.offset != body.size) {
        throw IllegalArgumentException("BCF record genotype block does not match l_indiv")
    }
    return Position(
        chromosome = Chromosome.of(chromosome),
        position = position,
        snpId = snpId,
        knownVariants = variants,
    )
}

/**
 * The ASCII htslib prints for one sample.
 *
 * End-of-vector stops the sample, so a haploid call stored in a diploid slot is `0` and not `0/.`. The separator
 * before every later allele is `|` when that value's low bit is set. [call] then reads the text.
 */
private fun genotypeText(values: IntArray, from: Int, until: Int): ByteArray {
    val text = StringBuilder()
    var written = 0
    var index = from
    while (index < until) {
        val value = values[index]
        if (value == BCF_INT_VECTOR_END) break
        if (written > 0) text.append(if ((value and 1) != 0) '|' else '/')
        if (value == BCF_INT_MISSING || value == 0) {
            text.append('.')
        } else {
            text.append((value shr 1) - 1)
        }
        written += 1
        index += 1
    }
    return text.toString().encodeToByteArray()
}

private fun contigName(body: ByteArray, header: BcfHeader): String {
    val index = body.int32Le(CHROM_OFFSET)
    if (index < 0 || index >= header.contigs.size) {
        throw IllegalArgumentException("BCF record names a contig the header does not list: $index")
    }
    val contig = header.contigs[index]
        ?: throw IllegalArgumentException("BCF record names a contig the header does not list: $index")
    return contig.id
}

private fun elementCount(count: Int, times: Int): Int {
    val product = count.toLong() * times.toLong()
    if (times < 0 || product > Int.MAX_VALUE) throw IllegalArgumentException("BCF typed value is too long")
    return product.toInt()
}

private fun width(type: Int): Int = when (type) {
    TYPE_NULL -> 0
    TYPE_INT8, TYPE_CHAR -> 1
    TYPE_INT16 -> 2
    TYPE_INT32, TYPE_FLOAT -> 4
    else -> throw IllegalArgumentException("BCF typed value has type $type")
}

private fun ByteArray.int32Le(at: Int): Int = uint32Le(at).toInt()

/** `l_shared` and `l_indiv`, the 8 bytes in front of every record body. */
internal const val BCF_RECORD_PREFIX_BYTES = 8
private const val CHROM_OFFSET = 0
private const val POS_OFFSET = 4
private const val ALLELE_INFO_OFFSET = 16
private const val SAMPLE_OFFSET = 20
private const val SHARED_PREFIX = 24

private const val TYPE_NULL = 0
private const val TYPE_INT8 = 1
private const val TYPE_INT16 = 2
private const val TYPE_INT32 = 3
private const val TYPE_FLOAT = 5
private const val TYPE_CHAR = 7
private const val OVERFLOW_COUNT = 15

private const val INT8_MISSING = -128
private const val INT8_VECTOR_END = -127
private const val INT16_MISSING = -32768
private const val INT16_VECTOR_END = -32767
