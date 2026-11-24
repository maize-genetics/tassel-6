# ADR-0022: BCF reading and writing, and BGZF writing

- Status: Accepted
- Date: 2026-11-17

## Context

TASSEL 6 reads VCF as plain text, as BGZF, and as one region of a CSI- or
TBI-indexed file (ADR-0021, epic tassel-6#78). It cannot read or write BCF,
and it cannot write VCF. `HapMap.write` is the only genotype writer. The
CLI and the app export HapMap only. Epic tassel-6#116 adds a BCF 2.2
reader, CSI region queries, and writers for VCF text, bgzipped VCF, and
BCF. This ADR is the decision those tasks implement. It supersedes
ADR-0021 as a whole.

ADR-0021 decided "Pure Kotlin, reading only" and listed writing bgzip,
building CSI or TBI, and BCF as out. An ADR is superseded only as a whole
(`docs/adr/README.md`). Every other ADR-0021 rule is restated below so a
later task does not re-open it. The rationale for those rules stays in
ADR-0021.

TASSEL 5 has no BCF reader and no BCF writer. htslib 1.25 is the BCF
oracle, as it already is for compressed and indexed VCF. VCF text export
is `dna/snp/ExportUtils.java` `writeToVCF`, from line 219, in the pinned
TASSEL 5 `v5.2.98` tree. That method writes `##fileformat=VCFv4.0`, the
`##Tassel` line when the reference allele is unknown, the `GT`, `AD`,
`DP`, `GQ`, and `PL` FORMAT lines, the `NS`, `DP`, and `AF` INFO lines,
then the `#CHROM` line. A record has QUAL `.`, FILTER `PASS`, and a `GT`
field joined with `/`. TASSEL 6 stores no reference allele, no depth, and
no taxon annotations (`docs/porting-ledger.md`, tassel-6#14), so the
writers emit the header and the REF/ALT rules `writeToVCF` produces for
that table. They do not invent `##SAMPLE` lines the model dropped. The
text oracle for those rules is the VCF-export goldens (tassel-6#119), not
a second copy of the Java branches in this ADR.

`Vcf.read` already accepts a phased diploid call. `VcfRecordParser.call`
ignores the separator between the two allele indexes, so `0|1` and `0/1`
become the same genotype. `VcfTest.readsUnphasedAndPhasedAndMissingCalls`
locks that. TASSEL 5 reads the characters at offsets 0 and 2 and also
ignores the separator. `Genotype` has no phase bit. Issue tassel-6#117
recommended refusing a phased genotype, on the premise that `Vcf.read`
rejects `|`. That premise is wrong, and the recommendation is not taken.
The reason is in the decision and under alternatives.

A BCF CSI carries no contig names (`l_aux = 0`).
`GenomicIndex.contigNames` rejects that shape, which is correct for VCF
and is the gap the BCF index path has to fill. Wheat's VCF header
declares `chr6D` and the VCF CSI omits it (`IndexedVcf.headerContigIds`).
tassel-6#118 records whether a BCF header and the sibling `.vcf.gz` index
list the same contigs. The UI rule below does not wait on that
measurement.

## Decision

### Carried forward from ADR-0021

These rules still hold. The paragraphs that replace "reading only" and
"BCF and bgzip writing are out" are in the next section.

- **Pure Kotlin, no new dependencies.** DEFLATE (RFC 1951), gzip
  (RFC 1952), BGZF (SAMv1 section 4.1), CSI v1, and TBI decoding stay in
  `tassel-core` `commonMain`. The implementation uses no zlib, libdeflate,
  htsjdk, or `java.util.zip`. A platform inflate (`java.util.zip`,
  `platform.zlib`) needs a spike with numbers first, the rule ADR-0020 set
  for BLAS, and it is not part of this epic. No new dependencies.
  `GenotypeTable`, the operation DSL, and `Engine` do not change.

- **Positional reads stay a source, not a path.** `RandomAccessSource` in
  `commonMain` has a `size: Long`, a suspending read at a `Long` position
  into a caller-owned `ByteArray`, and `close`. Concurrent reads are safe.
  The in-memory implementation lives in `commonMain`. File implementations
  live in `jvmMain`, in `nativeMain`, and in a Node-only source set. They
  do not live in `wasmJsMain`, because that source set is linked into the
  browser. The browser `File` / `Blob` implementation lives in
  `tassel-app`'s `wasmJsMain`. Readers never take a path. Sequential
  readers still take a `kotlinx.io.Source`. Plain sequential gzip stays a
  wrapping `RawSource` (ADR-0017).

- **Index discovery stays outside the core.** For a VCF, the CLI and the
  desktop app look beside the file for `.csi`, then `.tbi`. An explicit
  `--index` overrides that search. The browser has no sibling-file lookup,
  so the user supplies the index. Core readers take an already-opened
  index source.

- **Coordinates are `Long` until the table boundary.** BGZF, CSI, and TBI
  use `Long`. Values become `Position.position: Int` at the
  `GenotypeTable` boundary (ADR-0012). A position past `Int.MAX_VALUE` is
  an `IllegalArgumentException`. `Position` itself does not change.

- **The two oracles do not mix.** Decoded tutorial content stays
  `testdata/golden/` and TASSEL 5 `v5.2.98` (ADR-0009). Compressed and
  indexed behavior is htslib 1.25 in `testdata/htslib/`, with a `MANIFEST`
  of the version and the commands. Those files follow the same
  never-edit-by-hand rule and code-owner review as `testdata/golden/`.
  They are not TASSEL 5 goldens, and a change to them does not use the
  `golden-update` label. How a BCF fixture relates to a golden is below.

- **CRC32 and ISIZE are verified** for every gzip member and every BGZF
  block, as htslib does.

- **Parallelism is not a user-facing parameter.** Output is identical for
  any worker count. There is no thread parameter on the operation spec
  (ADR-0020). This applies to BGZF decoding and, below, to BGZF
  compression.

- **Region overlap is `record`.** A record is returned when its reference
  span overlaps the query. The span starts at `POS` and covers the longest
  of the `REF` allele, a valid `INFO/END`, the `INFO/SVLEN` of symbolic
  alleles such as `<DEL>`, and `FORMAT/LEN` for gVCF `<*>` alleles.
  `bcftools view -r` is the oracle. Query coordinates are 1-based and
  inclusive. `filter-sites` is unchanged: it keeps a site only when its
  `POS` is inside the requested range.

- **One region per query.** `chr3B`, `chr3B:1000-2000`, and `chr3B:1000-`
  mean the whole contig, an inclusive range, and from that position
  through the end of the contig. `chr3B:600000000` means that single
  position (`HTS_PARSE_ONE_COORD`). Commas that sit between digits are
  thousands separators. A comma that is not between digits is an error.
  `{name}:beg-end` is the form for a contig name that contains `:`.

- **Contig names match the index, then `Chromosome.of`.** The query is
  compared to the index names exactly, case-sensitive, first. If that
  misses, both sides are compared with `Chromosome.of`. One normalized
  hit is used. Zero hits is an unknown contig. Two or more is an error
  that names the ambiguous contigs. `Chromosome.of` itself does not
  change. For a BCF the names in this comparison are the header names, as
  the index section below requires.

- **A typo is an error. An empty range is an empty table.** An unknown
  contig is an error listing the contig names as written. An inverted
  region (`beg > end`) is an error. A known contig with nothing in range,
  or a region past the index's coordinate limit, returns a table with no
  sites.

- **A missing BGZF EOF marker is a warning.** The warning text is
  `EOF marker is absent. The input may be truncated` (`MISSING_BGZF_EOF`).
  The read succeeds. The core does not throw, and it does not grow
  `ProgressReporter`. The warning is carried on the read result, as
  `VcfRead.warnings` already does, and the caller shows it. The CLI prints
  it with `CliLog.warn`. A truncated block inside the file is an
  `IllegalArgumentException`.

- **`Vcf.read` stays the VCF reader.** The `Source` overload reads plain
  text. The `RandomAccessSource` overload reads a bgzipped VCF. An indexed
  VCF is `IndexedVcf.open`. None of the three detects BCF.

Reads still suspend at progress reports (ADR-0019).

### BCF reading and writing

- **Read parity.** A BCF read equals `Vcf.read` of the text `bcftools view`
  prints for that file: the same taxa, sites, alleles, and genotypes,
  including the quirks `VcfRecordParser` already has. Phase is dropped, as
  below. An allele index of 10 or more is unsupported: TASSEL 5 and
  `Vcf.read` read the characters at offsets 0 and 2, and the BCF reader
  returns that same result. A call that `Vcf.read` would reject is
  rejected. A call that `Vcf.read` would accept returns the same table.

- **Phase is dropped, not stored.** A phased diploid call is read as the
  unphased genotype of the same two alleles. `0|1` and `0/1` are the same
  `Genotype` afterwards. The reader does not throw. `Vcf.read` is left as
  it is. Refusing phase in the BCF reader would make that read differ from
  `Vcf.read` of the `bcftools view` text. Rejecting `|` inside `Vcf.read`
  would change a behavior TASSEL 5 and `VcfTest` already lock, and
  `Genotype` has no phase bit to store (ADR-0012). Writers emit unphased
  `GT` with `/`, which is the separator `writeToVCF` writes. Storing
  phase, round-tripping `|`, and adding a phase bit are out of this epic.

- **Framing.** BCF 2.2 inside BGZF is read sequentially and through a CSI.
  Raw BCF, the bytes `bgzip -dc` prints, with the magic at byte 0, is read
  sequentially only. It has no index. BCF 2.1 is refused. The error names
  htslib's reason: `only BCFv2.2 is supported`. The magic is the BCF 2.2
  magic from VCF 4.3 section 6, read from the start of the inflated BGZF
  payload or from byte 0 of a raw file.

- **Detection is one shared interface.** The CLI and the app open an
  indexed file through one interface in `tassel-core` `commonMain`. The
  interface reads the magic bytes above and chooses the VCF or BCF
  reader. It takes an already-opened data source and an already-opened
  index source. It does not search a filesystem. A sequential open uses
  that same magic test. `Vcf.read(RandomAccessSource)` and
  `IndexedVcf.open` stay VCF entry points. tassel-6#125 and tassel-6#126
  wire the CLI and the app to this interface. They do not add a second
  sniffer. VCF index discovery stays `.csi`, then `.tbi`, outside the
  core.

- **Header dictionaries.** The string dictionary is keyed by `IDX`. `IDX`
  values may have gaps. PASS is 0. A dictionary with no `IDX` numbers its
  strings in order of first appearance. The contig dictionary is separate
  from the string dictionary. Decoding trusts the typed descriptor on the
  record. It does not trust the declared `Type=` string in the header.

- **BCF indexes are CSI only.** The index is `<file>.csi`, or the file
  given by `--index`. Contig names come from the BCF header, because a
  BCF CSI has `l_aux = 0`. The VCF path keeps rejecting a CSI with no
  contig names. A TBI given for a BCF is an error. A CSI that carries
  tabix names, given for a BCF, is an error.

- **The BCF contig list is the header dictionary.** The UI shows that
  list, in header order. tassel-6#118 records whether it matches the
  sibling `.vcf.gz` index, including wheat `chr6D`. The rule is the same
  either way: the BCF list is the header, and the VCF list stays the index
  list (`IndexedVcf.contigs`).

- **Writers take a `Sink`.** `Vcf.write` writes VCF text. `Bcf.write`
  writes BCF 2.2 records. Both take a `kotlinx.io.Sink` (ADR-0017). A
  plain `.vcf` is `Vcf.write` on the file sink. A `.vcf.gz` is `Vcf.write`
  on a `BgzfWriter`. A `.bcf` is `Bcf.write` on a `BgzfWriter`, which is
  BCF inside BGZF. The BCF header is the header lines `writeToVCF` writes
  for this table, plus only what htslib adds when it encodes that VCF: the
  PASS filter line, `IDX=` on the dictionary lines, and the `##contig`
  lines. At BGZF stored blocks the header and the records are
  byte-identical to htslib for the same VCF. At the fixed compressed level
  below, tests compare inflated bytes and the BGZF block-size list.
  Compressed bytes are not required to match htslib.

- **BGZF writing.** `BgzfWriter` is a suspending writer over a `Sink`. It
  is not a `RawSink`. `write` does not suspend: wasm has no `runBlocking`.
  A caller's `flush()` does not end a block. `endHeader()` ends the header
  block. `close()` writes the 28-byte EOF block and leaves the underlying
  sink open. `tell()` reports the BGZF virtual offset for a later index
  writer. This epic does not write an index. The gzip header is pinned:
  MTIME 0, XFL 0, OS 255, XLEN 6, and the BC subfield. BCF blocks hold
  65,280 uncompressed bytes, and a record may straddle two blocks. VCF
  text blocks end at line breaks, with the header in its own block. The
  caller selects the BCF layout or the VCF text layout. There is no level
  argument and no thread argument (ADR-0020). tassel-6#128 may build the
  writer with stored blocks first. Users do not see stored-only output.
  The writer the CLI and the app ship uses one fixed DEFLATE level
  (tassel-6#131). The same bytes come out on every target and for any
  worker count.

- **`.vcf.gz` is BGZF.** TASSEL 5 writes a plain `GZIPOutputStream`
  (`Utils.getBufferedWriter`). BGZF is valid multi-member gzip, and it can
  be indexed. Parity for `.vcf.gz` is the inflated text against the VCF
  golden. It is not byte parity with TASSEL 5's gzip. When the writer
  lands, the porting ledger records `writeToVCF` of a `.vcf.gz` as a
  redesign. This decision does not add that row. The HapMap row for
  `ExportUtils` stays a port.

- **The `--output` suffix selects the writer.** The suffixes are
  `.hmp.txt`, `.vcf`, `.vcf.gz`, and `.bcf`. Matching is case-insensitive,
  as `Utils.addSuffixIfNeeded` lowercases the name before it compares
  suffixes. A name with none of those suffixes stays HapMap, which is
  what the CLI writes today. Standard output stays HapMap. There is no
  suffix for `.hmp.txt.gz`.

- **There is no `convert` operation.** Existing commands keep writing
  through `--output`. That is already how `filter-sites` converts a VCF
  input to HapMap (`VcfEndToEndTest`). tassel-6#133 adds the new suffixes
  to that path. It does not add a `convert` operation. `tassel translate`
  (ADR-0006) will map `-exportType VCF` onto the `.vcf` suffix. The
  translator is still the stretch row in `docs/porting-ledger.md`. This
  ADR records the mapping. It does not build the translator.

- **A golden and an htslib fixture stay different oracles.** A golden is
  TASSEL 5 text: the meaning of the table (ADR-0009). An htslib fixture is
  htslib 1.25's framing of a file: the bytes of BCF or BGZF. tassel-6#129
  is htslib's encoding of a TASSEL 5 VCF-export golden. The `MANIFEST`
  stores that golden's SHA-256 so the fixture can be regenerated when the
  golden changes. The fixture does not replace the golden, and editing the
  fixture does not use the `golden-update` label.

Still out of this epic: writing a CSI or a TBI, FORMAT fields other than
`GT` on a record, phase as stored data, BCF 2.1, and `.hmp.txt.gz`. Also
still out, carried forward from ADR-0021: remote sources, `.gzi` indexes,
and a platform inflate without a measured spike.

## Consequences

- tassel-6#120 through tassel-6#134 can add public API against these
  rules. tassel-6#118 and tassel-6#119 produce fixtures and goldens and do
  not wait on this ADR.
- The browser, the CLI, the desktop app, and the Native bindings share one
  BCF decoder and one BGZF writer. ADR-0008 is unchanged.
- `filter-sites` and `testdata/golden/` stay the TASSEL 5 contracts they
  are. An indexed BCF query follows the `record` overlap rule, so it will
  disagree with `filter-sites` on a record whose span overlaps the region
  and whose `POS` does not. That difference is ADR-0021, carried forward
  here.
- `.vcf.gz` from TASSEL 6 will not match TASSEL 5's `GZIPOutputStream`
  bytes. The inflated text will match the VCF golden. The ledger row for
  that writer is a redesign, added by the writer task.
- `tassel translate` has a recorded mapping for `-exportType VCF` and
  still has no implementation.
- Host core count does not change BGZF output bytes. A missing platform
  deflate is a performance question for a later spike, not a correctness
  question.
- This ADR is Proposed. The maintainer sets the status to Accepted.
  tassel-6#117 stays open until that happens.

## Alternatives considered

- **Refuse a phased genotype.** The recommendation in tassel-6#117,
  written on the premise that `Vcf.read` rejects `|`. `Vcf.read` accepts
  `0|1` and drops the phase, and `VcfTest` locks that. A BCF reader that
  threw would break the read-parity rule. Changing `Vcf.read` to reject
  `|` would change TASSEL 5 behavior. Rejected.
- **Detect BCF inside `Vcf.read(RandomAccessSource)` and
  `IndexedVcf.open`.** Those two would stop being VCF readers, and the
  CLI and the app would each grow a sniffer. One shared detector keeps the
  VCF entry points stable. tassel-6#125 and tassel-6#126 are wiring.
  Rejected.
- **Show the `.vcf.gz` index contig list for a BCF.** A BCF CSI has no
  names of its own. The header dictionary is the list the file actually
  declares. tassel-6#118 measures whether that list differs from the VCF
  index. It does not choose a different list. Rejected.
- **Add a `convert` operation.** `filter-sites` with `--output` already
  converts. A new operation would be a second way to do the same write.
  Rejected for this epic. tassel-6#133 builds one only if a later ADR
  chooses it.
- **Require compressed bytes to match `bgzip` or `bcftools`.** htslib may
  use libdeflate. Inflated bytes and the block-size list are the portable
  check. Level 0 remains the byte-identical check for header and records.
  Rejected for the compressed level.
- **Write `.vcf.gz` with a plain `GZIPOutputStream`, matching TASSEL 5
  bytes.** A single gzip member cannot be indexed. BGZF is valid
  multi-member gzip. Parity is the inflated text. Rejected.
- **Keep ADR-0021 and write a narrow ADR that only adds BCF.** ADRs are
  superseded only as a whole. A narrow ADR would leave "reading only" in
  force beside the new writers. Rejected.
