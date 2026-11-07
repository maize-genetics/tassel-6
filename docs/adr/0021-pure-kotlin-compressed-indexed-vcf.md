# ADR-0021: Pure-Kotlin compressed and indexed VCF reading

- Status: Accepted
- Date: 2026-11-07

## Context

TASSEL 6 reads VCF only as uncompressed text, start to finish, into an
in-memory `GenotypeTable` (`io.Vcf.read`, tassel-6#14). That serves the
tutorial file. It does not serve a wheat panel: a bgzip-compressed VCF of
tens to hundreds of gigabytes, indexed with CSI, where the question is one
region of one chromosome. Epic tassel-6#78 builds that read path. This ADR
fixes the choices the tasks with public API or user-visible behavior wait
on (tassel-6#86, #88, #90, #91, #93). It records the decision asked for in
tassel-6#79.

TBI and BAI use a fixed binning scheme (`min_shift` 14, `depth` 5) that
cannot represent a coordinate at or past 2^29 = 536,870,912. Most wheat
chromosomes are longer than that; 3B is about 830 Mbp. CSI stores
`min_shift` and `depth` in the file, so a coordinate past 512 Mbp can be
indexed. Both index formats are themselves BGZF-compressed and point into
the data with BGZF virtual offsets, so the stack rests on a DEFLATE decoder.

TASSEL 5 opens `.gz` text through `GZIPInputStream` in
`util/Utils.getBufferedReader`, chosen by the file extension. That also
decodes BGZF, because BGZF is valid multi-member gzip. Its only index is a
`.lix` file of BGZF virtual offsets every 10 lines of a bgzipped HapMap
(`dna/snp/io/LineIndexBuilder`, `BuilderFromHapMapLIX`, through htsjdk).
That index finds lines by number. TASSEL 5 has no region queries and no CSI
or TBI support. `dna/snp/io/BuilderFromVCF` remains the reader whose
genotype semantics stay; see `io.Vcf`.

The roadmap's Phase 3 put "bgzip- and tabix-indexed VCF" on htsjdk, on the
JVM only, and left "Gzip handling per platform" undecided. ADR-0001 keeps
`commonMain` free of JVM-only dependencies, and ADR-0015 puts readers for
the lab's pangenome VCFs in `tassel-core` so the browser has them. A
JVM-only htsjdk path would leave the browser and the Native bindings
without the reader.

ADR-0017 still holds and is extended here, not superseded. Readers take a
`kotlinx.io.Source` and never a path, and compression for a sequential read
stays a wrapping `RawSource`. kotlinx-io 0.9.2 has no compression and no
positional reads. Its consequence that "gzip can be added for the JVM first"
is narrowed below: the compressed and indexed read path is pure Kotlin on
every target. ADR-0012 and `GenotypeTable` are unchanged. Region reads
return in-memory tables. `Position.position` stays `Int`. ADR-0019 still
requires reads to suspend at progress reports. ADR-0020 is the precedent
for two rules reused here: a platform fast path is spiked first and adopted
only on a measured win, and results must not depend on the host's core
count. ADR-0009 is unchanged: `testdata/golden/` stays the TASSEL 5 parity
oracle.

htslib 1.25 is the oracle for compressed and indexed behavior, and its own
tools disagree in places this ADR has to pick. `tabix` reads
`chr3B:600000000` as "to the end of the contig", while `bcftools view -r`
reads it as one position (`HTS_PARSE_ONE_COORD`). `tabix` treats commas in
`chr3B:599,999,999-600,000,001` as thousands separators
(`HTS_PARSE_THOUSANDS_SEP`); `bcftools view -r` treats them as separators
in a list of regions and returns nothing. The default `record` overlap and
the `pos` overlap give different answers for a deletion that starts before
the region. Checked with htslib 1.25, `bcftools view -H -r
chr3B:536870002-536870010` returns the record at 536,870,000 (`ACGT` to
`A`) under `--regions-overlap record` and returns nothing under
`--regions-overlap pos`. htslib also prints nothing and exits 0 for an
unknown contig, an inverted region, and a region past the index's
coordinate limit.

## Decision

- **Pure Kotlin, reading only.** DEFLATE (RFC 1951), gzip (RFC 1952), BGZF
  (SAMv1 section 4.1), CSI v1, and TBI decoding all live in `tassel-core`
  `commonMain`. The read path uses no zlib, libdeflate, htsjdk, or
  `java.util.zip`. A platform inflate (`java.util.zip`, `platform.zlib`)
  needs a spike with numbers first, the same rule ADR-0020 set for BLAS,
  and it is not part of this epic. Writing bgzip, building CSI or TBI,
  BCF, `.gzi` indexes, and remote sources are out. No new dependencies.
  `.vcf.gz`, `.csi`, and `.tbi` files come from htslib.

- **Positional reads are a new source, not a path.** A `RandomAccessSource`
  in `commonMain` has three members: a `size: Long`, a `suspend` read at a
  `Long` position into a caller-owned `ByteArray` (the caller names the
  destination range inside that array), and `close`. The read suspends
  because the browser's `Blob` reads are asynchronous. Concurrent reads
  must be safe, because parallel decoding needs them. An in-memory
  implementation lives in `commonMain`. File implementations live in
  `jvmMain`, in `nativeMain` (macOS arm64 and Linux x64), and in a
  Node-only source set. They do not live in `wasmJsMain`: that source set
  is compiled for both `browser()` and `nodejs()`, and `tassel-app` links
  it into the browser, so a Node `fs` implementation there would ship in
  the browser bundle. The browser `File` / `Blob` implementation lives in
  `tassel-app`'s `wasmJsMain`. Readers themselves still never take a path.
  This extends ADR-0017. It does not supersede it. Sequential readers still
  take a `kotlinx.io.Source`, and plain sequential gzip stays a wrapping
  `RawSource`.

- **Index discovery is outside the core.** The CLI and the desktop app look
  beside the VCF for `.csi`, then `.tbi`, and an explicit `--index`
  overrides that search. CSI is first because it is the format that can
  describe contigs past 512 Mbp. The browser has no sibling-file lookup, so
  the user supplies the index file. Core readers take an already-opened
  index source. They do not search a filesystem.

- **Coordinates are `Long` until the table boundary.** BGZF, CSI, and TBI
  use `Long` throughout. Values are converted to `Position.position: Int`
  at the `GenotypeTable` boundary. A position past `Int.MAX_VALUE` is an
  `IllegalArgumentException`, the same exception `Vcf.read` already uses
  for a position that is not an integer. `Position` itself does not change.
  Wheat fits in `Int`. A record at 3,000,000,000 is a fixture that checks
  only the boundary refuses it.

- **Two oracles, and they do not mix.** Decoded tutorial content stays
  `testdata/golden/` and TASSEL 5 `v5.2.98` (ADR-0009). A bgzipped tutorial
  VCF must produce exactly
  `testdata/golden/vcf/mdp_genotype_w_depth.hmp.txt`. No file under
  `testdata/golden/` changes for this work. Compressed and indexed behavior
  is htslib 1.25 (`bgzip`, `tabix`, `bcftools`) in `testdata/htslib/`, with
  its own `MANIFEST` recording the htslib version and the exact commands.
  Those files follow the same never-edit-by-hand rule and code-owner review
  as `testdata/golden/`. They are not TASSEL 5 goldens, and a change to
  them does not use the `golden-update` label. That label means the TASSEL
  5 oracle changed. A CI check for `testdata/htslib/` waits until the
  directory exists (tassel-6#80).

- **CRC32 and ISIZE are verified** for every gzip member and every BGZF
  block, as htslib does.

- **Parallelism is not a user-facing parameter.** BGZF blocks may decode in
  parallel on the JVM and on Native (`Dispatchers.Default`), with one
  worker on wasmJs. Output is identical for any worker count. There is no
  thread parameter on the operation spec (ADR-0020).

- **Region overlap is `record`.** A record is returned when its reference
  span overlaps the query. The span starts at `POS` and covers the longest
  of the `REF` allele, a valid `INFO/END`, the `INFO/SVLEN` of symbolic
  alleles such as `<DEL>`, and `FORMAT/LEN` for gVCF `<*>` alleles. That is
  `tbx_parse1` in htslib `tbx.c`, kept in sync with `get_rlen` in `vcf.c`.
  `tabix` and `bcftools view -r` behave this way by default, and
  `bcftools view -r` is the oracle. Coordinates in the query are 1-based
  and inclusive, as in VCF and in tabix. This is why the deletion at
  536,870,000 (`ACGT` to `A`) is inside `chr3B:536870002-536870010`.
  `filter-sites` is unchanged: it still keeps a site only when its `POS`
  is inside the requested range. The two operations answer different
  questions.

- **One region per query.** `chr3B`, `chr3B:1000-2000`, and `chr3B:1000-`
  mean the same thing they mean to `tabix` and `bcftools view -r`: the
  whole contig, an inclusive range, and from that position through the end
  of the contig. `chr3B:600000000` means that single position, as in
  `bcftools` (`HTS_PARSE_ONE_COORD`), not tabix's reading of it as "through
  the end of the contig". Commas that sit between digits are thousands
  separators, so `chr3B:599,999,999-600,000,001` is one region. A comma
  that is not between digits is an error, so a list of regions is not
  silently accepted. `{name}:beg-end` is the form for a contig name that
  contains `:`, as in htslib `hts_parse_region`.

- **Contig names match the index, then `Chromosome.of`.** The index keeps
  VCF names as written (`chr3B`). `Chromosome.of` normalizes them (`3B`):
  it trims, upper-cases, and drops a leading `CHROMOSOME` or `CHR`. The
  query is compared to the index's own names exactly, case-sensitive,
  first. If that misses, both sides are compared with `Chromosome.of`. One
  normalized hit is used. Zero hits is an unknown contig. Two or more is
  an error that names the ambiguous index contigs (an index that contains
  both `chr3B` and `3B` is ambiguous for a query of `3B`). `Chromosome.of`
  itself does not change.

- **A typo is an error. An empty range is an empty table.** An unknown
  contig is an error listing the index's contig names as written. An
  inverted region (`beg > end`, such as `chr3B:2000-1000`) is an error.
  htslib 1.25 prints nothing and exits 0 for both; a typo must not look
  like an empty region. A known contig with nothing in range, or a region
  past the index's coordinate limit, returns a table with no sites.

- **A missing BGZF EOF marker is a warning.** htslib warns "EOF marker is
  absent. The input may be truncated" and keeps going. The core does the
  same: it does not throw, and it does not grow `ProgressReporter`. The
  caller shows the warning. The CLI prints it with the existing
  `CliLog.warn`. A truncated block inside the file is an
  `IllegalArgumentException`. How the warning is carried on the new read
  result is tassel-6#84's choice, as long as the read succeeds and the
  user sees that sentence.

- **`Vcf.read` stays the sequential plain-text reader.** It still takes a
  `kotlinx.io.Source`. The indexed entry point is a new suspending
  function that takes a `RandomAccessSource` plus an already-opened index
  and returns an in-memory `GenotypeTable` (ADR-0012). Its name belongs to
  tassel-6#90.

## Consequences

- tassel-6#86, #88, #90, #91, and #93 can add public API against these
  rules. tassel-6#80, #81, #82, and #83 did not wait on this ADR.
- The browser, the CLI, the desktop app, and the Native bindings share one
  decoder. ADR-0008 (BioKotlin for other work) is unchanged. Only the
  indexed-VCF line of the Phase 3 roadmap moves.
- `filter-sites` and `testdata/golden/` stay the TASSEL 5 contracts they
  are. Indexed queries will not match `filter-sites` on a record whose
  span overlaps the region but whose `POS` does not, and that difference
  is this ADR, not a porting bug.
- `commonMain` grows a DEFLATE implementation and the index decoders. The
  cost is code htslib and `java.util.zip` already contain. The benefit is
  one result on every target, with no JVM-only dependency.
- Host core count does not change indexed-read bytes. A missing platform
  inflate is a performance question for a later spike, not a correctness
  question.
- Code owners review `testdata/htslib/` the way they review
  `testdata/golden/`. Agents are told so in `AGENTS.md` before the
  directory exists.

## Alternatives considered

- **htsjdk or BioKotlin on the JVM only.** The Phase 3 roadmap line. It
  leaves the browser and Native without the reader, against ADR-0001 and
  ADR-0015. Rejected.
- **`pos` overlap, matching `filter-sites`.** Users comparing a region
  query with `tabix` and `bcftools view -r` would see a missing deletion.
  `filter-sites` remains available for the `POS`-inside-range question.
  Rejected for the indexed query.
- **tabix's reading of `chr:pos` as "through the end of the contig".**
  `bcftools view -r` is the tool users name for a VCF region, and a single
  coordinate is the narrower, less surprising reading. Rejected.
- **An unknown contig or an inverted region returns an empty table.** That
  is what htslib 1.25 does, and a typo then looks like a region with no
  sites. Rejected.
- **Look for `.tbi` before `.csi`.** TBI cannot describe wheat contigs past
  512 Mbp. CSI is the index that can. Rejected.
- **Fail the file when the BGZF EOF marker is absent.** htslib warns and
  continues. A file that htslib reads would not load. Rejected. A truncated
  block inside the file is still an error.
- **Put the Node file source in `wasmJsMain`.** That source set is also the
  browser compilation of `tassel-core`. Rejected.
