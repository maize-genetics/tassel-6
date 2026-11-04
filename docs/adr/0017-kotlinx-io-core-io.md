# ADR-0017: kotlinx-io for core I/O

- Status: Accepted (maintainer decision on the `hapmap-io` task, 2026-11-04)
- Date: 2026-11-04

## Context

`tassel-core` has to read and write genotype files on four targets: the JVM, Wasm
in Node and in the browser, macOS arm64, and Linux x64. ADR-0001 keeps
`commonMain` free of JVM-only dependencies, so `java.io` and htsjdk are out, and
Kotlin's standard library has no multiplatform byte stream.

The files are large. The tutorial VCF that the VCF reader (#14) must handle is
18 MB, and real datasets are much larger. ADR-0012 makes streaming the whole
point of the genotype API, so the I/O layer under it has to stream too.

The browser has no file system: the web app's file input (#21) hands over file
contents as bytes already in memory. The same reader must serve that case and a
large file on disk.

## Decision

- `tassel-core` depends on `org.jetbrains.kotlinx:kotlinx-io-core` `0.9.2`,
  declared in `gradle/libs.versions.toml` and exposed as an `api` dependency of
  `commonMain` because it appears in public signatures. The version is published
  for jvm, wasm-js, macosArm64, and linuxX64, and is built against the Kotlin
  2.3.0 standard library, below this project's 2.4.21.
- Readers take a `kotlinx.io.Source` and writers take a `kotlinx.io.Sink`.
- The core reader and writer API stays free of file paths. Callers with a file
  open it themselves (`kotlinx.io.files.SystemFileSystem` covers the JVM, Native,
  and Node); callers with bytes in memory, such as the browser, wrap them in a
  `kotlinx.io.Buffer`.

## Consequences

- One reader implementation serves files, in-memory bytes, and later network or
  compressed streams, because they are all `RawSource`s.
- Tests feed readers a `Buffer` over bytes from the `TestData` helper, which
  keeps `commonTest` working on every target without a file-system `expect`.
- Compression stays a wrapping `RawSource`, which is how gzip can be added for
  the JVM first without changing the readers.
- A new runtime dependency ships in every artifact, including the Wasm bundle.
  It is small and has no transitive dependencies beyond the Kotlin standard
  library.
- kotlinx-io is pre-1.0, so its API can still change. Keeping `Source` and `Sink`
  at the edge of `tassel-core` (parsers work on lines and bytes) limits what a
  breaking release would touch.

## Alternatives considered

- Readers taking a whole `ByteArray`: no new dependency, but it cannot stream the
  18 MB VCF that #14 needs, and on Wasm the bytes would also have to cross the JS
  boundary as a 36 MB string.
- A hand-rolled line-source interface in `commonMain`: no new dependency, but
  every platform still needs its own file plumbing behind it, which is the part
  kotlinx-io already provides and tests.
- Okio: the same design and a more settled API, but kotlinx-io is the Kotlin
  team's successor to it and matches the kotlinx libraries `commonMain` already
  uses.
