# Architecture Decision Records

Each ADR records one decision: the context, what was decided, and what follows
from it. ADRs are immutable once accepted. To change a decision, write a new ADR
that supersedes the old one and update the old one's status line only.

Copy [`0000-template.md`](0000-template.md) to start a new one. Number
sequentially; use a short kebab-case title.

| ADR | Title | Status |
|-----|-------|--------|
| [0001](0001-kotlin-multiplatform-pure-core.md) | Kotlin Multiplatform with a pure-Kotlin common core | Accepted |
| [0002](0002-operation-specs-single-source.md) | Operation specs as the single source of truth | Accepted |
| [0003](0003-compose-multiplatform-ui.md) | Compose Multiplatform for desktop and web | Accepted |
| [0004](0004-hybrid-web-compute.md) | Hybrid web compute: in-browser Wasm first, JVM server later | Accepted |
| [0005](0005-legacy-native-ports-plus-bridge.md) | Legacy analyses: native ports plus a TASSEL 5 bridge | Accepted |
| [0006](0006-modern-cli-with-legacy-translator.md) | Modern CLI with a TASSEL 5 pipeline translator | Accepted |
| [0007](0007-native-bindings-r-python.md) | R and Python bindings through a Kotlin/Native C API | Accepted |
| [0008](0008-biokotlin-core.md) | Reuse BioKotlin through a multiplatform `biokotlin-core` | Accepted |
| [0009](0009-strict-golden-parity.md) | Strict numerical parity through TASSEL 5 golden tests | Accepted |
| [0010](0010-git-flow-and-versioning.md) | Git flow branching, Conventional Commits, SemVer | Accepted |
| [0011](0011-jdk-25.md) | JDK 25 for JVM targets | Accepted |
| [0012](0012-streaming-genotype-access.md) | Streaming site-block genotype access; out-of-core storage later | Accepted |
| [0013](0013-teaching-first-release.md) | Teaching-first public release | Accepted |
| [0014](0014-agentic-mcp-support.md) | Agent support: local MCP server, in-app assistant, CLI skill | Accepted |
| [0015](0015-pangenome-complements.md) | Pangenome complements: consume grits/ropebwt3-phg outputs and add a viewer | Accepted |
| [0016](0016-agentic-development-workflow.md) | Agent-driven development workflow | Accepted |
