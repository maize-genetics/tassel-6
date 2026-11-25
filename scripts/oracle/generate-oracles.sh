#!/usr/bin/env bash
# Generate the independent oracle fixtures in testdata/oracle/ (ADR-0029 kind 2, tassel-6#183).
#
# Rebuilds testdata/oracle/ and its MANIFEST from scratch with pinned tools run in containers:
#   * NIST StRD linear-regression datasets (public domain), normalized by nist.py;
#   * mpmath log-survival grids for F, t, and chi-square, from tails.py;
#   * R lm()/anova()/drop1() fixtures, from lm-fixtures.R (seed set; tassel-6#186 extends it).
# Output is deterministic: running the script twice gives byte-identical files. Nothing outside
# testdata/oracle/ is modified. Run by people regenerating fixtures, never by CI or Gradle.
#
# Requires docker (a running daemon), curl, and a shasum or sha256sum. The first run pulls the two
# images by digest and downloads the NIST files into build/oracle/ (git-ignored); every later run
# verifies the NIST checksums. pip installs the pinned mpmath inside the container on each run.
# Containers run as linux/arm64 so R's floating-point output is the same on every host.
#
# Usage: scripts/oracle/generate-oracles.sh
set -euo pipefail

PLATFORM="linux/arm64"
PYTHON_IMAGE="python@sha256:2199a62885a12290dc9c5be3ca0681d367576ab7bf037da120e564723292a2f0"
PYTHON_TAG="python:3.12.8-slim"
R_IMAGE="rocker/r-ver@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c"
R_TAG="rocker/r-ver:4.4.2"
MPMATH_VERSION="1.3.0"

NIST_BASE="https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA"
# name sha256 of the NIST .dat file, as downloaded.
NIST_FILES=(
  "Norris 5ab6906c68a9c4b2eb1661c0cd7d2a09f5a469a6d4e05b6b691f190f5f7b944f"
  "Pontius de60baa5dc66fd804f9316f895bc5aa2b8f791c09d1ecc1ec39f8b199b8c15c7"
  "NoInt1 d8234c428dedaa9de23b67df7d04aeee83befae73ceb6b38321b224955a9b667"
  "NoInt2 a71c2810b1edcef8fec6b864c20c5f9d9908e530282a174772ce610999ca7c69"
  "Filip 403b34689e401d915cb28bd69311d3ca5c3007d88d3b670b00e3203042315072"
  "Longley fc4b0c824f8f0b834eb8e7d1c4c7e0b3cebdd7987bcb085e8f5e1007582480ea"
  "Wampler1 8bcd6f00dbe8c307ebc5c971e12b950eb3cdb75afbc2b38573082b5fcb0d7fbc"
  "Wampler2 01ba5f08b5b9f84431f3754216e006939a80264a1caa69fe08b3cbc543f9ee3e"
  "Wampler3 c5b025f6ebf0365de92d80e478ad826b01a6dfdafb71964ea8bae7eec993ce53"
  "Wampler4 dcff32075a52ecf31af70dad340776237544b1a29b5034a40f8b3e6b607cf9c1"
  "Wampler5 c4e260a6638db61287bfea75a2dd3b324c86c5d40931ddc6d8fab0c0305ffa39"
)

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SCRIPTS="$ROOT/scripts/oracle"
CACHE="$ROOT/build/oracle"
STAGE="$CACHE/stage"
ORACLE="$ROOT/testdata/oracle"
MANIFEST="$STAGE/MANIFEST"

die() {
  echo "generate-oracles: $*" >&2
  exit 1
}

sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

manifest() { printf '%s\n' "$@" >>"$MANIFEST"; }

docker info >/dev/null 2>&1 || die "the docker daemon is not running"
command -v curl >/dev/null || die "curl is required"

# --- Inputs -------------------------------------------------------------------

mkdir -p "$CACHE/nist"
for entry in "${NIST_FILES[@]}"; do
  read -r name hash <<<"$entry"
  file="$CACHE/nist/$name.dat"
  if [[ ! -f "$file" ]]; then
    echo "Downloading $NIST_BASE/$name.dat"
    curl -fsSL -o "$file.part" "$NIST_BASE/$name.dat"
    mv "$file.part" "$file"
  fi
  [[ "$(sha256 "$file")" == "$hash" ]] || die "checksum mismatch for $file; delete it and rerun"
done

docker pull -q --platform "$PLATFORM" "$PYTHON_IMAGE" >/dev/null
docker pull -q --platform "$PLATFORM" "$R_IMAGE" >/dev/null

rm -rf "$STAGE"
mkdir -p "$STAGE/nist" "$STAGE/tails" "$STAGE/r"

run_python() {
  docker run --rm --platform "$PLATFORM" --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e PYTHONUSERBASE=/tmp/py -e PYTHONHASHSEED=0 \
    -v "$SCRIPTS:/work/scripts:ro" -v "$CACHE/nist:/work/in:ro" -v "$STAGE:/work/out" \
    "$PYTHON_IMAGE" sh -c "pip install -q --disable-pip-version-check --no-cache-dir --user mpmath==$MPMATH_VERSION && $1"
}

run_r() {
  docker run --rm --platform "$PLATFORM" --user "$(id -u):$(id -g)" \
    -v "$SCRIPTS:/work/scripts:ro" -v "$STAGE:/work/out" \
    "$R_IMAGE" "$@"
}

# --- Fixtures -----------------------------------------------------------------

run_python "python /work/scripts/nist.py /work/in /work/out/nist"
run_python "python /work/scripts/tails.py /work/out/tails/log-sf.tsv"
run_r Rscript /work/scripts/lm-fixtures.R /work/out/r

PYTHON_VERSION="$(run_python "python --version" | tr -d '\r')"
MPMATH_INSTALLED="$(run_python "python -c 'import mpmath; print(mpmath.__version__)'" | tr -d '\r')"
[[ "$MPMATH_INSTALLED" == "$MPMATH_VERSION" ]] || die "mpmath $MPMATH_INSTALLED installed, expected $MPMATH_VERSION"
R_VERSION="$(run_r R --version | head -n 1)"

# --- MANIFEST -----------------------------------------------------------------

# Starting tolerances for a port to meet, in the syntax of testdata/golden/MANIFEST (ADR-0029).
# The issue that consumes a fixture may tighten its line and records why in its PR.
NIST_EASY_TOLERANCE="rel 1e-10"
NIST_HARD_TOLERANCE="rel 1e-8"   # Wampler4 and Wampler5: large noise, condition number 1e12 or so
NIST_FILIP_TOLERANCE="rel 1e-6"  # a degree-10 polynomial; double-precision QR reaches about 7 digits
TAIL_TOLERANCE="rel 1e-10"       # on log_sf, natural log
R_TOLERANCE="rel 1e-10"

: >"$MANIFEST"
manifest \
  "# Independent oracle fixtures (ADR-0029 kind 2, tassel-6#183). Written by" \
  "# scripts/oracle/generate-oracles.sh from pinned containers, never edited by hand, never" \
  "# regenerated by CI. Not TASSEL 5 goldens: no golden-update label; owner review via CODEOWNERS." \
  "#" \
  "# Format: [section] headers, then 'key = value' lines. A key that starts with oracle/ is a" \
  "# fixture file (relative to testdata/) and its value is its tolerance, in the syntax of" \
  "# testdata/golden/MANIFEST: exact | rel X | abs X, optionally followed by up-to-sign or" \
  "# rows-unordered. Tolerances are starting points for the port to meet, set per file here," \
  "# never in test code. 'command' may repeat. Tool versions and image digests are recorded." \
  "" \
  "[tools]" \
  "platform = $PLATFORM" \
  "python_image = $PYTHON_IMAGE" \
  "python_tag = $PYTHON_TAG" \
  "python_version = $PYTHON_VERSION" \
  "mpmath_version = $MPMATH_INSTALLED" \
  "r_image = $R_IMAGE" \
  "r_tag = $R_TAG" \
  "r_version = $R_VERSION" \
  "" \
  "[nist]" \
  "source = https://www.itl.nist.gov/div898/strd/lls/lls.shtml" \
  "citation = NIST Statistical Reference Datasets (StRD), Information Technology Laboratory, linear regression datasets" \
  "license = public domain (work of the US government)" \
  "base_url = $NIST_BASE" \
  "format = <Name>.data.tsv holds the data as NIST prints it; <Name>.certified.txt holds 'key = value' certified values copied as printed; B<i> = estimate standard_error" \
  "command = python /work/scripts/nist.py /work/in /work/out/nist"
for entry in "${NIST_FILES[@]}"; do
  read -r name hash <<<"$entry"
  manifest "source_sha256.$name = $hash"
done
for entry in "${NIST_FILES[@]}"; do
  read -r name hash <<<"$entry"
  case "$name" in
    Wampler4 | Wampler5) tolerance="$NIST_HARD_TOLERANCE" ;;
    Filip) tolerance="$NIST_FILIP_TOLERANCE" ;;
    *) tolerance="$NIST_EASY_TOLERANCE" ;;
  esac
  manifest "oracle/nist/$name.certified.txt = $tolerance" "oracle/nist/$name.data.tsv = exact"
done

manifest \
  "" \
  "[tails]" \
  "description = natural log of the survival function P(X > x) for F, Student's t, and chi-square" \
  "precision = mpmath at 80 digits, printed to 50 significant digits" \
  "columns = dist df1 df2 x log_sf (df2 is - except for f; x is an exact double printed with repr)" \
  "range = log p from -1 down to -50000, plus small-x and small-|t|, large-df rows that scripts/golden/Distributions.java drops as ill-conditioned; rows with x above 1e300 are not representable and are omitted" \
  "command = python /work/scripts/tails.py /work/out/tails/log-sf.tsv" \
  "oracle/tails/log-sf.tsv = $TAIL_TOLERANCE" \
  "" \
  "[r-lm]" \
  "description = R lm(), anova(), drop1(F) on a small covariate-plus-factor design" \
  "status = seed fixture only; tassel-6#186 adds the designs, contrasts, collinearity, and missing-row cases" \
  "command = Rscript /work/scripts/lm-fixtures.R /work/out/r" \
  "oracle/r/lm-smoke.data.tsv = exact"
for part in coefficients anova drop1; do
  manifest "oracle/r/lm-smoke.$part.tsv = $R_TOLERANCE"
done

# --- Publish ------------------------------------------------------------------

rm -rf "$ORACLE"
mkdir -p "$ORACLE"
cp -R "$STAGE/nist" "$STAGE/tails" "$STAGE/r" "$ORACLE/"
cp "$MANIFEST" "$ORACLE/MANIFEST"
echo "Wrote $(find "$ORACLE" -type f | wc -l | tr -d ' ') files under testdata/oracle/"
