#!/usr/bin/env bash
# Time TASSEL 5 against the TASSEL 6 JVM build on the same inputs and the same analyses
# (tassel-6#24), and print the Markdown table docs/poc/benchmarks.md carries.
#
# Downloads the pinned TASSEL 5 standalone release into build/bench/ (git-ignored), builds the
# TASSEL 6 CLI, generates the synthetic inputs with `tassel dev simulate`, then runs each analysis
# on both tools under /usr/bin/time, reporting the median wall time and the peak resident set over
# --runs runs after one warm-up.
#
# Both tools get the same -Xmx so neither wins on heap alone, and the pipeline rows pin TASSEL 5
# at -maxThreads 4 for the same reason scripts/generate-goldens.sh does: it fixes Centered IBS at
# 1,000-site blocks. kinship-kernel sweeps 1, 4, and 8 threads instead, and times TASSEL 6's
# centeredIbs directly.
#
# Requires curl, tar, /usr/bin/time, and a JDK 21+ for TASSEL 5. The TASSEL 6 CLI brings its own
# JDK 25 through the Gradle toolchain.
#
# Usage: scripts/bench/run-benchmarks.sh [--runs N] [--heap 8g] [--million]
#                                        [--only ANALYSIS] [--dataset NAME]
#
# ANALYSIS includes the pipeline rows (load, summary, filter-sites, filter-taxa, kinship, pca)
# and kinship-kernel, which times centeredIbs alone and sweeps TASSEL 5 at 1, 4, and 8 threads.
# synthetic-5k (5,000 taxa x 100,000 sites) and, with --million, synthetic-1m (1,000 x 1,000,000)
# are part of that sweep. Both need --heap 16g. A smaller heap skips them, unless --dataset names
# one, in which case the script stops. The 1,000,000-site HapMap is about 2 GB of text.
set -euo pipefail

TASSEL5_VERSION="5.2.98"
TASSEL5_TAG="v5.2.98"
TASSEL5_URL="https://github.com/maize-genetics/tassel/releases/download/v5.2.98/tassel-5-standalone-v5.2.98.tar.gz"
TASSEL5_SHA256="52e7675aec4c7ab34bb6e0cfe94b26709f8fd8dd761b4d7892b3aa0e38df3c6e"

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CACHE="$ROOT/build/bench"
WORK="$CACHE/work"
LOG="$CACHE/bench.log"

RUNS=3
HEAP="8g"
MILLION=0
ONLY=""
DATASET=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS="$2"; shift 2 ;;
    --heap) HEAP="$2"; shift 2 ;;
    --million) MILLION=1; shift ;;
    --only) ONLY="$2"; shift 2 ;;
    --dataset) DATASET="$2"; shift 2 ;;
    -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

die() {
  echo "run-benchmarks: $*" >&2
  exit 1
}

sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

# --- Toolchain ----------------------------------------------------------------

JAVA="${JAVA_HOME:+$JAVA_HOME/bin/}java"
command -v "$JAVA" >/dev/null || die "java not found; set JAVA_HOME or put a JDK 21+ on PATH"
JAVA_MAJOR="$("$JAVA" -XshowSettings:properties -version 2>&1 |
  awk -F' = ' '/java\.specification\.version/ { print $2 }')"
[[ "$JAVA_MAJOR" =~ ^[0-9]+$ ]] && ((JAVA_MAJOR >= 21)) ||
  die "TASSEL $TASSEL5_VERSION needs Java 21 or newer; found '$JAVA_MAJOR'"
# TASSEL 6 class files are JDK 25 (ADR-0011). The launcher uses JAVA_HOME and otherwise
# `java` on PATH. TASSEL 5 keeps the JDK chosen above, which only has to be 21 or newer.
if ((JAVA_MAJOR < 25)); then
  JDK25_HOME=""
  for candidate in "$HOME"/.gradle/jdks/*/jdk-*/Contents/Home; do
    [[ -x "$candidate/bin/java" ]] || continue
    major=$("$candidate/bin/java" -XshowSettings:properties -version 2>&1 |
      awk -F' = ' '/java\.specification\.version/ { print $2; exit }')
    if [[ "$major" == "25" ]]; then
      JDK25_HOME="$candidate"
      break
    fi
  done
  [[ -n "$JDK25_HOME" ]] ||
    die "TASSEL 6 needs JDK 25; found Java $JAVA_MAJOR and no toolchain JDK under ~/.gradle/jdks"
  export JAVA_HOME="$JDK25_HOME"
fi
command -v /usr/bin/time >/dev/null || die "/usr/bin/time not found; it is what reports peak memory"

mkdir -p "$CACHE" "$WORK"
: >"$LOG"

if [[ ! -f "$CACHE/tassel5.tar.gz" ]]; then
  echo "Downloading TASSEL $TASSEL5_VERSION"
  curl -fsSL -o "$CACHE/tassel5.tar.gz.part" "$TASSEL5_URL"
  mv "$CACHE/tassel5.tar.gz.part" "$CACHE/tassel5.tar.gz"
fi
[[ "$(sha256 "$CACHE/tassel5.tar.gz")" == "$TASSEL5_SHA256" ]] ||
  die "checksum mismatch for $CACHE/tassel5.tar.gz; delete it and rerun"

T5="$CACHE/tassel-5-standalone-$TASSEL5_TAG"
[[ -f "$T5/sTASSEL.jar" ]] || tar -xzf "$CACHE/tassel5.tar.gz" -C "$CACHE"

echo "Building the TASSEL 6 CLI"
(cd "$ROOT" && ./gradlew --quiet :tassel-cli:installDist) || die "could not build tassel-cli"
TASSEL6="$ROOT/tassel-cli/build/install/tassel/bin/tassel"
[[ -x "$TASSEL6" ]] || die "no CLI at $TASSEL6"

# --- Measurement ---------------------------------------------------------------
#
# /usr/bin/time spells its output differently on each platform, so the parsing is split and the
# units are normalized here: seconds as a decimal, peak resident set in mebibytes.

case "$(uname -s)" in
  Darwin) TIME_FLAG="-l" ;;
  Linux) TIME_FLAG="-v" ;;
  *) die "unsupported platform $(uname -s); /usr/bin/time differs" ;;
esac

# Set while timing a command whose exit code is not a result: TasselPipeline with no plugin prints
# its usage and exits 1, which is exactly the "started and did nothing" the startup row wants.
ALLOW_FAILURE=0

# Wall seconds and peak RSS in MiB for one run, as "SECONDS MIB". Output goes to the log.
measure() {
  local stats="$WORK/.time"
  local start end
  start=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
  /usr/bin/time $TIME_FLAG "$@" >>"$LOG" 2>"$stats" || ((ALLOW_FAILURE)) || {
    cat "$stats" >>"$LOG"
    die "command failed: $*; see $LOG"
  }
  end=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
  local seconds rss_mib
  seconds=$(awk -v a="$start" -v b="$end" 'BEGIN { printf "%.2f", b - a }')
  if [[ "$TIME_FLAG" == "-l" ]]; then
    # macOS reports "<bytes>  maximum resident set size".
    rss_mib=$(awk '/maximum resident set size/ { printf "%.0f", $1 / 1048576 }' "$stats")
  else
    # GNU reports "Maximum resident set size (kbytes): <kib>".
    rss_mib=$(awk -F': ' '/Maximum resident set size/ { printf "%.0f", $2 / 1024 }' "$stats")
  fi
  cat "$stats" >>"$LOG"
  echo "$seconds ${rss_mib:-0}"
}

# Median of RUNS measurements of a command, after one warm-up, as "SECONDS MIB".
repeat_measure() {
  measure "$@" >/dev/null
  local times=() rss=() line
  local i
  for ((i = 0; i < RUNS; i++)); do
    line=$(measure "$@")
    times+=("${line%% *}")
    rss+=("${line##* }")
  done
  local median_time median_rss
  median_time=$(printf '%s\n' "${times[@]}" | sort -g | awk -v n="$RUNS" 'NR == int((n + 1) / 2)')
  median_rss=$(printf '%s\n' "${rss[@]}" | sort -g | awk -v n="$RUNS" 'NR == int((n + 1) / 2)')
  echo "$median_time $median_rss"
}

TASSEL5_CMD=("$JAVA" "-Xmx$HEAP" -cp "$T5/sTASSEL.jar" net.maizegenetics.pipeline.TasselPipeline -maxThreads 4)
export JAVA_OPTS="-Xmx$HEAP --enable-native-access=ALL-UNNAMED"

# --- Inputs ---------------------------------------------------------------------

TUTORIAL="$ROOT/testdata/tutorial/mdp_genotype.hmp.txt"
[[ -f "$TUTORIAL" ]] || die "no tutorial input at $TUTORIAL; run scripts/generate-goldens.sh first"

simulate() {
  local taxa="$1" sites="$2" path="$3"
  [[ -f "$path" ]] && return 0
  echo "Simulating $taxa taxa x $sites sites"
  "$TASSEL6" --quiet dev simulate --taxa "$taxa" --sites "$sites" --seed 0 --output "$path" >>"$LOG" 2>&1 ||
    die "could not simulate $taxa x $sites; see $LOG"
}

heap_mib() {
  local spec="$1"
  local number="${spec%[GgMm]}"
  case "$spec" in
    *g|*G) echo $((number * 1024)) ;;
    *m|*M) echo "$number" ;;
    *) die "unrecognized heap '$spec'; use 8g or 512m" ;;
  esac
}

require_heap() {
  local need="$1" why="$2" have
  have=$(heap_mib "$HEAP")
  ((have >= need)) ||
    die "$why needs --heap $((need / 1024))g or more (currently --heap $HEAP)"
}

# Times centeredIbs, then TASSEL 5's kinship pipeline at 1, 4, and 8 threads.
# TASSEL 5's numbers include JVM startup and the file read; TASSEL 6's do not.
run_kinship_kernel() {
  local sets=() name
  kernel_wanted() { [[ -z "$DATASET" || "$DATASET" == "$1" ]]; }
  kernel_wanted tutorial && sets+=("tutorial|$TUTORIAL|1024")
  if kernel_wanted synthetic-100k; then
    local path="$WORK/sim_1000x100000.hmp.txt"
    simulate 1000 100000 "$path"
    sets+=("synthetic-100k|$path|1024")
  fi
  if kernel_wanted synthetic-5k; then
    if (( $(heap_mib "$HEAP") < 16384 )); then
      [[ -z "$DATASET" ]] || die "synthetic-5k needs --heap 16g (5,000 taxa x 100,000 sites; currently --heap $HEAP)"
      echo "skipping synthetic-5k: needs --heap 16g (currently --heap $HEAP)"
    else
      local path="$WORK/sim_5000x100000.hmp.txt"
      simulate 5000 100000 "$path"
      sets+=("synthetic-5k|$path|16384")
    fi
  fi
  if ((MILLION)) && kernel_wanted synthetic-1m; then
    if (( $(heap_mib "$HEAP") < 16384 )); then
      [[ -z "$DATASET" ]] ||
        die "synthetic-1m needs --heap 16g (1,000 taxa x 1,000,000 sites, about 2 GB of HapMap text; currently --heap $HEAP)"
      echo "skipping synthetic-1m: needs --heap 16g (currently --heap $HEAP)"
    else
      local path="$WORK/sim_1000x1000000.hmp.txt"
      simulate 1000 1000000 "$path"
      sets+=("synthetic-1m|$path|16384")
    fi
  fi
  ((${#sets[@]})) || die "no kinship-kernel datasets selected; --dataset was '$DATASET'"

  echo
  echo "Host:    $(uname -srm)"
  echo "Java:    $("$JAVA" -version 2>&1 | head -1)"
  echo "TASSEL5: $TASSEL5_VERSION at -maxThreads 1, 4, and 8 (pipeline, startup included)"
  echo "TASSEL6: kinship-bench median of $RUNS after one warm-up, file read not included"
  echo "Heap:    -Xmx$HEAP"
  echo
  printf '| Dataset | T6 kernel median | T6 peak | T5 1 thread | T5 4 threads | T5 8 threads | Kernel line |\n'
  printf '|---------|-----------------:|--------:|------------:|-------------:|-------------:|:------------|\n'

  local entry dataset input need stats out rss median threads five line
  for entry in "${sets[@]}"; do
    dataset="${entry%%|*}"
    input="${entry#*|}"; input="${input%|*}"
    need="${entry##*|}"
    require_heap "$need" "$dataset"
    stats="$WORK/.time"
    out="$WORK/kernel.out"
    /usr/bin/time $TIME_FLAG "$TASSEL6" --quiet dev kinship-bench --input "$input" --runs "$RUNS" \
      >"$out" 2>"$stats" || die "kinship-bench failed on $dataset; see $LOG"
    cat "$stats" >>"$LOG"
    if [[ "$TIME_FLAG" == "-l" ]]; then
      rss=$(awk '/maximum resident set size/ { printf "%.0f", $1 / 1048576 }' "$stats")
    else
      rss=$(awk -F': ' '/Maximum resident set size/ { printf "%.0f", $2 / 1024 }' "$stats")
    fi
    line=$(tr -d '\r' <"$out" | awk 'NF { last = $0 } END { print last }')
    median=$(awk '{ for (i = 1; i <= NF; i++) { split($i, p, "="); if (p[1] == "median") print p[2] } }' <<<"$line")
    [[ -n "$median" ]] || die "kinship-bench wrote no median for $dataset: $(cat "$out")"
    local -a thread_times=()
    for threads in 1 4 8; do
      CMD=("$JAVA" "-Xmx$HEAP" -cp "$T5/sTASSEL.jar" net.maizegenetics.pipeline.TasselPipeline \
        -maxThreads "$threads" -h "$input" -KinshipPlugin -method Centered_IBS -endPlugin \
        -export "$WORK/t5_${dataset}_k$threads" -exportType SqrMatrix)
      five=$(repeat_measure "${CMD[@]}")
      thread_times+=("${five%% *}")
    done
    printf '| %s | %s s | %s MiB | %s s | %s s | %s s | `%s` |\n' \
      "$dataset" "$median" "${rss:-0}" "${thread_times[0]}" "${thread_times[1]}" "${thread_times[2]}" "$line"
  done
  echo
  echo "TASSEL 5 columns are the full pipeline (startup included). Full tool output: $LOG"
}

if [[ "$ONLY" == "kinship-kernel" ]]; then
  run_kinship_kernel
  exit 0
fi

DATASETS=()
wanted() { [[ -z "$DATASET" || "$DATASET" == "$1" ]]; }

wanted tutorial && DATASETS+=("tutorial:$TUTORIAL")

if wanted synthetic-100k; then
  SYNTHETIC_100K="$WORK/sim_1000x100000.hmp.txt"
  simulate 1000 100000 "$SYNTHETIC_100K"
  DATASETS+=("synthetic-100k:$SYNTHETIC_100K")
fi

# The 1M-site set is the roadmap's second cut: a 1000 x 1,000,000 HapMap is about 2 GB of text and
# 1 GB of genotypes in memory, so it needs a large heap and is opt-in rather than part of a routine
# run.
if ((MILLION)) && wanted synthetic-1m; then
  SYNTHETIC_1M="$WORK/sim_1000x1000000.hmp.txt"
  simulate 1000 1000000 "$SYNTHETIC_1M"
  DATASETS+=("synthetic-1m:$SYNTHETIC_1M")
fi

((${#DATASETS[@]})) || die "no datasets selected; --dataset was '$DATASET'"

# --- Analyses --------------------------------------------------------------------
#
# Each analysis names the TASSEL 5 pipeline and the TASSEL 6 subcommand that do the same work.
# "load" is read-then-write, because neither tool exposes a read on its own and a read whose result
# is thrown away is not a thing either tool can be asked for.

ANALYSES=(load summary filter-sites filter-taxa kinship pca)

# Fills CMD with the TASSEL 5 invocation for one analysis.
tassel5_command() {
  local analysis="$1" input="$2" out="$3"
  if [[ "$analysis" == "startup" ]]; then
    CMD=("${TASSEL5_CMD[@]}")
    return
  fi
  CMD=("${TASSEL5_CMD[@]}" -h "$input")
  case "$analysis" in
    load) CMD+=(-export "$out" -exportType Hapmap) ;;
    summary) CMD+=(-GenotypeSummaryPlugin -endPlugin
      -export "$out.overview,$out.alleles,$out.sites,$out.taxa" -exportType Table) ;;
    filter-sites) CMD+=(-FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin
      -export "$out,$out.filter" -exportType Hapmap) ;;
    filter-taxa) CMD+=(-FilterTaxaBuilderPlugin -minNotMissing 0.9 -endPlugin
      -export "$out,$out.filter" -exportType Hapmap) ;;
    kinship) CMD+=(-KinshipPlugin -method Centered_IBS -endPlugin -export "$out" -exportType SqrMatrix) ;;
    pca) CMD+=(-PrincipalComponentsPlugin -covariance true -endPlugin
      -export "$out.pcs,$out.eigenvalues,$out.eigenvectors" -exportType Table) ;;
  esac
}

# Fills CMD with the TASSEL 6 invocation for one analysis.
tassel6_command() {
  local analysis="$1" input="$2" out="$3"
  CMD=("$TASSEL6" --quiet)
  case "$analysis" in
    startup) CMD+=(--version) ;;
    # No TASSEL 6 subcommand only reads, so "load" is a filter that keeps everything, which is
    # read plus write and so comparable with TASSEL 5's bare export.
    load) CMD+=(filter-sites --input "$input" --output "$out") ;;
    summary) CMD+=(summary --input "$input" --output "$out") ;;
    filter-sites) CMD+=(filter-sites --input "$input" --min-maf 0.05 --output "$out") ;;
    filter-taxa) CMD+=(filter-taxa --input "$input" --min-not-missing 0.9 --output "$out") ;;
    kinship) CMD+=(kinship --input "$input" --output "$out") ;;
    pca) CMD+=(pca --input "$input" --output "$out") ;;
  esac
}

# --- Run ---------------------------------------------------------------------------

echo
echo "Host:    $(uname -srm)"
echo "Java:    $("$JAVA" -version 2>&1 | head -1)"
echo "TASSEL5: $TASSEL5_VERSION"
echo "TASSEL6: $(cd "$ROOT" && awk -F= '/^version=/ { print $2 }' gradle.properties)"
echo "Heap:    -Xmx$HEAP, median of $RUNS runs after one warm-up"
echo

printf '| Dataset | Analysis | TASSEL 5 time | TASSEL 6 time | Ratio | TASSEL 5 peak | TASSEL 6 peak |\n'
printf '|---------|----------|--------------:|--------------:|------:|--------------:|--------------:|\n'

# One row per analysis, plus a first row that does no work at all. The startup row is here because
# TASSEL 5 spends several seconds starting its pipeline before it reads anything, which on a
# tutorial-sized input is most of the measurement; without it the table reads as a statement about
# the analyses, and it is not.
row() {
  local dataset="$1" analysis="$2" five="$3" six="$4"
  local ratio
  ratio=$(awk -v a="${five%% *}" -v b="${six%% *}" 'BEGIN { if (b > 0) printf "%.2f", a / b; else print "-" }')
  printf '| %s | %s | %s s | %s s | %sx | %s MiB | %s MiB |\n' \
    "$dataset" "$analysis" "${five%% *}" "${six%% *}" "$ratio" "${five##* }" "${six##* }"
}

if [[ -z "$ONLY" || "$ONLY" == "startup" ]]; then
  ALLOW_FAILURE=1
  tassel5_command startup "" ""
  startup5=$(repeat_measure "${CMD[@]}")
  tassel6_command startup "" ""
  startup6=$(repeat_measure "${CMD[@]}")
  ALLOW_FAILURE=0
  row "none" "startup only" "$startup5" "$startup6"
fi

for entry in "${DATASETS[@]}"; do
  dataset="${entry%%:*}"
  input="${entry#*:}"
  for analysis in "${ANALYSES[@]}"; do
    [[ -n "$ONLY" && "$ONLY" != "$analysis" ]] && continue
    tassel5_command "$analysis" "$input" "$WORK/t5_${dataset}_$analysis"
    five=$(repeat_measure "${CMD[@]}")
    tassel6_command "$analysis" "$input" "$WORK/t6_${dataset}_$analysis"
    six=$(repeat_measure "${CMD[@]}")
    row "$dataset" "$analysis" "$five" "$six"
  done
done

echo
echo "Ratio above 1 means TASSEL 6 is faster. Full tool output: $LOG"
