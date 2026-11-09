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
# --only vcf times wheat-scale VCF reads (tassel-6#81, tassel-6#89). It generates a loadable plain
# VCF and a bgzipped file larger than 2 GiB with scripts/bench/generate-wheat-vcf.sh, then times the
# TASSEL 6 CLI on the plain and bgzipped files, TASSEL 5 -vcf on both, BGZF inflate at 1 worker and
# at the host worker count, htslib, and indexed region reads of chr3B (tassel-6#90). vcf needs
# --heap 8g, plus bgzip and bcftools from htslib 1.25.
#
# Requires curl, tar, /usr/bin/time, and a JDK 21+ for TASSEL 5. The TASSEL 6 CLI brings its own
# JDK 25 through the Gradle toolchain. --only vcf also requires htslib 1.25.
#
# Usage: scripts/bench/run-benchmarks.sh [--runs N] [--heap 8g] [--million]
#                                        [--only ANALYSIS] [--dataset NAME]
#
# ANALYSIS includes the pipeline rows (load, summary, filter-sites, filter-taxa, kinship, pca),
# kinship-kernel, and vcf. kinship-kernel times centeredIbs alone and sweeps TASSEL 5 at 1, 4, and
# 8 threads. synthetic-5k (5,000 taxa x 100,000 sites) and, with --million, synthetic-1m
# (1,000 x 1,000,000) are part of that sweep. Both need --heap 16g. A smaller heap skips them,
# unless --dataset names one, in which case the script stops. The 1,000,000-site HapMap is about
# 2 GB of text. vcf needs --heap 8g and does not use --dataset or --million.
set -euo pipefail

TASSEL5_VERSION="5.2.98"
TASSEL5_TAG="v5.2.98"
TASSEL5_URL="https://github.com/maize-genetics/tassel/releases/download/v5.2.98/tassel-5-standalone-v5.2.98.tar.gz"
TASSEL5_SHA256="52e7675aec4c7ab34bb6e0cfe94b26709f8fd8dd761b4d7892b3aa0e38df3c6e"

# Wheat VCF (tassel-6#81). 1,000 x 1,000,000 is the loadable size both tools are asked to read
# at --heap 8g. The region file is streamed and indexed; 5,110,342 sites is the count that
# landed on about 2.5 GiB of bgzip. Calibration: 100,000 sites, htslib 1.25 `bgzip -c -@ 1`,
# 52,527,890 bytes (525.28 bytes/site). 5,110,342 * 525.28 is 2.50 GiB.
VCF_TAXA=1000
VCF_SITES=1000000
VCF_REGION_SITES=5110342
VCF_SEED=81
VCF_BGZIP_THREADS=4
VCF_HTSLIB_VERSION="1.25"

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
    -h|--help) sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
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

# Whole-file VCF reads, BGZF inflate, then htslib on the same loadable file and on a >2 GiB region
# file. TASSEL 6's whole-file row is filter-sites with defaults: read the VCF and write HapMap.
# The .vcf.gz row uses the parallel BGZF reader. Inflate rows are `dev bgzf-bench`: the scan only,
# after one warm-up, JVM startup excluded. htslib inflate is the whole `bgzip -d` process.
run_vcf() {
  require_heap 8192 "vcf"
  command -v bgzip >/dev/null || die "bgzip not found; brew install htslib"
  command -v bcftools >/dev/null || die "bcftools not found; brew install bcftools"
  local bgzip_out bgzip_line bcftools_out bcftools_first bcftools_second
  bgzip_out=$(bgzip --version 2>&1) || die "bgzip --version failed"
  bgzip_line=${bgzip_out%%$'\n'*}
  [[ "$bgzip_line" == "bgzip (htslib) $VCF_HTSLIB_VERSION" ]] ||
    die "bgzip must be htslib $VCF_HTSLIB_VERSION; found '$bgzip_line'"
  bcftools_out=$(bcftools --version 2>&1) || die "bcftools --version failed"
  bcftools_first=${bcftools_out%%$'\n'*}
  bcftools_second=${bcftools_out#*$'\n'}
  bcftools_second=${bcftools_second%%$'\n'*}
  [[ "$bcftools_first" == "bcftools $VCF_HTSLIB_VERSION" && "$bcftools_second" == "Using htslib $VCF_HTSLIB_VERSION" ]] ||
    die "bcftools must be $VCF_HTSLIB_VERSION built with htslib $VCF_HTSLIB_VERSION; found '$bcftools_first' / '$bcftools_second'"

  local gen="$ROOT/scripts/bench/generate-wheat-vcf.sh"
  [[ -x "$gen" ]] || die "no generator at $gen"
  local load="$WORK/wheat-load-${VCF_TAXA}x${VCF_SITES}.vcf"
  local region="$WORK/wheat-region-${VCF_TAXA}x${VCF_REGION_SITES}.vcf.gz"
  "$gen" plain --taxa "$VCF_TAXA" --sites "$VCF_SITES" --seed "$VCF_SEED" --output "$load"
  "$gen" region --taxa "$VCF_TAXA" --sites "$VCF_REGION_SITES" --seed "$VCF_SEED" --output "$region"
  local load_gz="${load}.gz"

  local plain_bytes gz_bytes region_bytes
  plain_bytes=$(file_bytes "$load")
  gz_bytes=$(file_bytes "$load_gz")
  region_bytes=$(file_bytes "$region")

  local t5_java t6_java t6_bin
  t5_java=$("$JAVA" -version 2>&1 || true)
  t6_bin="${JAVA_HOME:+$JAVA_HOME/bin/}java"
  t6_java=$("$t6_bin" -version 2>&1 || true)
  echo
  echo "Host:    $(uname -srm)"
  echo "TASSEL5 JDK: ${t5_java%%$'\n'*}"
  echo "TASSEL6 JDK: ${t6_java%%$'\n'*}"
  echo "TASSEL5: $TASSEL5_VERSION at -maxThreads 4"
  echo "TASSEL6: $(cd "$ROOT" && awk -F= '/^version=/ { print $2 }' gradle.properties)"
  echo "htslib:  $VCF_HTSLIB_VERSION"
  echo "Heap:    -Xmx$HEAP, median of $RUNS runs after one warm-up"
  echo "Loadable: $VCF_TAXA taxa x $VCF_SITES sites, seed $VCF_SEED"
  echo "  plain $(file_bytes_human "$plain_bytes") ($plain_bytes bytes) $load"
  echo "  bgzip $(file_bytes_human "$gz_bytes") ($gz_bytes bytes) $load_gz"
  echo "Region: $VCF_TAXA taxa x $VCF_REGION_SITES sites, seed $VCF_SEED"
  echo "  bgzip $(file_bytes_human "$region_bytes") ($region_bytes bytes) $region"
  echo "  index: bcftools index -c"

  printf '| Tool | Read | Input | Median | Peak | Command |\n'
  printf '|------|------|-------|-------:|-----:|:--------|\n'

  local result command
  command="$TASSEL6 --quiet filter-sites --input $load --output $WORK/t6_vcf_load"
  result=$(repeat_measure "$TASSEL6" --quiet filter-sites --input "$load" --output "$WORK/t6_vcf_load")
  vcf_row "TASSEL 6" "whole-file read" "loadable .vcf" "$result" "$command"

  command="$TASSEL6 --quiet filter-sites --input $load_gz --output $WORK/t6_vcf_gz"
  result=$(repeat_measure "$TASSEL6" --quiet filter-sites --input "$load_gz" --output "$WORK/t6_vcf_gz")
  vcf_row "TASSEL 6" "whole-file read" "loadable .vcf.gz" "$result" "$command"

  bgzf_bench_row "$load_gz" 1
  bgzf_bench_row "$load_gz" ""

  command="$JAVA -Xmx$HEAP -cp $T5/sTASSEL.jar net.maizegenetics.pipeline.TasselPipeline -maxThreads 4 -vcf $load -export $WORK/t5_vcf_load -exportType Hapmap"
  result=$(repeat_measure "$JAVA" "-Xmx$HEAP" -cp "$T5/sTASSEL.jar" net.maizegenetics.pipeline.TasselPipeline \
    -maxThreads 4 -vcf "$load" -export "$WORK/t5_vcf_load" -exportType Hapmap)
  vcf_row "TASSEL 5" "whole-file read" "loadable .vcf" "$result" "$command"

  command="$JAVA -Xmx$HEAP -cp $T5/sTASSEL.jar net.maizegenetics.pipeline.TasselPipeline -maxThreads 4 -vcf $load_gz -export $WORK/t5_vcf_gz -exportType Hapmap"
  result=$(repeat_measure "$JAVA" "-Xmx$HEAP" -cp "$T5/sTASSEL.jar" net.maizegenetics.pipeline.TasselPipeline \
    -maxThreads 4 -vcf "$load_gz" -export "$WORK/t5_vcf_gz" -exportType Hapmap)
  vcf_row "TASSEL 5" "whole-file read" "loadable .vcf.gz" "$result" "$command"

  command="bgzip -d -@ 1 -o /dev/null $load_gz"
  result=$(repeat_measure bgzip -d -@ 1 -o /dev/null "$load_gz")
  vcf_row "htslib" "decompress, 1 thread" "loadable .vcf.gz" "$result" "$command"

  command="bgzip -d -@ $VCF_BGZIP_THREADS -o /dev/null $load_gz"
  result=$(repeat_measure bgzip -d -@ "$VCF_BGZIP_THREADS" -o /dev/null "$load_gz")
  vcf_row "htslib" "decompress, $VCF_BGZIP_THREADS threads" "loadable .vcf.gz" "$result" "$command"

  command="bcftools view -Ou -o /dev/null $load_gz"
  result=$(repeat_measure bcftools view -Ou -o /dev/null "$load_gz")
  vcf_row "htslib" "decompress and parse" "loadable .vcf.gz" "$result" "$command"

  local region_query region_index="${region}.csi"
  [[ -f "$region_index" ]] || die "no CSI beside $region"
  for region_query in "chr3B:1-1000000" "chr3B:1-50000000" "chr3B"; do
    vcf_region_row "$region" "$region_index" "$region_query"
    command="bcftools view -H -r $region_query -o /dev/null $region"
    result=$(repeat_measure bcftools view -H -r "$region_query" -o /dev/null "$region")
    vcf_row "htslib" "region $region_query" "region .vcf.gz" "$result" "$command"
  done

  echo
  echo "TASSEL whole-file rows are the CLI path: read the VCF, write HapMap, JVM startup included."
  echo "TASSEL 6 decompress rows are dev bgzf-bench: inflate only, one warm-up discarded, startup excluded."
  echo "TASSEL 6 region rows are dev vcf-region-bench: one region, one warm-up discarded, startup excluded."
  echo "Full tool output: $LOG"
}

# One dev bgzf-bench invocation. [workers] empty means the host count. Prints a table row.
bgzf_bench_row() {
  local input="$1" workers="$2"
  local -a args=(--quiet dev bgzf-bench --input "$input" --runs "$RUNS")
  if [[ -n "$workers" ]]; then
    args+=(--workers "$workers")
  fi
  local stats="$WORK/.time" out="$WORK/bgzf-bench.out" rss line median count label command
  /usr/bin/time $TIME_FLAG "$TASSEL6" "${args[@]}" >"$out" 2>"$stats" || die "bgzf-bench failed; see $LOG"
  cat "$stats" >>"$LOG"
  if [[ "$TIME_FLAG" == "-l" ]]; then
    rss=$(awk '/maximum resident set size/ { printf "%.0f", $1 / 1048576 }' "$stats")
  else
    rss=$(awk -F': ' '/Maximum resident set size/ { printf "%.0f", $2 / 1024 }' "$stats")
  fi
  line=$(tr -d '\r' <"$out" | awk 'NF { last = $0 } END { print last }')
  median=$(awk '{ for (i = 1; i <= NF; i++) { split($i, p, "="); if (p[1] == "median") print p[2] } }' <<<"$line")
  count=$(awk '{ for (i = 1; i <= NF; i++) { split($i, p, "="); if (p[1] == "workers") print p[2] } }' <<<"$line")
  [[ -n "$median" && -n "$count" ]] || die "bgzf-bench wrote no median: $(cat "$out")"
  label="decompress, $count worker"
  if [[ "$count" != "1" ]]; then
    label="decompress, $count workers"
  fi
  command="$TASSEL6 ${args[*]}"
  printf '| %s | %s | %s | %s s | %s MiB | `%s` |\n' \
    "TASSEL 6" "$label" "loadable .vcf.gz" "$median" "${rss:-0}" "$command"
}

# One dev vcf-region-bench invocation. The median excludes JVM startup. Prints a table row.
vcf_region_row() {
  local input="$1" index="$2" query="$3"
  local -a args=(--quiet dev vcf-region-bench --input "$input" --index "$index" --region "$query" --runs "$RUNS")
  local stats="$WORK/.time" out="$WORK/vcf-region-bench.out" rss line median command
  /usr/bin/time $TIME_FLAG "$TASSEL6" "${args[@]}" >"$out" 2>"$stats" || die "vcf-region-bench failed; see $LOG"
  cat "$stats" >>"$LOG"
  if [[ "$TIME_FLAG" == "-l" ]]; then
    rss=$(awk '/maximum resident set size/ { printf "%.0f", $1 / 1048576 }' "$stats")
  else
    rss=$(awk -F': ' '/Maximum resident set size/ { printf "%.0f", $2 / 1024 }' "$stats")
  fi
  line=$(tr -d '\r' <"$out" | awk 'NF { last = $0 } END { print last }')
  median=$(awk '{ for (i = 1; i <= NF; i++) { split($i, p, "="); if (p[1] == "median") print p[2] } }' <<<"$line")
  [[ -n "$median" ]] || die "vcf-region-bench wrote no median: $(cat "$out")"
  command="$TASSEL6 ${args[*]}"
  printf '| %s | %s | %s | %s s | %s MiB | `%s` |\n' \
    "TASSEL 6" "region $query" "region .vcf.gz" "$median" "${rss:-0}" "$command"
}

file_bytes() {
  if stat -f%z "$1" >/dev/null 2>&1; then
    stat -f%z "$1"
  else
    stat -c%s "$1"
  fi
}

file_bytes_human() {
  awk -v n="$1" 'BEGIN { printf "%.2f GiB", n / 1024 / 1024 / 1024 }'
}

vcf_row() {
  local tool="$1" what="$2" input="$3" result="$4" command="$5"
  printf '| %s | %s | %s | %s s | %s MiB | `%s` |\n' \
    "$tool" "$what" "$input" "${result%% *}" "${result##* }" "$command"
}

if [[ "$ONLY" == "kinship-kernel" ]]; then
  run_kinship_kernel
  exit 0
fi

if [[ "$ONLY" == "vcf" ]]; then
  run_vcf
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
