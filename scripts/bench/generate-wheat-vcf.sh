#!/usr/bin/env bash
# Deterministic wheat-scale VCF for scripts/bench/run-benchmarks.sh (tassel-6#81).
#
# Contig lengths are IWGSC RefSeq v2.1 Table 4 Length (bp), Zhu et al., The Plant
# Genome 14, e20112 (2021), the same 21 lengths as scripts/generate-htslib-fixtures.sh.
# Every chromosome, including chr6D, gets at least one site. Extra sites are split
# in proportion to length (largest remainder, ties keep chromosome order). Positions
# on a chromosome run from 1 through its length. Genotypes are biallelic A/C calls
# from a Park-Miller LCG, in chromosome order and then taxon order, so the same
# taxa, sites, and seed always write the same bytes.
#
# plain writes PATH.vcf and a single-threaded bgzip PATH.vcf.gz.
# region streams into bgzip (the plain text is never stored), writes a CSI with
# `bcftools index -c`, and refuses a file that is not larger than 2 GiB.
#
# A stamp next to the output records taxa, sites, and seed. A matching stamp
# reuses the files. Writes go to a .part name and are renamed after success.
#
# Requires bgzip and bcftools from htslib 1.25. No other tools.
#
# Usage: scripts/bench/generate-wheat-vcf.sh plain|region \
#          --taxa N --sites N --seed N --output PATH
set -euo pipefail

HTSLIB_VERSION="1.25"
# format is part of the stamp so a generator change does not reuse an old file.
FORMAT=1
# Compressed offsets past 2^31. The region file has to clear this.
REGION_MIN_BYTES=2147483648
# Products at or above 2^53 are not exact in awk's double, so positions would drift.
EXACT_INT_LIMIT=9007199254740992

CHROM_NAMES=(
  chr1A chr1B chr1D
  chr2A chr2B chr2D
  chr3A chr3B chr3D
  chr4A chr4B chr4D
  chr5A chr5B chr5D
  chr6A chr6B chr6D
  chr7A chr7B chr7D
)
CHROM_LENGTHS=(
  598660471 700547350 498638509
  787782082 812755788 656544405
  754128162 851934019 619618552
  754227511 673810255 518332611
  713360525 714805278 569951140
  622669697 731188232 495380293
  744491536 764081788 642921167
)

MODE=""
TAXA=""
SITES=""
SEED=""
OUTPUT=""

die() {
  printf 'generate-wheat-vcf: %s\n' "$*" >&2
  exit 1
}

file_bytes() {
  if stat -f%z "$1" >/dev/null 2>&1; then
    stat -f%z "$1"
  else
    stat -c%s "$1"
  fi
}

first_line() {
  local text="$1"
  printf '%s\n' "${text%%$'\n'*}"
}

require_htslib() {
  command -v bgzip >/dev/null || die "bgzip not found; brew install htslib"
  command -v bcftools >/dev/null || die "bcftools not found; brew install bcftools"
  local bgzip_line bcftools_out bcftools_first bcftools_second
  bgzip_line=$(first_line "$(bgzip --version 2>&1)") || die "bgzip --version failed"
  [[ "$bgzip_line" == "bgzip (htslib) $HTSLIB_VERSION" ]] ||
    die "bgzip must be htslib $HTSLIB_VERSION; found '$bgzip_line'"
  bcftools_out=$(bcftools --version 2>&1) || die "bcftools --version failed"
  bcftools_first=$(first_line "$bcftools_out")
  bcftools_second=$(first_line "${bcftools_out#*$'\n'}")
  [[ "$bcftools_first" == "bcftools $HTSLIB_VERSION" && "$bcftools_second" == "Using htslib $HTSLIB_VERSION" ]] ||
    die "bcftools must be $HTSLIB_VERSION built with htslib $HTSLIB_VERSION; found '$bcftools_first' / '$bcftools_second'"
}

usage() {
  sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    plain|region)
      [[ -z "$MODE" ]] || die "mode given twice"
      MODE="$1"
      shift
      ;;
    --taxa) TAXA="$2"; shift 2 ;;
    --sites) SITES="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done

[[ "$MODE" == plain || "$MODE" == region ]] || die "first argument must be plain or region"
[[ "$TAXA" =~ ^[0-9]+$ && "$TAXA" -ge 1 ]] || die "--taxa must be a positive integer"
[[ "$SITES" =~ ^[0-9]+$ && "$SITES" -ge 1 ]] || die "--sites must be a positive integer"
[[ "$SEED" =~ ^[0-9]+$ && "$SEED" -ge 1 && "$SEED" -lt 2147483647 ]] ||
  die "--seed must be in 1..2147483646 (0 locks the generator)"
[[ -n "$OUTPUT" ]] || die "--output is required"
if [[ "$MODE" == plain ]]; then
  [[ "$OUTPUT" == *.vcf && "$OUTPUT" != *.vcf.gz ]] || die "plain --output must end in .vcf"
else
  [[ "$OUTPUT" == *.vcf.gz ]] || die "region --output must end in .vcf.gz"
fi

# Largest-remainder split. One site is reserved for every chromosome first, so a
# 21-site file still touches chr6D. The rest follow length.
assign_counts() {
  local n=${#CHROM_NAMES[@]}
  local total=0 len extra assigned left best i add rem sum prod
  (( SITES >= n )) || die "--sites $SITES is below $n, the number of wheat chromosomes"
  for len in "${CHROM_LENGTHS[@]}"; do
    total=$((total + len))
  done
  extra=$((SITES - n))
  assigned=0
  COUNT_SPEC=""
  COUNTS=()
  REMAINDERS=()
  for i in "${!CHROM_LENGTHS[@]}"; do
    len=${CHROM_LENGTHS[$i]}
    add=$((extra * len / total))
    rem=$((extra * len % total))
    COUNTS[$i]=$((1 + add))
    REMAINDERS[$i]=$rem
    assigned=$((assigned + add))
  done
  left=$((extra - assigned))
  USED=()
  for i in "${!CHROM_LENGTHS[@]}"; do
    USED[$i]=0
  done
  while (( left > 0 )); do
    best=-1
    for i in "${!CHROM_LENGTHS[@]}"; do
      (( USED[i] )) && continue
      if (( best < 0 )); then
        best=$i
        continue
      fi
      if (( REMAINDERS[i] > REMAINDERS[best] || (REMAINDERS[i] == REMAINDERS[best] && i < best) )); then
        best=$i
      fi
    done
    USED[$best]=1
    COUNTS[$best]=$((COUNTS[best] + 1))
    left=$((left - 1))
  done
  sum=0
  for i in "${!COUNTS[@]}"; do
    (( COUNTS[i] >= 1 && COUNTS[i] <= CHROM_LENGTHS[i] )) ||
      die "${CHROM_NAMES[$i]} was given ${COUNTS[$i]} sites"
    if (( COUNTS[i] > 1 )); then
      prod=$(( (COUNTS[i] - 1) * (CHROM_LENGTHS[i] - 1) ))
      (( prod < EXACT_INT_LIMIT )) ||
        die "${CHROM_NAMES[$i]} positions would not be exact integers; lower --sites"
    fi
    sum=$((sum + COUNTS[i]))
    COUNT_SPEC+="${CHROM_NAMES[$i]}:${CHROM_LENGTHS[$i]}:${COUNTS[$i]},"
  done
  (( sum == SITES )) || die "allocated $sum sites, expected $SITES"
}

stamp_text() {
  printf 'format=%s\ntaxa=%s\nsites=%s\nseed=%s\n' "$FORMAT" "$TAXA" "$SITES" "$SEED"
}

outputs_ready() {
  local bytes
  if [[ "$MODE" == plain ]]; then
    [[ -f "$OUTPUT" && -f "${OUTPUT}.gz" ]] || return 1
  else
    [[ -f "$OUTPUT" && -f "${OUTPUT}.csi" ]] || return 1
    bytes=$(file_bytes "$OUTPUT")
    (( bytes > REGION_MIN_BYTES )) || return 1
  fi
}

# Writes the VCF text to stdout. Header, then one sorted record per site.
write_vcf() {
  awk -v taxa="$TAXA" -v seed="$SEED" -v spec="$COUNT_SPEC" '
    function emit(chrom, pos,    i, nbuf, buf, r) {
      printf "%s\t%d\t.\tA\tC\t.\tPASS\t.\tGT", chrom, pos
      buf = ""
      nbuf = 0
      for (i = 0; i < taxa; i++) {
        seed = (seed * 48271) % 2147483647
        r = int(seed % 5)
        if (r == 0) buf = buf "\t0/0"
        else if (r == 1) buf = buf "\t0/1"
        else if (r == 2) buf = buf "\t1/1"
        else if (r == 3) buf = buf "\t0|1"
        else buf = buf "\t./."
        nbuf++
        if (nbuf == 40) {
          printf "%s", buf
          buf = ""
          nbuf = 0
        }
      }
      printf "%s\n", buf
      emitted++
      if (emitted % 100000 == 0) printf "generate-wheat-vcf: %d sites\n", emitted > "/dev/stderr"
    }
    BEGIN {
      n = split(spec, rec, ",")
      if (rec[n] == "") n--
      printf "##fileformat=VCFv4.2\n"
      for (i = 1; i <= n; i++) {
        split(rec[i], p, ":")
        name[i] = p[1]
        len[i] = p[2] + 0
        cnt[i] = p[3] + 0
        printf "##contig=<ID=%s,length=%d>\n", name[i], len[i]
      }
      printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
      for (t = 1; t <= taxa; t++) printf "\tt%d", t
      printf "\n"
      emitted = 0
      for (i = 1; i <= n; i++) {
        if (cnt[i] == 1) {
          emit(name[i], 1)
          continue
        }
        span = len[i] - 1
        steps = cnt[i] - 1
        if (steps * span >= 9007199254740992) {
          printf "position product exceeds exact integer range\n" > "/dev/stderr"
          exit 2
        }
        for (s = 0; s < cnt[i]; s++) emit(name[i], 1 + int(s * span / steps))
      }
      if (emitted % 100000 != 0) printf "generate-wheat-vcf: %d sites\n", emitted > "/dev/stderr"
    }
  '
}

require_htslib
assign_counts

STAMP="${OUTPUT}.stamp"
mkdir -p "$(dirname "$OUTPUT")"

if [[ -f "$STAMP" ]] && [[ "$(cat "$STAMP")" == "$(stamp_text)" ]] && outputs_ready; then
  printf 'generate-wheat-vcf: reusing %s\n' "$OUTPUT"
  exit 0
fi

rm -f "$STAMP"
if [[ "$MODE" == plain ]]; then
  printf 'generate-wheat-vcf: writing %s taxa x %s sites\n' "$TAXA" "$SITES"
  write_vcf >"${OUTPUT}.part"
  mv "${OUTPUT}.part" "$OUTPUT"
  printf 'generate-wheat-vcf: bgzip -c -@ 1 %s\n' "$OUTPUT"
  bgzip -c -@ 1 "$OUTPUT" >"${OUTPUT}.gz.part"
  mv "${OUTPUT}.gz.part" "${OUTPUT}.gz"
else
  printf 'generate-wheat-vcf: streaming %s taxa x %s sites into %s\n' "$TAXA" "$SITES" "$OUTPUT"
  write_vcf | bgzip -c -@ 1 >"${OUTPUT}.part"
  mv "${OUTPUT}.part" "$OUTPUT"
  bytes=$(file_bytes "$OUTPUT")
  (( bytes > REGION_MIN_BYTES )) ||
    die "$(basename "$OUTPUT") is $bytes bytes; need more than $REGION_MIN_BYTES. Raise --sites."
  rm -f "${OUTPUT}.csi"
  bcftools index -c "$OUTPUT"
fi
stamp_text >"$STAMP"
printf 'generate-wheat-vcf: wrote %s\n' "$OUTPUT"
