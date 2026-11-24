#!/usr/bin/env bash
# Generate the htslib 1.25 oracle in testdata/htslib/ (ADR-0021, tassel-6#80;
# BCF and BGZF layout fixtures, ADR-0022, tassel-6#118; BCF encodings of the
# TASSEL 5 VCF-export goldens, tassel-6#129).
#
# Rebuilds testdata/htslib/ and its MANIFEST from scratch. Output is deterministic:
# running the script twice gives byte-identical files. Nothing under testdata/golden/
# or testdata/tutorial/ is modified.
#
# Requires bgzip, tabix, and bcftools built on htslib 1.25, Apple gzip 487.0.1
# (the plain-gzip fixtures are pinned to that binary), python3 (the edge-case BCF
# type check only; not a Gradle dependency), and xxd, od, dd, and awk.
# This is the maintainer's Mac (`brew install htslib bcftools`). Do not regenerate
# the plain-gzip fixtures with GNU gzip.
# bcftools view embeds Date= in the header unless --no-version is set. Every view
# in the BCF section passes it, or two runs would not match.
#
# Usage: scripts/generate-htslib-fixtures.sh
set -euo pipefail

HTSLIB_VERSION="1.25"
GZIP_PIN="Apple gzip 487.0.1"
SEED=80
FILLER_COUNT=32000
# A BGZF block holds at most 65280 uncompressed bytes. A contiguous run of records
# longer than two blocks cannot fit in two blocks, so the query spans several.
SPAN_BYTES=130560
# IWGSC RefSeq v2.1 pseudomolecule lengths, Table 4 "Length (bp)" of Zhu et al.,
# The Plant Genome 14, e20112 (2021), https://pmc.ncbi.nlm.nih.gov/articles/PMC8360199/
# chr6D is declared and then left with no records.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TESTDATA="$ROOT/testdata"
HTSLIB="$TESTDATA/htslib"
SCRATCH="$ROOT/build/htslib-fixtures"
MANIFEST="$HTSLIB/MANIFEST"

die() {
  printf 'generate-htslib-fixtures: %s\n' "$*" >&2
  exit 1
}

manifest() { printf '%s\n' "$@" >>"$MANIFEST"; }

note() { manifest "# $*"; }

first_line() {
  local text="$1"
  printf '%s\n' "${text%%$'\n'*}"
}

# Unsigned little-endian int32s at the CSI header. Prints "min_shift depth".
csi_shift_depth() {
  bgzip -dc "$1" | dd bs=1 count=12 2>/dev/null | od -An -tu1 | awk '
    { for (i = 1; i <= NF; i++) b[++n] = $i }
    END {
      if (n < 12) exit 2
      if (b[1] != 67 || b[2] != 83 || b[3] != 73 || b[4] != 1) exit 3
      shift = b[5] + b[6] * 256 + b[7] * 65536 + b[8] * 16777216
      depth = b[9] + b[10] * 256 + b[11] * 65536 + b[12] * 16777216
      printf "%d %d\n", shift, depth
    }'
}

# Prints "members eof_blocks" for a file that is a sequence of BGZF blocks.
# Walking the BSIZE field avoids treating a copy of the gzip magic inside a
# deflate stream as another member. bgzip 1.25 text mode also flushes after a
# VCF header, so a small VCF is a header block, a body block, and an EOF block.
walk_bgzf() {
  od -An -tu1 -v "$1" | awk '
    {
      for (i = 1; i <= NF; i++) raw[++n] = $i + 0
    }
    END {
      off = 1
      members = 0
      eofs = 0
      while (off <= n) {
        if (off + 17 > n) exit 2
        if (raw[off] != 31 || raw[off + 1] != 139 || raw[off + 2] != 8 || raw[off + 3] != 4)
          exit 3
        bsize = raw[off + 16] + raw[off + 17] * 256
        size = bsize + 1
        if (size < 28 || off + size - 1 > n) exit 4
        if (size == 28) eofs++
        off += size
        members++
      }
      if (off != n + 1) exit 5
      printf "%d %d\n", members, eofs
    }
  '
}

# One decimal on-disk member length per line (BSIZE+1), including the EOF member.
# The walk is the one walk_bgzf already uses. Compression level can change these
# lengths; do not substitute the uncompressed ISIZE.
bgzf_member_sizes() {
  od -An -tu1 -v "$1" | awk '
    {
      for (i = 1; i <= NF; i++) raw[++n] = $i + 0
    }
    END {
      off = 1
      while (off <= n) {
        if (off + 17 > n) exit 2
        if (raw[off] != 31 || raw[off + 1] != 139 || raw[off + 2] != 8 || raw[off + 3] != 4)
          exit 3
        bsize = raw[off + 16] + raw[off + 17] * 256
        size = bsize + 1
        if (size < 28 || off + size - 1 > n) exit 4
        print size
        off += size
      }
      if (off != n + 1) exit 5
    }
  '
}

# Lowercase SHA-256 hex. sha256sum on Linux, shasum on the maintainer's Mac.
sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

# Data blocks of a BGZF file that ends with exactly one EOF block.
bgzf_data_blocks() {
  local got members eofs
  got=$(walk_bgzf "$1") || return 1
  members=${got%% *}
  eofs=${got##* }
  [[ "$eofs" == "1" ]] || return 1
  echo $((members - 1))
}

# Little-endian uint32 at a byte offset.
le_u32_at() {
  dd bs=1 skip="$2" count=4 <"$1" 2>/dev/null | od -An -tu1 | awk '
    { for (i = 1; i <= NF; i++) if (n < 4) v[++n] = $i }
    END { print v[1] + v[2] * 256 + v[3] * 65536 + v[4] * 16777216 }
  '
}

# l_aux from a CSI header. A BCF CSI is 0: the index stores no contig names.
csi_l_aux() {
  bgzip -dc "$1" >"$SCRATCH/csi.plain"
  le_u32_at "$SCRATCH/csi.plain" 12
}

# Contig names from a VCF CSI aux block. The aux starts with the tabix preset
# (six int32s) and a name-blob length, then NUL-terminated names.
csi_contig_names() {
  local plain="$SCRATCH/csi.plain" l_aux names_len
  bgzip -dc "$1" >"$plain"
  l_aux=$(le_u32_at "$plain" 12)
  [[ "$l_aux" == "0" ]] && return 0
  names_len=$(le_u32_at "$plain" 40)
  dd bs=1 skip=44 count="$names_len" <"$plain" 2>/dev/null | tr '\0' ' ' | awk '{
    sub(/ +$/, "")
    print
  }'
}

# Space-separated unsigned bytes at the start of a file.
prefix_bytes() {
  dd bs=1 count="$2" <"$1" 2>/dev/null | od -An -tu1 | awk '{
    for (i = 1; i <= NF; i++) printf "%s%s", (n++ ? " " : ""), $i
  }'
}

# One little-endian byte of a file.
le_byte() {
  dd bs=1 skip="$2" count=1 <"$1" 2>/dev/null | od -An -tu1 | awk '
    { for (i = 1; i <= NF; i++) if (n < 1) v[++n] = $i }
    END { print v[1] }
  '
}

# Read min_shift and depth from a CSI and record them. Depth must be expect.
assert_csi() {
  local csi="$1" expect="$2" label="$3"
  local got shift depth
  got=$(csi_shift_depth "$csi") || die "could not read the CSI header of $label"
  shift=${got%% *}
  depth=${got##* }
  [[ -n "$shift" && -n "$depth" ]] || die "empty CSI header for $label"
  [[ "$depth" == "$expect" ]] || die "$label CSI depth is $depth, expected $expect (min_shift $shift)"
  manifest "csi_min_shift = $shift" "csi_depth = $depth"
}

# tabix must refuse a file whose coordinates do not fit in a TBI. The matching
# stderr line is left in TABIX_LINE. Any partial .tbi is removed.
TABIX_LINE=
refuse_tabix() {
  local gz="$1"
  local err="$SCRATCH/tabix.err"
  local status line
  set +e
  tabix -p vcf "$gz" >"$SCRATCH/tabix.out" 2>"$err"
  status=$?
  set -e
  rm -f "${gz}.tbi"
  [[ "$status" -ne 0 ]] || die "tabix indexed $gz; a coordinate past 2^29 must not fit in a TBI"
  line=$(grep -F "cannot be stored in a tbi index" "$err" || true)
  line=${line%%$'\n'*}
  [[ -n "$line" ]] || die "tabix refused $gz without the expected message: $(cat "$err")"
  TABIX_LINE=$line
}

# region_query GZ REGION DEST_DIR EXPECT SPAN STDERR
# EXPECT is "empty" or "records". SPAN is "span" to require the record-mode
# output to cover more than two BGZF blocks. STDERR is "allow-stderr" for the
# tutorial, whose header makes htslib warn about PL and GQ.
region_query() {
  local gz="$1" region="$2" dest_dir="$3" expect="$4" span="$5" stderr_mode="$6"
  local mode slug dest err status bytes
  slug=$(printf '%s' "$region" | tr ':-' '__')
  err="$SCRATCH/view.err"
  for mode in record pos; do
    dest="$dest_dir/${slug}.${mode}.txt"
    set +e
    bcftools view -H --regions-overlap "$mode" -r "$region" "$gz" >"$dest" 2>"$err"
    status=$?
    set -e
    [[ "$status" -eq 0 ]] || die "bcftools view -r $region --regions-overlap $mode exited $status: $(cat "$err")"
    if [[ "$stderr_mode" != "allow-stderr" && -s "$err" ]]; then
      die "bcftools view -r $region ($mode) wrote to stderr: $(cat "$err")"
    fi
    if [[ "$expect" == "empty" && -s "$dest" ]]; then
      die "expected no records for $region ($mode)"
    fi
    if [[ "$span" == "span" && "$mode" == "record" ]]; then
      bytes=$(wc -c <"$dest" | tr -d ' ')
      [[ "$bytes" -gt "$SPAN_BYTES" ]] || die "region $region is $bytes bytes; it does not span several BGZF blocks"
    fi
    manifest "command = bcftools view -H --regions-overlap $mode -r $region $gz"
    manifest "file = $dest"
    REGION_QUERIES+=("${gz}"$'\t'"${region}"$'\t'"${mode}"$'\t'"${dest}")
  done
}

# bcf_run LABEL MODE DEST -- command...
# MODE is "strict" (stderr must be empty) or "require-sanity" (only the PL/GQ warning).
# DEST is a file, or "-" to leave stdout alone.
bcf_run() {
  local label="$1" mode="$2" dest="$3"
  shift 3
  local err="$SCRATCH/bcf.err" status other
  set +e
  if [[ "$dest" == "-" ]]; then
    "$@" 2>"$err"
  else
    "$@" >"$dest" 2>"$err"
  fi
  status=$?
  set -e
  [[ "$status" -eq 0 ]] || die "$label exited $status: $(cat "$err")"
  case "$mode" in
    strict)
      [[ ! -s "$err" ]] || die "$label wrote to stderr: $(cat "$err")"
      ;;
    require-sanity)
      grep -q "bcf_hdr_check_sanity" "$err" ||
        die "$label did not warn about PL/GQ: $(cat "$err")"
      other=$(grep -v "bcf_hdr_check_sanity" "$err" || true)
      [[ -z "$other" ]] || die "$label wrote unexpected stderr: $(cat "$err")"
      ;;
    *) die "bad stderr mode $mode" ;;
  esac
}

# Read min_shift, depth, and l_aux. Die unless they equal the expected triple.
# Records prefix_csi_min_shift, prefix_csi_depth, and prefix_csi_l_aux.
record_csi() {
  local csi="$1" prefix="$2" expect_shift="$3" expect_depth="$4" expect_aux="$5"
  local got shift depth l_aux
  got=$(csi_shift_depth "$csi") || die "could not read the CSI header of $prefix"
  shift=${got%% *}
  depth=${got##* }
  l_aux=$(csi_l_aux "$csi")
  [[ "$shift" == "$expect_shift" && "$depth" == "$expect_depth" && "$l_aux" == "$expect_aux" ]] ||
    die "$prefix CSI is min_shift $shift depth $depth l_aux $l_aux, expected $expect_shift $expect_depth $expect_aux"
  manifest \
    "${prefix}_csi_min_shift = $shift" \
    "${prefix}_csi_depth = $depth" \
    "${prefix}_csi_l_aux = $l_aux"
}

# --- Toolchain ----------------------------------------------------------------

export LC_ALL=C

command -v bgzip >/dev/null || die "bgzip not found; brew install htslib"
command -v tabix >/dev/null || die "tabix not found; brew install htslib"
command -v bcftools >/dev/null || die "bcftools not found; brew install bcftools"
command -v gzip >/dev/null || die "gzip not found"
command -v xxd >/dev/null || die "xxd not found"
command -v od >/dev/null || die "od not found"
command -v python3 >/dev/null || die "python3 not found; the edge-case BCF type check needs it"

bgzip_out=$(bgzip --version 2>&1) || die "bgzip --version failed"
bgzip_line=$(first_line "$bgzip_out")
[[ "$bgzip_line" == "bgzip (htslib) $HTSLIB_VERSION" ]] ||
  die "bgzip must be htslib $HTSLIB_VERSION; found '$bgzip_line'"

tabix_out=$(tabix --version 2>&1) || die "tabix --version failed"
tabix_line=$(first_line "$tabix_out")
[[ "$tabix_line" == "tabix (htslib) $HTSLIB_VERSION" ]] ||
  die "tabix must be htslib $HTSLIB_VERSION; found '$tabix_line'"

bcftools_out=$(bcftools --version 2>&1) || die "bcftools --version failed"
bcftools_first=$(first_line "$bcftools_out")
bcftools_rest=${bcftools_out#*$'\n'}
bcftools_second=$(first_line "$bcftools_rest")
[[ "$bcftools_first" == "bcftools $HTSLIB_VERSION" && "$bcftools_second" == "Using htslib $HTSLIB_VERSION" ]] ||
  die "bcftools must be $HTSLIB_VERSION built with htslib $HTSLIB_VERSION; found '$bcftools_first' / '$bcftools_second'"

gzip_out=$(gzip --version 2>&1) || die "gzip --version failed"
gzip_line=$(first_line "$gzip_out")
[[ "$gzip_line" == "$GZIP_PIN" ]] ||
  die "plain gzip fixtures are pinned to $GZIP_PIN; found '$gzip_line'"

# --- Output directory ----------------------------------------------------------

rm -rf "$HTSLIB" "$SCRATCH"
mkdir -p "$SCRATCH" \
  "$HTSLIB/tutorial/regions" \
  "$HTSLIB/wheat/regions" \
  "$HTSLIB/large-pos/regions" \
  "$HTSLIB/variants"
: >"$MANIFEST"

cd "$TESTDATA"

manifest \
  "# htslib $HTSLIB_VERSION fixtures for compressed and indexed VCF (ADR-0021, tassel-6#80)" \
  "# and for BCF and BGZF layout (ADR-0022, tassel-6#118)." \
  "# [golden-bcf] is htslib's encoding of the TASSEL 5 VCF-export goldens (tassel-6#129)." \
  "# Generated by scripts/generate-htslib-fixtures.sh. Do not edit by hand." \
  "#" \
  "# Paths are relative to testdata/, and every command runs from there." \
  "#   file       may repeat; each value is one output file" \
  "#   command    may repeat; the invocation that produced the following file" \
  "#   identical  may repeat; bcftools view -H --regions-overlap on the BCF was" \
  "#              identical to this file, so the body is not stored again" \
  "# CSI min_shift and depth are read back from the index header." \
  "# Both bcftools region-overlap modes are kept. ADR-0021 chooses record for the" \
  "# reader; this oracle does not make that choice." \
  "#" \
  "[htslib]" \
  "version = $HTSLIB_VERSION" \
  "bgzip = $bgzip_line" \
  "tabix = $tabix_line" \
  "bcftools = $bcftools_first" \
  "bcftools_htslib = $bcftools_second" \
  "gzip = $gzip_line"

echo "Generating htslib $HTSLIB_VERSION fixtures"

REGION_QUERIES=()

# --- Tutorial ------------------------------------------------------------------

manifest "" "[tutorial]"
note "bgzip of tutorial/mdp_genotype_w_depth.vcf. That file has no ##contig lines."
note "htslib $HTSLIB_VERSION writes this CSI with depth 8. Readers must take min_shift and depth from the file."
note "bcftools index -s prints [W::bcf_hdr_check_sanity] warnings about PL and GQ on stderr; only stdout is kept."
note "Chromosome 1 of the tutorial is 3.18 MB of records. Dumping it whole, twice, would put the tree over the 10 MB cap,"
note "so the whole-contig oracle is wheat chr3B. The tutorial range below still crosses several BGZF blocks."

bgzip -c tutorial/mdp_genotype_w_depth.vcf >htslib/tutorial/mdp_genotype_w_depth.vcf.gz
manifest "command = bgzip -c tutorial/mdp_genotype_w_depth.vcf"
manifest "file = htslib/tutorial/mdp_genotype_w_depth.vcf.gz"

bcftools index -c htslib/tutorial/mdp_genotype_w_depth.vcf.gz
manifest "command = bcftools index -c htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
manifest "file = htslib/tutorial/mdp_genotype_w_depth.vcf.gz.csi"
assert_csi htslib/tutorial/mdp_genotype_w_depth.vcf.gz.csi 8 "tutorial"

tabix -p vcf htslib/tutorial/mdp_genotype_w_depth.vcf.gz
manifest "command = tabix -p vcf htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
manifest "file = htslib/tutorial/mdp_genotype_w_depth.vcf.gz.tbi"

bcftools index -s htslib/tutorial/mdp_genotype_w_depth.vcf.gz \
  >htslib/tutorial/index-s.txt 2>"$SCRATCH/index-s.err" ||
  die "bcftools index -s failed on the tutorial: $(cat "$SCRATCH/index-s.err")"
grep -q "bcf_hdr_check_sanity" "$SCRATCH/index-s.err" ||
  die "expected the tutorial PL/GQ sanity warning on stderr"
[[ -s htslib/tutorial/index-s.txt ]] || die "tutorial index -s wrote no stdout"
manifest "command = bcftools index -s htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
manifest "file = htslib/tutorial/index-s.txt"

tabix -l htslib/tutorial/mdp_genotype_w_depth.vcf.gz >htslib/tutorial/contigs.txt
manifest "command = tabix -l htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
manifest "file = htslib/tutorial/contigs.txt"

tutorial_blocks=$(bgzf_data_blocks htslib/tutorial/mdp_genotype_w_depth.vcf.gz) ||
  die "could not count tutorial BGZF blocks"
manifest "bgzf_data_blocks = $tutorial_blocks"

manifest "" "[tutorial-regions]"
note "Queries are beg-end, and a single position is pos-pos, so this oracle does not depend on how chr:pos is parsed."
note "Stderr is discarded: bcftools view repeats the PL/GQ header warning."
region_query htslib/tutorial/mdp_genotype_w_depth.vcf.gz "1:157104-157104" \
  htslib/tutorial/regions records nospan allow-stderr
region_query htslib/tutorial/mdp_genotype_w_depth.vcf.gz "1:157104-10068726" \
  htslib/tutorial/regions records span allow-stderr
region_query htslib/tutorial/mdp_genotype_w_depth.vcf.gz "1:299170078-299999999" \
  htslib/tutorial/regions empty nospan allow-stderr
region_query htslib/tutorial/mdp_genotype_w_depth.vcf.gz "99:1-1000" \
  htslib/tutorial/regions empty nospan allow-stderr
region_query htslib/tutorial/mdp_genotype_w_depth.vcf.gz "1:10068726-157104" \
  htslib/tutorial/regions empty nospan allow-stderr
manifest "unknown_contig_exit = 0" "inverted_region_exit = 0"

# --- Wheat ---------------------------------------------------------------------

manifest "" "[wheat]"
note "Synthetic wheat-like VCF. Contig lengths are IWGSC RefSeq v2.1 Table 4 Length (bp),"
note "Zhu et al., The Plant Genome 14, e20112 (2021), PMC8360199."
note "Names are chr1A through chr7D. chr6D is declared and has no records."
note "bcftools index -s omits a contig with no records, so chr6D is not in that table."
note "Four samples. Filler SNPs are on chr3B from a fixed LCG seed so the file spans 20 to 50 BGZF blocks."
note "A sequence deletion at chr3B:536870000 (ACGT to A) is the ADR-0021 example."
note "One <DEL> carries INFO/END across several 16384 bp CSI bins. Another <DEL> carries only INFO/SVLEN."
note "tabix -p vcf refuses the file because a coordinate past 2^29 cannot be stored in a TBI. There is no .tbi."
manifest \
  "lengths_source = Zhu et al., The Plant Genome 14, e20112 (2021), PMC8360199, Table 4, Length (bp)" \
  "empty_contig = chr6D" \
  "samples = s1 s2 s3 s4" \
  "seed = $SEED" \
  "filler_count = $FILLER_COUNT"

cat >"$SCRATCH/wheat.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1A,length=598660471>
##contig=<ID=chr1B,length=700547350>
##contig=<ID=chr1D,length=498638509>
##contig=<ID=chr2A,length=787782082>
##contig=<ID=chr2B,length=812755788>
##contig=<ID=chr2D,length=656544405>
##contig=<ID=chr3A,length=754128162>
##contig=<ID=chr3B,length=851934019>
##contig=<ID=chr3D,length=619618552>
##contig=<ID=chr4A,length=754227511>
##contig=<ID=chr4B,length=673810255>
##contig=<ID=chr4D,length=518332611>
##contig=<ID=chr5A,length=713360525>
##contig=<ID=chr5B,length=714805278>
##contig=<ID=chr5D,length=569951140>
##contig=<ID=chr6A,length=622669697>
##contig=<ID=chr6B,length=731188232>
##contig=<ID=chr6D,length=495380293>
##contig=<ID=chr7A,length=744491536>
##contig=<ID=chr7B,length=764081788>
##contig=<ID=chr7D,length=642921167>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n' >>"$SCRATCH/wheat.vcf"

awk -v seed="$SEED" -v count="$FILLER_COUNT" '
  function rnd() {
    seed = (seed * 48271) % 2147483647
    return seed
  }
  function gt(    r) {
    r = int(rnd() % 5)
    if (r == 0) return "0/0"
    if (r == 1) return "0/1"
    if (r == 2) return "1/1"
    if (r == 3) return "1/0"
    return "0|1"
  }
  function snp(chrom, pos) {
    printf "%s\t%d\t.\tA\tC\t.\tPASS\t.\tGT\t%s\t%s\t%s\t%s\n", chrom, pos, gt(), gt(), gt(), gt()
  }
  function insertion() {
    printf "chr3B\t50000\t.\tA\tACG\t.\tPASS\t.\tGT\t0|1\t0/1\t1/1\t0/0\n"
  }
  function multi() {
    printf "chr3B\t60000\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/1\t1/2\t0/2\t0/0\n"
  }
  function missing() {
    printf "chr3B\t70000\t.\tA\tC\t.\tPASS\t.\tGT\t./.\t0/1\t0/0\t1/1\n"
  }
  function phased() {
    printf "chr3B\t80000\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\t0|1\t0/0\n"
  }
  function del_end() {
    printf "chr3B\t163840\t.\tN\t<DEL>\t.\tPASS\tEND=229376;SVTYPE=DEL\tGT\t0/1\t0/0\t0/0\t0/1\n"
  }
  function deletion() {
    printf "chr3B\t536870000\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
  }
  function svlen() {
    printf "chr3B\t600000000\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-40000;SVTYPE=DEL\tGT\t0/1\t0/0\t0/1\t0/0\n"
  }
  BEGIN {
    snp("chr1A", 590000000)
    snp("chr2A", 600000000)
    end_done = 0
    for (i = 1; i <= count; i++) {
      pos = i * 1000
      if (!end_done && pos > 163840) {
        del_end()
        end_done = 1
      }
      if (pos == 50000) { insertion(); continue }
      if (pos == 60000) { multi(); continue }
      if (pos == 70000) { missing(); continue }
      if (pos == 80000) { phased(); continue }
      snp("chr3B", pos)
    }
    if (!end_done) del_end()
    deletion()
    svlen()
    snp("chr3B", 700000000)
    snp("chr7B", 600000000)
    snp("chr7D", 600000000)
  }
' >>"$SCRATCH/wheat.vcf"

bgzip -c "$SCRATCH/wheat.vcf" >htslib/wheat/wheat.vcf.gz
manifest "command = bgzip -c (synthetic wheat VCF, seed $SEED, $FILLER_COUNT chr3B filler sites)"
manifest "file = htslib/wheat/wheat.vcf.gz"

wheat_blocks=$(bgzf_data_blocks htslib/wheat/wheat.vcf.gz) || die "could not count wheat BGZF blocks"
[[ "$wheat_blocks" -ge 20 && "$wheat_blocks" -le 50 ]] ||
  die "wheat has $wheat_blocks BGZF data blocks; expected 20 to 50"
manifest "bgzf_data_blocks = $wheat_blocks"

bcftools index -c htslib/wheat/wheat.vcf.gz
manifest "command = bcftools index -c htslib/wheat/wheat.vcf.gz"
manifest "file = htslib/wheat/wheat.vcf.gz.csi"
assert_csi htslib/wheat/wheat.vcf.gz.csi 6 "wheat"

refuse_tabix htslib/wheat/wheat.vcf.gz
manifest "command = tabix -p vcf htslib/wheat/wheat.vcf.gz"
manifest "tabix_refused = $TABIX_LINE"

bcftools index -s htslib/wheat/wheat.vcf.gz >htslib/wheat/index-s.txt 2>"$SCRATCH/index-s.err" ||
  die "bcftools index -s failed on wheat: $(cat "$SCRATCH/index-s.err")"
[[ ! -s "$SCRATCH/index-s.err" ]] || die "bcftools index -s wrote to stderr for wheat: $(cat "$SCRATCH/index-s.err")"
wheat_contigs=$(wc -l <htslib/wheat/index-s.txt | tr -d ' ')
[[ "$wheat_contigs" == "5" ]] ||
  die "wheat index -s has $wheat_contigs lines, expected the 5 contigs that have records"
for contig in chr1A chr2A chr3B chr7B chr7D; do
  grep -q "^${contig}"$'\t' htslib/wheat/index-s.txt || die "wheat index -s is missing $contig"
done
grep -q '^chr6D'$'\t' htslib/wheat/index-s.txt &&
  die "index -s listed chr6D; htslib 1.25 omits a contig with no records"
grep -q '##contig=<ID=chr6D,length=495380293>' "$SCRATCH/wheat.vcf" ||
  die "chr6D is missing from the wheat VCF header"
manifest "index_s_omits_empty_contigs = true"
manifest "command = bcftools index -s htslib/wheat/wheat.vcf.gz"
manifest "file = htslib/wheat/index-s.txt"
note "No tabix -l output: the file has no TBI. Contig names and counts are the index -s table."

manifest "" "[wheat-regions]"
note "A region that starts inside the sequence deletion, the END deletion, or the SVLEN deletion"
note "is returned for --regions-overlap record and omitted for --regions-overlap pos."
note "An unknown contig and an inverted region print nothing and exit 0."
region_query htslib/wheat/wheat.vcf.gz "chr3B:1000-1000" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:1000-5000000" \
  htslib/wheat/regions records span strict
region_query htslib/wheat/wheat.vcf.gz "chr3B" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:599000000-601000000" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:700000000-700000000" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:536870002-536870010" \
  htslib/wheat/regions records nospan strict
# 180100-180200 sits inside the END deletion and between filler SNPs, which are every 1000 bp.
region_query htslib/wheat/wheat.vcf.gz "chr3B:180100-180200" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:600010000-600020000" \
  htslib/wheat/regions records nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:800000000-851934019" \
  htslib/wheat/regions empty nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr6D:1-495380293" \
  htslib/wheat/regions empty nospan strict
region_query htslib/wheat/wheat.vcf.gz "chrUn:1-1000" \
  htslib/wheat/regions empty nospan strict
region_query htslib/wheat/wheat.vcf.gz "chr3B:5000-1000" \
  htslib/wheat/regions empty nospan strict
manifest "unknown_contig_exit = 0" "inverted_region_exit = 0"

# record overlap returns the spanning alleles; pos overlap does not.
for slug in chr3B_536870002_536870010 chr3B_180100_180200 chr3B_600010000_600020000; do
  [[ -s "htslib/wheat/regions/${slug}.record.txt" ]] || die "$slug record overlap returned nothing"
  [[ ! -s "htslib/wheat/regions/${slug}.pos.txt" ]] || die "$slug pos overlap returned a record"
done
[[ -s htslib/wheat/regions/chr3B_700000000_700000000.record.txt ]] ||
  die "the site past 2^29 was not returned"
[[ -s htslib/wheat/regions/chr3B.record.txt ]] || die "whole chr3B returned nothing"

# --- Position past Int.MAX_VALUE ----------------------------------------------

manifest "" "[large-pos]"
note "One contig of length 4000000000 and one record at 3000000000, past Int.MAX_VALUE."
note "htslib 1.25 indexes it (CSI depth 6) and answers a region query. The GenotypeTable boundary is what cannot hold the position."
note "tabix refuses it. There is no .tbi, so there is no tabix -l output."

{
  printf '%s\n' '##fileformat=VCFv4.2'
  printf '%s\n' '##contig=<ID=chr1,length=4000000000>'
  printf '%s\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr1\t3000000000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n'
} >"$SCRATCH/large_pos.vcf"

bgzip -c "$SCRATCH/large_pos.vcf" >htslib/large-pos/large_pos.vcf.gz
manifest "command = bgzip -c (one record at 3000000000 on chr1 of length 4000000000)"
manifest "file = htslib/large-pos/large_pos.vcf.gz"

bcftools index -c htslib/large-pos/large_pos.vcf.gz
manifest "command = bcftools index -c htslib/large-pos/large_pos.vcf.gz"
manifest "file = htslib/large-pos/large_pos.vcf.gz.csi"
assert_csi htslib/large-pos/large_pos.vcf.gz.csi 6 "large-pos"

refuse_tabix htslib/large-pos/large_pos.vcf.gz
manifest "command = tabix -p vcf htslib/large-pos/large_pos.vcf.gz"
manifest "tabix_refused = $TABIX_LINE"

bcftools index -s htslib/large-pos/large_pos.vcf.gz >htslib/large-pos/index-s.txt 2>"$SCRATCH/index-s.err" ||
  die "bcftools index -s failed on large-pos: $(cat "$SCRATCH/index-s.err")"
[[ ! -s "$SCRATCH/index-s.err" ]] ||
  die "bcftools index -s wrote to stderr for large-pos: $(cat "$SCRATCH/index-s.err")"
manifest "command = bcftools index -s htslib/large-pos/large_pos.vcf.gz"
manifest "file = htslib/large-pos/index-s.txt"

manifest "" "[large-pos-regions]"
region_query htslib/large-pos/large_pos.vcf.gz "chr1:3000000000-3000000000" \
  htslib/large-pos/regions records nospan strict
[[ -s htslib/large-pos/regions/chr1_3000000000_3000000000.record.txt ]] ||
  die "the record at 3000000000 was not returned"

# --- Compression variants ------------------------------------------------------

manifest "" "[variants]"
note "One small VCF, then six compressions of those same bytes."
note "plain gzip is Apple gzip -n, not BGZF. Regenerating it with GNU gzip would not match."
note "concatenated is two gzip members. mid-eof is two bgzip members, so a 28-byte EOF block sits in the middle."
note "no-eof drops that trailing 28-byte block. truncated keeps the first 40 bytes, inside the first block."
note "bgzip 1.25 text mode flushes after the VCF header, so the bgzipped source is a header block plus a body block."
note "bad-crc flips the first CRC32 byte of the last data block, the four bytes immediately before that block ISIZE."

{
  printf '%s\n' '##fileformat=VCFv4.2'
  printf '%s\n' '##contig=<ID=chr1,length=1000>'
  printf '%s\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n'
  printf 'chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\n'
  printf 'chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\t./.\n'
} >htslib/variants/source.vcf
manifest "file = htslib/variants/source.vcf"

gzip -n -c htslib/variants/source.vcf >htslib/variants/plain.vcf.gz
plain_header=$(dd bs=1 count=4 <htslib/variants/plain.vcf.gz 2>/dev/null | od -An -tu1 | awk '
  { for (i = 1; i <= NF; i++) if (n < 4) v[++n] = $i }
  END { print v[1], v[2], v[3], v[4] }
')
[[ "$plain_header" == "31 139 8 0" ]] || die "plain gzip header is '$plain_header', expected gzip without FEXTRA"
manifest "command = gzip -n -c htslib/variants/source.vcf"
manifest "file = htslib/variants/plain.vcf.gz"

src_bytes=$(wc -c <htslib/variants/source.vcf | tr -d ' ')
split_at=$((src_bytes / 2))
[[ "$split_at" -gt 0 && "$split_at" -lt "$src_bytes" ]] || die "could not split source.vcf ($src_bytes bytes)"
dd bs=1 count="$split_at" <htslib/variants/source.vcf >"$SCRATCH/half1" 2>/dev/null
dd bs=1 skip="$split_at" <htslib/variants/source.vcf >"$SCRATCH/half2" 2>/dev/null
[[ -s "$SCRATCH/half1" && -s "$SCRATCH/half2" ]] || die "a split half is empty"

gzip -n -c "$SCRATCH/half1" >"$SCRATCH/g1"
gzip -n -c "$SCRATCH/half2" >"$SCRATCH/g2"
cat "$SCRATCH/g1" "$SCRATCH/g2" >htslib/variants/concatenated.vcf.gz
g1_bytes=$(wc -c <"$SCRATCH/g1" | tr -d ' ')
g2_bytes=$(wc -c <"$SCRATCH/g2" | tr -d ' ')
dd bs=1 count="$g1_bytes" <htslib/variants/concatenated.vcf.gz >"$SCRATCH/prefix" 2>/dev/null
dd bs=1 skip="$g1_bytes" <htslib/variants/concatenated.vcf.gz >"$SCRATCH/suffix" 2>/dev/null
cmp -s "$SCRATCH/g1" "$SCRATCH/prefix" || die "concatenated gzip prefix is not the first member"
cmp -s "$SCRATCH/g2" "$SCRATCH/suffix" || die "concatenated gzip suffix is not the second member"
[[ "$(le_byte "$SCRATCH/g1" 3)" != "4" && "$(le_byte "$SCRATCH/g2" 3)" != "4" ]] ||
  die "a concatenated gzip member is BGZF"
manifest "split_at = $split_at"
manifest "command = gzip -n -c (bytes 0..$((split_at - 1)) of htslib/variants/source.vcf)"
manifest "command = gzip -n -c (bytes ${split_at}..end of htslib/variants/source.vcf)"
manifest "command = cat"
manifest "file = htslib/variants/concatenated.vcf.gz"

bgzip -c "$SCRATCH/half1" >"$SCRATCH/b1"
bgzip -c "$SCRATCH/half2" >"$SCRATCH/b2"
cat "$SCRATCH/b1" "$SCRATCH/b2" >htslib/variants/mid-eof.vcf.gz
mid_walk=$(walk_bgzf htslib/variants/mid-eof.vcf.gz) || die "mid-eof is not a sequence of BGZF blocks"
[[ "${mid_walk##* }" == "2" ]] || die "mid-eof has ${mid_walk##* } EOF blocks, expected 2"
manifest "command = bgzip -c (bytes 0..$((split_at - 1)) of htslib/variants/source.vcf)"
manifest "command = bgzip -c (bytes ${split_at}..end of htslib/variants/source.vcf)"
manifest "command = cat"
manifest "file = htslib/variants/mid-eof.vcf.gz"

bgzip -c htslib/variants/source.vcf >"$SCRATCH/full.bgz"
full_walk=$(walk_bgzf "$SCRATCH/full.bgz") || die "could not walk the bgzipped variant source"
[[ "${full_walk##* }" == "1" && "${full_walk%% *}" -ge 2 ]] ||
  die "bgzipped variant source is '$full_walk', expected data blocks plus one EOF"
full_bytes=$(wc -c <"$SCRATCH/full.bgz" | tr -d ' ')
first_block=$(($(dd bs=1 skip=16 count=2 <"$SCRATCH/full.bgz" 2>/dev/null | od -An -tu1 | awk '
  { for (i = 1; i <= NF; i++) if (n < 2) v[++n] = $i }
  END { print v[1] + v[2] * 256 }
') + 1))

dd bs=1 count=$((full_bytes - 28)) <"$SCRATCH/full.bgz" >htslib/variants/no-eof.vcf.gz 2>/dev/null
no_walk=$(walk_bgzf htslib/variants/no-eof.vcf.gz) || die "no-eof is not a sequence of complete BGZF blocks"
[[ "${no_walk##* }" == "0" ]] || die "no-eof still has an EOF block"
manifest "command = bgzip -c htslib/variants/source.vcf"
manifest "command = drop the final 28-byte BGZF EOF block"
manifest "file = htslib/variants/no-eof.vcf.gz"

dd bs=1 count=40 <"$SCRATCH/full.bgz" >htslib/variants/truncated.vcf.gz 2>/dev/null
[[ "40" -gt 18 && "40" -lt "$first_block" ]] ||
  die "40 bytes is not inside the first BGZF block (block size $first_block)"
manifest "command = bgzip -c htslib/variants/source.vcf"
manifest "command = keep the first 40 bytes"
manifest "file = htslib/variants/truncated.vcf.gz"

crc_off=$((full_bytes - 28 - 8))
crc_byte=$(le_byte "$SCRATCH/full.bgz" "$crc_off")
flipped=$((crc_byte ^ 1))
{
  dd bs=1 count="$crc_off" <"$SCRATCH/full.bgz" 2>/dev/null
  printf "\\$(printf '%03o' "$flipped")"
  dd bs=1 skip=$((crc_off + 1)) <"$SCRATCH/full.bgz" 2>/dev/null
} >htslib/variants/bad-crc.vcf.gz
bad_bytes=$(wc -c <htslib/variants/bad-crc.vcf.gz | tr -d ' ')
[[ "$bad_bytes" == "$full_bytes" ]] || die "bad-crc size $bad_bytes != $full_bytes"
cmp -s "$SCRATCH/full.bgz" htslib/variants/bad-crc.vcf.gz && die "CRC flip did not change the file"
manifest "command = bgzip -c htslib/variants/source.vcf"
manifest "crc_byte_offset = $crc_off"
manifest "command = xor the first CRC32 byte of the last data block with 1"
manifest "file = htslib/variants/bad-crc.vcf.gz"

# --- BCF -----------------------------------------------------------------------

mkdir -p \
  htslib/bcf/tutorial \
  htslib/bcf/wheat \
  htslib/bcf/small \
  htslib/bcf/boundary \
  htslib/bcf/edge \
  htslib/bcf/corrupt \
  htslib/bgzf-layout

manifest "" "[bcf]"
note "BCF 2.2 from bcftools view --no-version -Ob of the bgzipped VCF, then bcftools index -c."
note "The plain tutorial VCF has no ##contig lines and bcftools view -Ob refuses it, so the .vcf.gz is the input."
note "That read adds ##contig lines for chromosomes 1 through 10. The source VCF header has none."
note "bcftools view without --no-version embeds Date= in the header, so every view below passes --no-version."
note "A BCF CSI has l_aux 0 and stores no contig names. Names for a BCF come from the header."
note "Wheat chr6D is in the BCF header, in neither index -s table, and not in the VCF CSI name list."
note "The tutorial VCF CSI depth is 8. The tutorial BCF CSI depth is 6. Readers take depth from the file."
note "small.raw.bcf is bgzip -dc, so the BCF magic is at byte 0."
note "small.stored.bcf is bcftools view -Ou. That is BGZF with stored blocks, not raw BCF."
note "bcftools view -Ou of the tutorial is about 18.7 MB and is not stored."
note "A region line identical to a stored .vcf.gz region is recorded with identical and not copied."
note "The default bcftools overlap is record, so the comparison passes --regions-overlap to match the pos files."

# Tutorial and wheat BCF, converted from the .vcf.gz this script just wrote.
bcf_run "tutorial BCF" require-sanity - \
  bcftools view --no-version -Ob \
  -o htslib/bcf/tutorial/mdp_genotype_w_depth.bcf \
  htslib/tutorial/mdp_genotype_w_depth.vcf.gz
manifest "command = bcftools view --no-version -Ob -o htslib/bcf/tutorial/mdp_genotype_w_depth.bcf htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
manifest "file = htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"

bcf_run "tutorial BCF index" require-sanity - \
  bcftools index -c htslib/bcf/tutorial/mdp_genotype_w_depth.bcf
manifest "command = bcftools index -c htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"
manifest "file = htslib/bcf/tutorial/mdp_genotype_w_depth.bcf.csi"
record_csi htslib/bcf/tutorial/mdp_genotype_w_depth.bcf.csi tutorial_bcf 14 6 0

bcf_run "tutorial BCF index -s" require-sanity htslib/bcf/tutorial/index-s.txt \
  bcftools index -s htslib/bcf/tutorial/mdp_genotype_w_depth.bcf
manifest "command = bcftools index -s htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"
manifest "file = htslib/bcf/tutorial/index-s.txt"
cmp -s htslib/bcf/tutorial/index-s.txt htslib/tutorial/index-s.txt ||
  die "tutorial BCF index -s differs from the VCF index -s"
manifest "tutorial_index_s_equal = true"

bcf_run "tutorial BCF header" require-sanity "$SCRATCH/tutorial-bcf-header.txt" \
  bcftools view --no-version -h htslib/bcf/tutorial/mdp_genotype_w_depth.bcf
tutorial_contigs=$(grep -c '^##contig=' "$SCRATCH/tutorial-bcf-header.txt" || true)
[[ "$tutorial_contigs" == "10" ]] ||
  die "tutorial BCF header has $tutorial_contigs contigs, expected 10"
manifest "tutorial_header_contigs = 10"
tutorial_vcf_aux=$(csi_l_aux htslib/tutorial/mdp_genotype_w_depth.vcf.gz.csi)
[[ "$tutorial_vcf_aux" == "49" ]] ||
  die "tutorial VCF CSI l_aux is $tutorial_vcf_aux, expected 49"
manifest "tutorial_vcf_csi_l_aux = 49"

bcf_run "wheat BCF" strict - \
  bcftools view --no-version -Ob \
  -o htslib/bcf/wheat/wheat.bcf \
  htslib/wheat/wheat.vcf.gz
manifest "command = bcftools view --no-version -Ob -o htslib/bcf/wheat/wheat.bcf htslib/wheat/wheat.vcf.gz"
manifest "file = htslib/bcf/wheat/wheat.bcf"

bcf_run "wheat BCF index" strict - \
  bcftools index -c htslib/bcf/wheat/wheat.bcf
manifest "command = bcftools index -c htslib/bcf/wheat/wheat.bcf"
manifest "file = htslib/bcf/wheat/wheat.bcf.csi"
record_csi htslib/bcf/wheat/wheat.bcf.csi wheat_bcf 14 6 0

bcf_run "wheat BCF index -s" strict htslib/bcf/wheat/index-s.txt \
  bcftools index -s htslib/bcf/wheat/wheat.bcf
manifest "command = bcftools index -s htslib/bcf/wheat/wheat.bcf"
manifest "file = htslib/bcf/wheat/index-s.txt"
cmp -s htslib/bcf/wheat/index-s.txt htslib/wheat/index-s.txt ||
  die "wheat BCF index -s differs from the VCF index -s"
manifest "wheat_index_s_equal = true"

bcf_run "wheat BCF header" strict "$SCRATCH/wheat-bcf-header.txt" \
  bcftools view --no-version -h htslib/bcf/wheat/wheat.bcf
grep -q '^##contig=<ID=chr6D,' "$SCRATCH/wheat-bcf-header.txt" ||
  die "chr6D is missing from the wheat BCF header"
wheat_header_contigs=$(grep -c '^##contig=' "$SCRATCH/wheat-bcf-header.txt" || true)
[[ "$wheat_header_contigs" == "21" ]] ||
  die "wheat BCF header has $wheat_header_contigs contigs, expected 21"
if grep -q '^chr6D'$'\t' htslib/bcf/wheat/index-s.txt; then
  die "wheat BCF index -s lists chr6D"
fi
if grep -q '^chr6D'$'\t' htslib/wheat/index-s.txt; then
  die "wheat VCF index -s lists chr6D"
fi
wheat_index_s_contigs=$(wc -l <htslib/bcf/wheat/index-s.txt | tr -d ' ')
[[ "$wheat_index_s_contigs" == "5" ]] ||
  die "wheat BCF index -s has $wheat_index_s_contigs lines, expected 5"
wheat_vcf_aux=$(csi_l_aux htslib/wheat/wheat.vcf.gz.csi)
[[ "$wheat_vcf_aux" == "58" ]] || die "wheat VCF CSI l_aux is $wheat_vcf_aux, expected 58"
wheat_vcf_names=$(csi_contig_names htslib/wheat/wheat.vcf.gz.csi)
[[ "$wheat_vcf_names" == "chr1A chr2A chr3B chr7B chr7D" ]] ||
  die "wheat VCF CSI contigs are '$wheat_vcf_names'"
manifest \
  "chr6d_in_bcf_header = true" \
  "chr6d_in_bcf_index_s = false" \
  "chr6d_in_vcf_index_s = false" \
  "chr6d_in_vcf_csi = false" \
  "header_matches_vcf_index = false" \
  "wheat_header_contigs = 21" \
  "wheat_index_s_contigs = 5" \
  "wheat_vcf_csi_l_aux = 58" \
  "wheat_vcf_csi_contigs = chr1A chr2A chr3B chr7B chr7D"

# Replay every tutorial and wheat region against the BCF. Large-pos is not converted:
# htslib refuses a position that does not fit in a BCF int32.
for entry in "${REGION_QUERIES[@]}"; do
  IFS=$'\t' read -r gz region overlap dest <<<"$entry"
  case "$gz" in
    htslib/tutorial/mdp_genotype_w_depth.vcf.gz)
      bcf=htslib/bcf/tutorial/mdp_genotype_w_depth.bcf
      stderr_mode=require-sanity
      group=tutorial
      ;;
    htslib/wheat/wheat.vcf.gz)
      bcf=htslib/bcf/wheat/wheat.bcf
      stderr_mode=strict
      group=wheat
      ;;
    *) continue ;;
  esac
  bcf_run "BCF region $region ($overlap)" "$stderr_mode" "$SCRATCH/bcf-region.txt" \
    bcftools view -H --regions-overlap "$overlap" -r "$region" "$bcf"
  manifest "command = bcftools view -H --regions-overlap $overlap -r $region $bcf"
  if cmp -s "$SCRATCH/bcf-region.txt" "$dest"; then
    manifest "identical = $dest"
  else
    out_bytes=$(wc -c <"$SCRATCH/bcf-region.txt" | tr -d ' ')
    if [[ "$out_bytes" -gt 100000 ]]; then
      die "BCF region $region ($overlap) differs from $dest ($out_bytes bytes). Not storing a second copy of a large region."
    fi
    mkdir -p htslib/bcf/regions
    base=$(basename "$dest")
    rel="htslib/bcf/regions/${group}-${base}"
    cp "$SCRATCH/bcf-region.txt" "$rel"
    manifest "file = $rel"
  fi
done
identical_count=$(awk -F ' = ' '$1 == "identical" { n++ } END { print n + 0 }' htslib/MANIFEST)
[[ "$identical_count" == "34" ]] ||
  die "expected 34 identical region lines, found $identical_count"

# htslib refuses the existing large-pos fixture. POS 3000000000 does not fit in BCF.
set +e
bcftools view --no-version -Ob -o "$SCRATCH/large-pos.bcf" \
  htslib/large-pos/large_pos.vcf.gz >"$SCRATCH/large-pos.out" 2>"$SCRATCH/large-pos.err"
large_pos_status=$?
set -e
[[ "$large_pos_status" -ne 0 ]] ||
  die "large-pos converted to BCF; htslib $HTSLIB_VERSION must refuse POS 3000000000"
rm -f "$SCRATCH/large-pos.bcf"
large_pos_line=$(grep -F "64-bit values not representable in BCF" "$SCRATCH/large-pos.err" || true)
large_pos_line=${large_pos_line%%$'\n'*}
[[ -n "$large_pos_line" ]] ||
  die "large-pos refusal did not name 64-bit values: $(cat "$SCRATCH/large-pos.err")"
manifest "command = bcftools view --no-version -Ob htslib/large-pos/large_pos.vcf.gz"
manifest "large_pos_bcf_refused = $large_pos_line"

# One small BCF, then the raw and stored-block views a reader has to tell apart.
{
  printf '%s\n' '##fileformat=VCFv4.2'
  printf '%s\n' '##contig=<ID=chr1,length=1000>'
  printf '%s\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n'
} >"$SCRATCH/small.vcf"
bcf_run "small BCF" strict - \
  bcftools view --no-version -Ob -o htslib/bcf/small/small.bcf "$SCRATCH/small.vcf"
manifest "command = bcftools view --no-version -Ob (one SNP on chr1)"
manifest "file = htslib/bcf/small/small.bcf"
bgzip -dc htslib/bcf/small/small.bcf >htslib/bcf/small/small.raw.bcf
[[ "$(prefix_bytes htslib/bcf/small/small.raw.bcf 5)" == "66 67 70 2 2" ]] ||
  die "raw BCF magic is '$(prefix_bytes htslib/bcf/small/small.raw.bcf 5)', expected 66 67 70 2 2"
manifest "command = bgzip -dc htslib/bcf/small/small.bcf"
manifest "file = htslib/bcf/small/small.raw.bcf"
bcf_run "small stored BCF" strict - \
  bcftools view --no-version -Ou -o htslib/bcf/small/small.stored.bcf "$SCRATCH/small.vcf"
[[ "$(prefix_bytes htslib/bcf/small/small.stored.bcf 4)" == "31 139 8 4" ]] ||
  die "stored BCF is not BGZF"
[[ "$(le_byte htslib/bcf/small/small.stored.bcf 18)" == "1" ]] ||
  die "stored BCF deflate block is not a stored block"
[[ "$(prefix_bytes htslib/bcf/small/small.bcf 4)" == "31 139 8 4" ]] ||
  die "compressed BCF is not BGZF"
[[ "$(le_byte htslib/bcf/small/small.bcf 18)" != "1" ]] ||
  die "compressed BCF looks like stored blocks"
bgzip -dc htslib/bcf/small/small.stored.bcf >"$SCRATCH/stored.raw"
bgzip -dc htslib/bcf/small/small.bcf >"$SCRATCH/compressed.raw"
cmp -s "$SCRATCH/stored.raw" htslib/bcf/small/small.raw.bcf ||
  die "stored-block payload is not the raw BCF"
cmp -s "$SCRATCH/compressed.raw" htslib/bcf/small/small.raw.bcf ||
  die "compressed payload is not the raw BCF"
manifest "command = bcftools view --no-version -Ou (one SNP on chr1)"
manifest "file = htslib/bcf/small/small.stored.bcf"
bcf_run "small BCF text" strict htslib/bcf/small/small.view.vcf \
  bcftools view --no-version htslib/bcf/small/small.bcf
manifest "command = bcftools view --no-version htslib/bcf/small/small.bcf"
manifest "file = htslib/bcf/small/small.view.vcf"

# The largest 1-based POS that fits in a BCF int32. 2147483648 is refused, as is large-pos.
{
  printf '%s\n' '##fileformat=VCFv4.2'
  printf '%s\n' '##contig=<ID=chr1,length=2147483647>'
  printf '%s\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr1\t2147483647\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n'
} >"$SCRATCH/boundary.vcf"
bcf_run "boundary BCF" strict - \
  bcftools view --no-version -Ob -o htslib/bcf/boundary/pos2147483647.bcf "$SCRATCH/boundary.vcf"
manifest "command = bcftools view --no-version -Ob (one record at POS 2147483647)"
manifest "file = htslib/bcf/boundary/pos2147483647.bcf"
bcf_run "boundary BCF index" strict - \
  bcftools index -c htslib/bcf/boundary/pos2147483647.bcf
manifest "command = bcftools index -c htslib/bcf/boundary/pos2147483647.bcf"
manifest "file = htslib/bcf/boundary/pos2147483647.bcf.csi"
record_csi htslib/bcf/boundary/pos2147483647.bcf.csi boundary 14 6 0
bcf_run "boundary BCF text" strict htslib/bcf/boundary/pos2147483647.view.vcf \
  bcftools view --no-version htslib/bcf/boundary/pos2147483647.bcf
grep -q $'^chr1\t2147483647\t' htslib/bcf/boundary/pos2147483647.view.vcf ||
  die "boundary view text does not contain POS 2147483647"
manifest "command = bcftools view --no-version htslib/bcf/boundary/pos2147483647.bcf"
manifest "file = htslib/bcf/boundary/pos2147483647.view.vcf"

# Edge-case records share one file. A no-sample record cannot: bcftools rejects the column count.
awk '
  BEGIN {
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=chr1,length=1000>"
    print "##FILTER=<ID=q10,Description=\"qual\">"
    print "##INFO=<ID=I8,Number=1,Type=Integer,Description=\"int8\">"
    print "##INFO=<ID=I16,Number=1,Type=Integer,Description=\"int16\">"
    print "##INFO=<ID=I32,Number=1,Type=Integer,Description=\"int32\">"
    print "##INFO=<ID=FL,Number=1,Type=Float,Description=\"float\">"
    for (i = 0; i < 130; i++)
      printf "##INFO=<ID=K%03d,Number=1,Type=Integer,Description=\"dict %d\">\n", i, i
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=F8,Number=1,Type=Integer,Description=\"int8\">"
    print "##FORMAT=<ID=F16,Number=1,Type=Integer,Description=\"int16\">"
    print "##FORMAT=<ID=F32,Number=1,Type=Integer,Description=\"int32\">"
    print "##FORMAT=<ID=FF,Number=1,Type=Float,Description=\"float\">"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    printf "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0\t1\n"
    printf "chr1\t20\tLONG_VARIANT_ID_16\tA\tC\t.\t.\tI8=100;I16=200;I32=40000;FL=1.5;K129=7\tGT:F8:F16:F32:FF\t0:100:200:40000:1.5\t0/1:101:201:40001:2.5\n"
    alts = ""
    for (i = 0; i < 64; i++) {
      allele = "A" substr("ACGT", i % 4 + 1, 1) i
      alts = alts (i ? "," : "") allele
    }
    printf "chr1\t30\t.\tA\t%s\t.\tPASS\t.\tGT\t10/0\t64/0\n", alts
  }
' >"$SCRATCH/edge.vcf"
bcf_run "edge BCF" strict - \
  bcftools view --no-version -Ob -o htslib/bcf/edge/edge.bcf "$SCRATCH/edge.vcf"
bgzip -dc htslib/bcf/edge/edge.bcf >"$SCRATCH/edge.raw"
python3 - "$SCRATCH/edge.raw" <<'PY'
import struct
import sys

path = sys.argv[1]
raw = open(path, "rb").read()
if raw[:5] != b"BCF\x02\x02":
    sys.exit(f"edge magic is {raw[:5]!r}")
l_text = struct.unpack_from("<I", raw, 5)[0]
header = raw[9:9 + l_text].split(b"\0", 1)[0].decode("ascii")
n_dict = sum(
    1
    for line in header.splitlines()
    if line.startswith(("##INFO=", "##FORMAT=", "##FILTER="))
)
if n_dict <= 127:
    sys.exit(f"dictionary has {n_dict} entries, want more than 127")


def read_typed(buf, i, repeat=1):
    type_byte = buf[i]
    i += 1
    typ = type_byte & 0xF
    size = type_byte >> 4
    if size == 15:
        size_byte = buf[i]
        i += 1
        size_type = size_byte & 0xF
        if size_byte >> 4 != 1:
            sys.exit(f"size descriptor {size_byte:#x}")
        if size_type == 1:
            size = struct.unpack_from("<b", buf, i)[0]
            i += 1
        elif size_type == 2:
            size = struct.unpack_from("<h", buf, i)[0]
            i += 2
        elif size_type == 3:
            size = struct.unpack_from("<i", buf, i)[0]
            i += 4
        else:
            sys.exit(f"bad size type {size_type}")
    width = {0: 0, 1: 1, 2: 2, 3: 4, 5: 4, 7: 1}.get(typ)
    if width is None:
        sys.exit(f"bad type {typ} (byte {type_byte:#x})")
    if size < 0:
        sys.exit(f"negative size {size}")
    nbytes = width * size * repeat
    data = buf[i:i + nbytes]
    if len(data) != nbytes:
        sys.exit("typed value extends past the buffer")
    return type_byte, typ, size, data, i + nbytes


def ints(typ, data):
    fmt, width = {1: ("<b", 1), 2: ("<h", 2), 3: ("<i", 4)}[typ]
    return [struct.unpack_from(fmt, data, off)[0] for off in range(0, len(data), width)]


def parse_record(shared, indiv):
    pos = struct.unpack_from("<i", shared, 4)[0]
    qual = struct.unpack_from("<I", shared, 12)[0]
    n_allele_info, n_fmt_sample = struct.unpack_from("<II", shared, 16)
    n_info = n_allele_info & 0xFFFF
    n_allele = n_allele_info >> 16
    n_fmt = n_fmt_sample >> 24
    n_sample = n_fmt_sample & 0xFFFFFF
    i = 24
    id_byte, id_typ, id_size, id_data, i = read_typed(shared, i)
    for _ in range(n_allele):
        _, _, _, _, i = read_typed(shared, i)
    filter_byte, filter_typ, filter_size, filter_data, i = read_typed(shared, i)
    infos = []
    for _ in range(n_info):
        key = read_typed(shared, i)
        i = key[4]
        val = read_typed(shared, i)
        i = val[4]
        infos.append((key, val))
    if i != len(shared):
        sys.exit(f"POS {pos + 1} shared has {len(shared) - i} leftover bytes")
    fmts = []
    ii = 0
    for _ in range(n_fmt):
        key = read_typed(indiv, ii)
        ii = key[4]
        val = read_typed(indiv, ii, n_sample)
        ii = val[4]
        fmts.append((key, val))
    if ii != len(indiv):
        sys.exit(f"POS {pos + 1} indiv has {len(indiv) - ii} leftover bytes")
    return {
        "pos": pos,
        "qual": qual,
        "n_allele": n_allele,
        "n_sample": n_sample,
        "id_byte": id_byte,
        "id_data": id_data,
        "filter_typ": filter_typ,
        "filter_size": filter_size,
        "infos": infos,
        "fmts": fmts,
    }


off = 9 + l_text
records = []
while off < len(raw):
    l_shared, l_indiv = struct.unpack_from("<II", raw, off)
    start = off + 8
    end = start + l_shared + l_indiv
    if end > len(raw):
        sys.exit("record extends past the raw BCF")
    records.append(parse_record(raw[start:start + l_shared], raw[start + l_shared:end]))
    off = end
if off != len(raw) or len(records) != 3:
    sys.exit(f"{len(records)} records, consumed {off} of {len(raw)}")

hap, mixed, wide = records
if hap["pos"] != 9 or hap["n_sample"] != 2:
    sys.exit(f"haploid record pos {hap['pos'] + 1} samples {hap['n_sample']}")
gt_key, gt_val = hap["fmts"][0]
if gt_key[1] != 2 or ints(gt_key[1], gt_key[3]) != [0x88]:
    sys.exit(f"haploid GT key {gt_key[1]} {gt_key[3].hex()}")
if gt_val[1] != 1 or gt_val[2] != 1 or gt_val[3] != bytes((0x03, 0x05)):
    sys.exit(f"haploid GT {gt_val[1]} {gt_val[2]} {gt_val[3].hex()}")

if mixed["pos"] != 19 or mixed["qual"] != 0x7F800001:
    sys.exit(f"mixed record pos {mixed['pos'] + 1} qual {mixed['qual']:#x}")
if mixed["id_byte"] != 0xF7 or mixed["id_data"] != b"LONG_VARIANT_ID_16":
    sys.exit(f"ID byte {mixed['id_byte']:#x} data {mixed['id_data']!r}")
if mixed["filter_typ"] != 0 or mixed["filter_size"] != 0:
    sys.exit(f"missing FILTER is type {mixed['filter_typ']} size {mixed['filter_size']}")
expect_info = [
    (1, 2, 1, [100]),
    (1, 3, 2, [200]),
    (1, 4, 3, [40000]),
    (1, 5, 5, None),
    (2, 0x87, 1, [7]),
]
if len(mixed["infos"]) != len(expect_info):
    sys.exit(f"{len(mixed['infos'])} INFO fields")
for (key, val), (key_type, key_value, val_type, val_ints) in zip(mixed["infos"], expect_info):
    if key[1] != key_type or ints(key[1], key[3]) != [key_value]:
        sys.exit(f"INFO key type {key[1]} bytes {key[3].hex()}")
    if val[1] != val_type:
        sys.exit(f"INFO value type {val[1]}, expected {val_type}")
    if val_ints is not None and ints(val[1], val[3]) != val_ints:
        sys.exit(f"INFO value {val[3].hex()}")
    if val_type == 5 and struct.unpack("<I", val[3])[0] != 0x3FC00000:
        sys.exit(f"float INFO {val[3].hex()}")
gt_key, gt_val = mixed["fmts"][0]
if gt_val[1] != 1 or gt_val[2] != 2 or gt_val[3] != bytes((0x03, 0x81, 0x02, 0x04)):
    sys.exit(f"mixed GT {gt_val[1]} {gt_val[2]} {gt_val[3].hex()}")
fmt_expect = [
    (1, [100, 101]),
    (2, [200, 201]),
    (3, [40000, 40001]),
    (5, None),
]
for (key, val), (val_type, val_ints) in zip(mixed["fmts"][1:], fmt_expect):
    if key[1] != 2:
        sys.exit(f"FORMAT key type {key[1]}, expected int16")
    if val[1] != val_type:
        sys.exit(f"FORMAT value type {val[1]}, expected {val_type}")
    if val_ints is not None and ints(val[1], val[3]) != val_ints:
        sys.exit(f"FORMAT value {val[3].hex()}")
float_fmt = mixed["fmts"][4][1][3]
if struct.unpack("<ff", float_fmt) != (1.5, 2.5):
    sys.exit(f"FORMAT float {float_fmt.hex()}")

if wide["pos"] != 29 or wide["n_allele"] != 65:
    sys.exit(f"wide record pos {wide['pos'] + 1} alleles {wide['n_allele']}")
gt_key, gt_val = wide["fmts"][0]
if gt_val[1] != 2 or gt_val[2] != 2 or ints(2, gt_val[3]) != [22, 2, 130, 2]:
    sys.exit(f"wide GT type {gt_val[1]} {gt_val[3].hex()}")
PY
manifest "command = bcftools view --no-version -Ob (edge-case VCF: haploid, mixed ploidy, 64 ALT alleles, IDX past 127)"
manifest "file = htslib/bcf/edge/edge.bcf"
bcf_run "edge BCF text" strict htslib/bcf/edge/edge.view.vcf \
  bcftools view --no-version htslib/bcf/edge/edge.bcf
grep -q $'LONG_VARIANT_ID_16' htslib/bcf/edge/edge.view.vcf ||
  die "edge view text is missing the long ID"
grep -q $'\t10/0\t' htslib/bcf/edge/edge.view.vcf ||
  die "edge view text is missing allele index 10"
manifest "command = bcftools view --no-version htslib/bcf/edge/edge.bcf"
manifest "file = htslib/bcf/edge/edge.view.vcf"

{
  printf '%s\n' '##fileformat=VCFv4.3'
  printf '%s\n' '##contig=<ID=chr1,length=1000>'
  printf '%s\n' '##INFO=<ID=NS,Number=1,Type=Integer,Description="samples">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr1\t1\t.\tA\tC\t.\t.\tNS=2\n'
} >"$SCRATCH/no-sample.vcf"
bcf_run "no-sample BCF" strict - \
  bcftools view --no-version -Ob -o htslib/bcf/edge/no-sample.bcf "$SCRATCH/no-sample.vcf"
bcf_run "no-sample header" strict "$SCRATCH/no-sample-header.txt" \
  bcftools view --no-version -h htslib/bcf/edge/no-sample.bcf
no_sample_cols=$(awk -F '\t' '/^#CHROM/ { print NF }' "$SCRATCH/no-sample-header.txt")
[[ "$no_sample_cols" == "8" ]] ||
  die "no-sample #CHROM line has $no_sample_cols columns, expected 8"
manifest "command = bcftools view --no-version -Ob (one record, no samples)"
manifest "file = htslib/bcf/edge/no-sample.bcf"
bcf_run "no-sample BCF text" strict htslib/bcf/edge/no-sample.view.vcf \
  bcftools view --no-version htslib/bcf/edge/no-sample.bcf
manifest "command = bcftools view --no-version htslib/bcf/edge/no-sample.bcf"
manifest "file = htslib/bcf/edge/no-sample.view.vcf"

# Corrupt files are valid BGZF. The damage is in the BCF payload, except no-eof,
# which drops the trailing 28-byte block the way variants/no-eof.vcf.gz does.
{
  printf '\000'
  dd bs=1 skip=1 <htslib/bcf/small/small.raw.bcf 2>/dev/null
} | bgzip -c >htslib/bcf/corrupt/bad-magic.bcf
bgzip -dc htslib/bcf/corrupt/bad-magic.bcf >"$SCRATCH/bad-magic.raw"
[[ "$(prefix_bytes "$SCRATCH/bad-magic.raw" 5)" == "0 67 70 2 2" ]] ||
  die "bad magic inflates to '$(prefix_bytes "$SCRATCH/bad-magic.raw" 5)'"
manifest "command = bgzip -dc htslib/bcf/small/small.bcf"
manifest "command = set byte 0 from 0x42 to 0x00"
manifest "command = bgzip -c"
manifest "file = htslib/bcf/corrupt/bad-magic.bcf"

{
  dd bs=1 count=4 <htslib/bcf/small/small.raw.bcf 2>/dev/null
  printf '\001'
  dd bs=1 skip=5 <htslib/bcf/small/small.raw.bcf 2>/dev/null
} | bgzip -c >htslib/bcf/corrupt/bcf21.bcf
bgzip -dc htslib/bcf/corrupt/bcf21.bcf >"$SCRATCH/bcf21.raw"
[[ "$(prefix_bytes "$SCRATCH/bcf21.raw" 5)" == "66 67 70 2 1" ]] ||
  die "BCF 2.1 magic inflates to '$(prefix_bytes "$SCRATCH/bcf21.raw" 5)'"
manifest "command = bgzip -dc htslib/bcf/small/small.bcf"
manifest "command = set byte 4 from 0x02 to 0x01"
manifest "command = bgzip -c"
manifest "file = htslib/bcf/corrupt/bcf21.bcf"

raw_l_text=$(le_u32_at htslib/bcf/small/small.raw.bcf 5)
keep=$((9 + raw_l_text + 12))
dd bs=1 count="$keep" <htslib/bcf/small/small.raw.bcf 2>/dev/null |
  bgzip -c >htslib/bcf/corrupt/truncated-record.bcf
bgzip -dc htslib/bcf/corrupt/truncated-record.bcf >"$SCRATCH/truncated.raw"
truncated_bytes=$(wc -c <"$SCRATCH/truncated.raw" | tr -d ' ')
[[ "$truncated_bytes" == "$keep" ]] ||
  die "truncated record inflates to $truncated_bytes bytes, expected $keep"
set +e
bcftools view --no-version htslib/bcf/corrupt/truncated-record.bcf \
  >"$SCRATCH/truncated.out" 2>"$SCRATCH/truncated.err"
truncated_status=$?
set -e
[[ "$truncated_status" -ne 0 ]] || die "truncated BCF record was readable"
manifest "command = bgzip -dc htslib/bcf/small/small.bcf"
manifest "command = keep the header, the 8-byte record lengths, and 4 bytes of the record"
manifest "command = bgzip -c"
manifest "file = htslib/bcf/corrupt/truncated-record.bcf"

small_walk=$(walk_bgzf htslib/bcf/small/small.bcf) || die "small.bcf is not BGZF"
[[ "${small_walk##* }" == "1" ]] || die "small.bcf has ${small_walk##* } EOF blocks, expected 1"
small_bytes=$(wc -c <htslib/bcf/small/small.bcf | tr -d ' ')
dd bs=1 count=$((small_bytes - 28)) <htslib/bcf/small/small.bcf \
  >htslib/bcf/corrupt/no-eof.bcf 2>/dev/null
no_eof_walk=$(walk_bgzf htslib/bcf/corrupt/no-eof.bcf) || die "no-eof BCF is not a sequence of BGZF blocks"
[[ "${no_eof_walk##* }" == "0" ]] || die "no-eof BCF still has an EOF block"
manifest "command = bcftools view --no-version -Ob (one SNP on chr1)"
manifest "command = drop the final 28-byte BGZF EOF block"
manifest "file = htslib/bcf/corrupt/no-eof.bcf"

# --- BGZF layout ---------------------------------------------------------------

manifest "" "[bgzf-layout]"
note "bgzip -c of a VCF line and a VCF header that do not fit in one 65280-byte BGZF block."
note "The plain VCFs are not stored. long_line_bytes is the data line without its newline."
note "long_header_bytes is the text before #CHROM."

awk '
  BEGIN {
    header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    prefix = "chr1\t100\t.\tA\t"
    suffix = "\t.\tPASS\t.\tGT\ts1"
    # length() of the record excludes the newline that terminates it.
    need = 65537 - length(prefix) - length(suffix)
    if (need <= 0) exit 2
    printf "%s%s", header, prefix
    for (i = 0; i < need; i++) printf "C"
    printf "%s\n", suffix
  }
' >"$SCRATCH/long-line.vcf" || die "could not build the long VCF line"
long_line_bytes=$(awk '/^chr1/ { print length($0); exit }' "$SCRATCH/long-line.vcf")
[[ "$long_line_bytes" == "65537" ]] ||
  die "long VCF line is $long_line_bytes bytes, expected 65537"
bgzip -c "$SCRATCH/long-line.vcf" >htslib/bgzf-layout/long-line.vcf.gz
long_line_blocks=$(bgzf_data_blocks htslib/bgzf-layout/long-line.vcf.gz) ||
  die "could not walk long-line BGZF blocks"
[[ "$long_line_blocks" -ge 2 ]] ||
  die "long VCF line fits in $long_line_blocks BGZF data blocks"
manifest "command = bgzip -c (VCF whose data line is 65537 bytes)"
manifest "long_line_bytes = $long_line_bytes"
manifest "long_line_data_blocks = $long_line_blocks"
manifest "file = htslib/bgzf-layout/long-line.vcf.gz"

awk '
  BEGIN {
    header = "##fileformat=VCFv4.2\n"
    i = 0
    while (length(header) <= 65536) {
      i++
      header = header sprintf("##INFO=<ID=P%04d,Number=0,Type=Flag,Description=\"pad\">\n", i)
    }
    header = header "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    printf "%s", header
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    printf "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n"
  }
' >"$SCRATCH/long-header.vcf" || die "could not build the long VCF header"
long_header_bytes=$(awk '
  /^#CHROM/ { print n; exit }
  { n += length($0) + 1 }
' "$SCRATCH/long-header.vcf")
[[ "$long_header_bytes" -gt 65536 ]] ||
  die "VCF header is $long_header_bytes bytes, expected more than 65536"
bgzip -c "$SCRATCH/long-header.vcf" >htslib/bgzf-layout/long-header.vcf.gz
long_header_blocks=$(bgzf_data_blocks htslib/bgzf-layout/long-header.vcf.gz) ||
  die "could not walk long-header BGZF blocks"
[[ "$long_header_blocks" -ge 2 ]] ||
  die "long VCF header fits in $long_header_blocks BGZF data blocks"
manifest "command = bgzip -c (VCF whose header is longer than 65536 bytes)"
manifest "long_header_bytes = $long_header_bytes"
manifest "long_header_data_blocks = $long_header_blocks"
manifest "file = htslib/bgzf-layout/long-header.vcf.gz"

# --- BCF encodings of the TASSEL 5 VCF-export goldens (tassel-6#129) ----------
# A golden is TASSEL 5 text. These files are htslib's BCF framing of that text.
# The .vcf.gz and its CSI stay in scratch: the plain VCF has no ##contig lines,
# and bcftools view -Ob refuses it until it is bgzipped and indexed. Level-0
# binaries are not stored. The GT-only tutorial at level 0 is about 1.9 MB.

mkdir -p htslib/bcf/golden

manifest "" "[golden-bcf]"
note "htslib $HTSLIB_VERSION BCF encodings of the TASSEL 5 VCF-export goldens (tassel-6#119)."
note "A golden is TASSEL 5 text. These files are htslib's framing of that text. The two oracles stay distinct (ADR-0022)."
note "Editing these files does not use the golden-update label, and testdata/golden/ is not modified."
note "Each golden has no ##contig lines, so bcftools view -Ob refuses the plain VCF."
note "The route is bgzip, bcftools index -c, then bcftools view --no-version -Ob of the bgzipped file. The .vcf.gz and .csi stay in scratch."
note "blocks.txt is one decimal BSIZE+1 per BGZF member, including the 28-byte EOF member, from the same walk as walk_bgzf."
note "Level-0 binaries are not stored. The full tutorial at level 0 is about 18.7 MB and the GT-only tutorial at level 0 is about 1.9 MB."
note "smallest_golden is the golden with the fewest bytes. level0_blocks_match says whether -Ob and -Ob -l 0 share a member-size list for that golden."
note "When they differ, the level-0 list is committed as text and the level-0 binary is deleted. #130's stored writer matches that list when the file exists, and the -Ob list otherwise."
note "Inflated bytes always come from the committed -Ob BCF."

expected_goldens=(
  golden/vcf-export/mdp_genotype.vcf
  golden/vcf-export/include-list.vcf
  golden/vcf-export/vcf_export_edges.vcf
)
golden_vcfs=()
in_vcf_export=0
while IFS= read -r line || [[ -n "$line" ]]; do
  case "$line" in
    "[vcf-export]")
      in_vcf_export=1
      continue
      ;;
    \[*\])
      if [[ "$in_vcf_export" == 1 ]]; then
        break
      fi
      ;;
  esac
  if [[ "$in_vcf_export" == 1 ]]; then
    case "$line" in
      golden/vcf-export/*.vcf" = exact")
        golden_vcfs+=("${line% = exact}")
        ;;
    esac
  fi
done < golden/MANIFEST

[[ "${#golden_vcfs[@]}" -eq "${#expected_goldens[@]}" ]] ||
  die "golden MANIFEST [vcf-export] has ${#golden_vcfs[@]} VCF files, expected ${#expected_goldens[@]}"
for i in 0 1 2; do
  [[ "${golden_vcfs[$i]}" == "${expected_goldens[$i]}" ]] ||
    die "vcf-export golden $((i + 1)) is ${golden_vcfs[$i]}, expected ${expected_goldens[$i]}"
done

smallest_name=
smallest_bytes=
for golden in "${golden_vcfs[@]}"; do
  [[ -f "$golden" ]] || die "missing VCF-export golden $golden"
  name=$(basename "$golden" .vcf)
  bytes=$(wc -c <"$golden" | tr -d ' ')
  if [[ -n "$smallest_bytes" && "$bytes" == "$smallest_bytes" ]]; then
    die "goldens $smallest_name and $name are both $bytes bytes"
  fi
  if [[ -z "$smallest_bytes" || "$bytes" -lt "$smallest_bytes" ]]; then
    smallest_name=$name
    smallest_bytes=$bytes
  fi
done
manifest "smallest_golden = $smallest_name"

for golden in "${golden_vcfs[@]}"; do
  name=$(basename "$golden" .vcf)
  gz="$SCRATCH/golden-${name}.vcf.gz"
  bcf="htslib/bcf/golden/${name}.bcf"
  blocks="htslib/bcf/golden/${name}.blocks.txt"
  bgzip -c "$golden" >"$gz"
  # index -c writes no PL/GQ warning. view does, because the golden header declares PL and GQ.
  bcf_run "$name golden index" strict - \
    bcftools index -c "$gz"
  bcf_run "$name golden BCF" require-sanity - \
    bcftools view --no-version -Ob -o "$bcf" "$gz"
  bgzf_member_sizes "$bcf" >"$blocks" || die "could not walk BGZF members of $bcf"
  walk=$(walk_bgzf "$bcf") || die "$bcf is not a sequence of BGZF blocks"
  [[ "${walk##* }" == "1" ]] || die "$bcf has ${walk##* } EOF blocks, expected 1"
  [[ "$(tail -n 1 "$blocks")" == "28" ]] || die "$blocks does not end with the 28-byte EOF member"
  sum=$(awk '{ s += $1 } END { print s + 0 }' "$blocks")
  bcf_bytes=$(wc -c <"$bcf" | tr -d ' ')
  [[ "$sum" == "$bcf_bytes" ]] || die "$blocks sums to $sum, but $bcf is $bcf_bytes bytes"
  digest=$(sha256 "$golden")
  manifest \
    "command = bgzip -c $golden" \
    "command = bcftools index -c" \
    "command = bcftools view --no-version -Ob -o $bcf" \
    "file = $bcf" \
    "command = BGZF member sizes of $bcf, one decimal BSIZE+1 per member, including the EOF block" \
    "file = $blocks" \
    "${name}_golden = $golden" \
    "${name}_sha256 = $digest"

  if [[ "$name" == "$smallest_name" ]]; then
    l0="$SCRATCH/golden-${name}.l0.bcf"
    l0_blocks="$SCRATCH/golden-${name}.blocks-l0.txt"
    bcf_run "$name golden BCF level 0" require-sanity - \
      bcftools view --no-version -Ob -l 0 -o "$l0" "$gz"
    bgzf_member_sizes "$l0" >"$l0_blocks" || die "could not walk BGZF members of the level-0 $name BCF"
    if cmp -s "$blocks" "$l0_blocks"; then
      manifest "level0_blocks_match = true"
    else
      cp "$l0_blocks" "htslib/bcf/golden/${name}.blocks-l0.txt"
      manifest \
        "command = bcftools view --no-version -Ob -l 0" \
        "command = BGZF member sizes of the level-0 file; the level-0 binary is not stored" \
        "file = htslib/bcf/golden/${name}.blocks-l0.txt" \
        "level0_blocks_match = false"
    fi
    rm -f "$l0"
    [[ ! -e "$l0" ]] || die "level-0 BCF $l0 was not deleted"
  fi
done

# --- Every file is listed, and the tree stays under 10 MB ---------------------

leftover=$(find htslib -type f ! -name MANIFEST | sort | while read -r f; do
  grep -qF "file = $f" htslib/MANIFEST || printf '%s\n' "$f"
done)
[[ -z "$leftover" ]] || die "unlisted files in testdata/htslib: $leftover"

missing=$(awk -F ' = ' '$1 == "file" { print $2 }' htslib/MANIFEST | while read -r f; do
  [[ -f "$f" ]] || printf '%s\n' "$f"
done)
[[ -z "$missing" ]] || die "MANIFEST lists missing files: $missing"

bytes=$(find htslib -type f -exec wc -c {} + | awk '/total$/ { print $1 }')
[[ -n "$bytes" && "$bytes" -lt 10485760 ]] || die "testdata/htslib is ${bytes:-unknown} bytes, over 10 MB"
# tassel-6#118 keeps this tree at or below 9.5 MiB so later fixtures still fit under
# the 10 MiB refusal above. tassel-6#129 added the golden BCFs and stayed under
# 9961472, so this gate remains.
[[ "$bytes" -le 9961472 ]] ||
  die "testdata/htslib is $bytes bytes, over 9.5 MiB (9961472)."

file_count=$(find htslib -type f | wc -l | tr -d ' ')
echo "Wrote $file_count files to testdata/htslib ($(du -sh "$HTSLIB" | awk '{ print $1 }'), $bytes bytes)"
