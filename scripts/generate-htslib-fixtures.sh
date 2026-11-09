#!/usr/bin/env bash
# Generate the htslib 1.25 oracle in testdata/htslib/ (ADR-0021, tassel-6#80).
#
# Rebuilds testdata/htslib/ and its MANIFEST from scratch. Output is deterministic:
# running the script twice gives byte-identical files. Nothing under testdata/golden/
# or testdata/tutorial/ is modified.
#
# Requires bgzip, tabix, and bcftools built on htslib 1.25, Apple gzip 487.0.1
# (the plain-gzip fixtures are pinned to that binary), and xxd, od, dd, and awk.
# This is the maintainer's Mac (`brew install htslib bcftools`). Do not regenerate
# the plain-gzip fixtures with GNU gzip.
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

# Data blocks of a BGZF file that ends with exactly one EOF block.
bgzf_data_blocks() {
  local got members eofs
  got=$(walk_bgzf "$1") || return 1
  members=${got%% *}
  eofs=${got##* }
  [[ "$eofs" == "1" ]] || return 1
  echo $((members - 1))
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
  done
}

# --- Toolchain ----------------------------------------------------------------

export LC_ALL=C

command -v bgzip >/dev/null || die "bgzip not found; brew install htslib"
command -v tabix >/dev/null || die "tabix not found; brew install htslib"
command -v bcftools >/dev/null || die "bcftools not found; brew install bcftools"
command -v gzip >/dev/null || die "gzip not found"
command -v xxd >/dev/null || die "xxd not found"
command -v od >/dev/null || die "od not found"

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
  "# htslib $HTSLIB_VERSION fixtures for compressed and indexed VCF (ADR-0021, tassel-6#80)." \
  "# Generated by scripts/generate-htslib-fixtures.sh. Do not edit by hand." \
  "#" \
  "# Paths are relative to testdata/, and every command runs from there." \
  "#   file     may repeat; each value is one output file" \
  "#   command  may repeat; the invocation that produced the following file" \
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

file_count=$(find htslib -type f | wc -l | tr -d ' ')
echo "Wrote $file_count files to testdata/htslib ($(du -sh "$HTSLIB" | awk '{ print $1 }'), $bytes bytes)"
