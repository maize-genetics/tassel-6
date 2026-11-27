#!/usr/bin/env bash
# Generate the TASSEL 5 golden outputs in testdata/golden/ (ADR-0009).
#
# Downloads the pinned TASSEL 5 standalone release and the TASSEL tutorial data into
# build/goldens/ (git-ignored), copies the tutorial inputs into testdata/tutorial/, writes
# the synthetic vcf_export_edges HapMap there, then rebuilds testdata/golden/ and its
# MANIFEST from scratch. Output is deterministic: running the script twice gives
# byte-identical files.
#
# Requires curl, tar, and a JDK 21 or newer (java and javac from JAVA_HOME or PATH).
#
# Usage: scripts/generate-goldens.sh
set -euo pipefail

TASSEL5_VERSION="5.2.98"
TASSEL5_TAG="v5.2.98"
TASSEL5_URL="https://github.com/maize-genetics/tassel/releases/download/v5.2.98/tassel-5-standalone-v5.2.98.tar.gz"
TASSEL5_SHA256="52e7675aec4c7ab34bb6e0cfe94b26709f8fd8dd761b4d7892b3aa0e38df3c6e"
TUTORIAL_URL="https://github.com/maize-genetics/tassel_test_data/releases/download/v1.0.0/tassel_test_data_v1.tar.gz"
TUTORIAL_SHA256="6a0c515f7c3e6053f6c4cc31705da3e2e655a055e8bf927859fd2ee2fb757424"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CACHE="$ROOT/build/goldens"
TESTDATA="$ROOT/testdata"
GOLDEN="$TESTDATA/golden"
MANIFEST="$GOLDEN/MANIFEST"
LOG="$CACHE/tassel5.log"

die() {
  echo "generate-goldens: $*" >&2
  exit 1
}

sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

# fetch URL SHA256 DEST: download once, verify every time.
fetch() {
  if [[ ! -f "$3" ]]; then
    echo "Downloading $1"
    curl -fsSL -o "$3.part" "$1"
    mv "$3.part" "$3"
  fi
  [[ "$(sha256 "$3")" == "$2" ]] || die "checksum mismatch for $3; delete it and rerun"
}

# --- Toolchain ----------------------------------------------------------------

if [[ -n "${JAVA_HOME:-}" ]]; then
  JAVA="$JAVA_HOME/bin/java"
  JAVAC="$JAVA_HOME/bin/javac"
else
  JAVA="java"
  JAVAC="javac"
fi
command -v "$JAVA" >/dev/null || die "java not found; set JAVA_HOME or put a JDK 21+ on PATH"
command -v "$JAVAC" >/dev/null || die "javac not found; a full JDK 21+ is required"
JAVA_MAJOR="$("$JAVA" -XshowSettings:properties -version 2>&1 |
  awk -F' = ' '/java\.specification\.version/ { print $2 }')"
[[ "$JAVA_MAJOR" =~ ^[0-9]+$ ]] && ((JAVA_MAJOR >= 21)) ||
  die "TASSEL $TASSEL5_VERSION needs Java 21 or newer; found '$JAVA_MAJOR'"

# --- Pinned TASSEL 5 and tutorial data ----------------------------------------

mkdir -p "$CACHE"
fetch "$TASSEL5_URL" "$TASSEL5_SHA256" "$CACHE/tassel5.tar.gz"
fetch "$TUTORIAL_URL" "$TUTORIAL_SHA256" "$CACHE/tutorial.tar.gz"

T5="$CACHE/tassel-5-standalone-$TASSEL5_TAG"
[[ -f "$T5/sTASSEL.jar" ]] || tar -xzf "$CACHE/tassel5.tar.gz" -C "$CACHE"
[[ -d "$CACHE/dataFiles/Tutorial" ]] || tar -xzf "$CACHE/tutorial.tar.gz" -C "$CACHE"

FULL_CLASSES="$CACHE/full-precision-classes"
rm -rf "$FULL_CLASSES"
"$JAVAC" --release 21 -d "$FULL_CLASSES" "$ROOT/scripts/golden/DoubleFormat.java"

# Colt 1.2.0 ships in the standalone's lib/; Distributions.java calls it directly.
COLT="$T5/lib/colt-1.2.0.jar"
[[ -f "$COLT" ]] || die "colt-1.2.0.jar not found in $T5/lib"
COMMONS_MATH="$T5/lib/commons-math3-3.4.1.jar"
[[ -f "$COMMONS_MATH" ]] || die "commons-math3-3.4.1.jar not found in $T5/lib"
DIST_CLASSES="$CACHE/distributions-classes"
rm -rf "$DIST_CLASSES"
"$JAVAC" --release 21 -cp "$COLT:$COMMONS_MATH" -d "$DIST_CLASSES" "$ROOT/scripts/golden/Distributions.java"

INPUTS=(mdp_genotype.hmp.txt mdp_genotype_w_depth.vcf mdp_traits.txt mdp_population_structure.txt
  mdp_genotype_chr9_10.hmp.txt)
mkdir -p "$TESTDATA/tutorial"
for f in "${INPUTS[@]}"; do
  cp "$CACHE/dataFiles/Tutorial/$f" "$TESTDATA/tutorial/$f"
done

# TASSEL 5's own phenotype test inputs, from the same tarball, as dataFiles-relative paths.
PHENOTYPE_INPUTS=(
  PhenotypeTests/basic_phenotype_test.txt CandidateTests/earht_with_rep.txt
  PhenotypeTests/version4format_test.txt PhenotypeTests/version4format_covariate.txt
  CandidateTests/mdp_traits_nomissing.txt
)
mkdir -p "$TESTDATA/phenotype"
for f in "${PHENOTYPE_INPUTS[@]}"; do
  cp "$CACHE/dataFiles/$f" "$TESTDATA/phenotype/$(basename "$f")"
done

# Synthetic HapMap for the VCF-export edge cases. This file is an input, not a golden;
# the VCF TASSEL 5 writes from it is the golden. Two characters per call, and TASSEL 5
# reads the second character as the first allele, so CA is genotype A/C. Sites stay in
# position order, and counts are unequal so the major allele is REF with no tie.
#   plus_minus    +/-    ++ ++ --    + and - listed and called
#   n_in_alleles  A/N/C  AA AA CC    N token stripped from the alleles column
#   triallelic    A/C/G  AA CA GG    three called alleles
#   no_calls      N      NN NN NN    no called allele, so REF is . (NA has length 2, so it would be kept)
#   unlisted      A/C    AA AA GG    G is called but not listed
#   n_prefix      +/A    ++ ++ --    unlisted - is written N; uncalled A is written NA
#   known_indel   AT/A   TT TT --    length>1 known variant (indelInKnownVariant)
printf '%s\n' \
  $'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\tt1\tt2\tt3' \
  $'plus_minus\t+/-\t1\t10\t+\tNA\tNA\tNA\tNA\tNA\tNA\t++\t++\t--' \
  $'n_in_alleles\tA/N/C\t1\t20\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tAA\tCC' \
  $'triallelic\tA/C/G\t1\t30\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tCA\tGG' \
  $'no_calls\tN\t1\t40\t+\tNA\tNA\tNA\tNA\tNA\tNA\tNN\tNN\tNN' \
  $'unlisted\tA/C\t1\t50\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tAA\tGG' \
  $'n_prefix\t+/A\t1\t60\t+\tNA\tNA\tNA\tNA\tNA\tNA\t++\t++\t--' \
  $'known_indel\tAT/A\t1\t70\t+\tNA\tNA\tNA\tNA\tNA\tNA\tTT\tTT\t--' \
  >"$TESTDATA/tutorial/vcf_export_edges.hmp.txt"

# --- Helpers that run TASSEL 5 and record what they ran -----------------------

# sTASSEL.jar's manifest Class-Path lists its lib/ jars in a fixed order.
CLASSPATH_PLAIN="$T5/sTASSEL.jar"
CLASSPATH_FULL="$FULL_CLASSES:$T5/sTASSEL.jar"

manifest() { printf '%s\n' "$@" >>"$MANIFEST"; }

run_tassel5() {
  local classpath="$1"
  shift
  (cd "$TESTDATA" &&
    "$JAVA" -Xmx4g -cp "$classpath" net.maizegenetics.pipeline.TasselPipeline -maxThreads 4 "$@") \
    >>"$LOG" 2>&1 || die "TASSEL 5 failed; see $LOG"
}

# section NAME DESCRIPTION: outputs go in golden/NAME/, which TASSEL 5 will not create.
section() {
  mkdir -p "$GOLDEN/$1"
  manifest "" "[$1]" "description = $2"
}

tassel5() {
  manifest "command = tassel5 $*"
  run_tassel5 "$CLASSPATH_PLAIN" "$@"
}

tassel5_full() {
  manifest "command = tassel5-full $*"
  run_tassel5 "$CLASSPATH_FULL" "$@"
}

# output PATH TOLERANCE: PATH is relative to testdata/.
output() {
  [[ -s "$TESTDATA/$1" ]] || die "TASSEL 5 did not write $1; see $LOG"
  manifest "$1 = $2"
}

# Filter plugins also emit their filter definition as a second dataset; it is not a golden.
discard() {
  rm -f "$TESTDATA/$1"
}

# expect_vcf_site FILE ID REF ALT GTS
# GTS is the three genotype columns joined by tabs. Fails when a writeToVCF branch is missing.
expect_vcf_site() {
  local file="$1" id="$2" ref="$3" alt="$4" gts="$5"
  local line got_ref got_alt g1 g2 g3 got_gts
  line="$(awk -F '\t' -v id="$id" '$3 == id { print; exit }' "$file")"
  [[ -n "$line" ]] || die "VCF export has no site $id in $file"
  IFS=$'\t' read -r _ _ _ got_ref got_alt _ _ _ _ g1 g2 g3 <<<"$line"
  got_gts="${g1}"$'\t'"${g2}"$'\t'"${g3}"
  [[ "$got_ref" == "$ref" && "$got_alt" == "$alt" && "$got_gts" == "$gts" ]] ||
    die "site $id is REF=$got_ref ALT=$got_alt GT=$got_gts; expected REF=$ref ALT=$alt GT=$gts"
}

# The edges VCF is the oracle for the indel, N, and no-allele branches. Check the REF/ALT
# (and the GT columns that distinguish them) rather than trusting that TASSEL wrote a file.
check_vcf_export_edges() {
  local file="$1" n
  [[ -s "$file" ]] || die "VCF export edges file is empty"
  grep -q '^##fileformat=VCFv4.0' "$file" || die "edges VCF is missing ##fileformat=VCFv4.0"
  grep -q 'Reference allele is not known' "$file" ||
    die "edges VCF is missing the unknown-reference ##Tassel line"
  if grep -q '^##contig=' "$file"; then die "edges VCF has a ##contig line"; fi
  if grep -q '^##SAMPLE=' "$file"; then die "edges VCF has a ##SAMPLE line"; fi
  n="$(awk '$0 !~ /^#/ { c++ } END { print c + 0 }' "$file")"
  [[ "$n" == 7 ]] || die "edges VCF has $n data lines, expected 7"
  expect_vcf_site "$file" plus_minus '+' '-' $'0/0\t0/0\t1/1'
  expect_vcf_site "$file" n_in_alleles A C $'0/0\t0/0\t1/1'
  expect_vcf_site "$file" triallelic A 'G,C' $'0/0\t0/2\t1/1'
  expect_vcf_site "$file" no_calls '.' '.' $'./.\t./.\t./.'
  expect_vcf_site "$file" unlisted A 'G,C' $'0/0\t0/0\t1/1'
  expect_vcf_site "$file" n_prefix '+' 'N,NA' $'0/0\t0/0\t1/1'
  expect_vcf_site "$file" known_indel AT A $'0/0\t0/0\t1/1'
}

# --- Goldens -------------------------------------------------------------------

rm -rf "$GOLDEN"
mkdir -p "$GOLDEN"
: >"$LOG"
H="tutorial/mdp_genotype.hmp.txt"
V="tutorial/mdp_genotype_w_depth.vcf"

manifest \
  "# TASSEL 5 golden outputs for TASSEL 6 parity tests (ADR-0009)." \
  "# Generated by scripts/generate-goldens.sh. Do not edit by hand." \
  "#" \
  "# Paths are relative to testdata/, and every command runs from there." \
  "#   tassel5      = java -Xmx4g -cp sTASSEL.jar net.maizegenetics.pipeline.TasselPipeline -maxThreads 4" \
  "#   tassel5-full = tassel5 with scripts/golden/DoubleFormat.java ahead of sTASSEL.jar on the" \
  "#                  classpath, so table doubles are written with Double.toString instead of" \
  "#                  TASSEL's five-significant-digit format. Nothing else changes." \
  "# -maxThreads 4 fixes Centered IBS kinship's site blocks at 1,000 sites, so its float" \
  "# summation order does not depend on the host's core count." \
  "#" \
  "# Tolerances, per output file:" \
  "#   exact          identical text after normalizing line endings" \
  "#   abs X          every number within X absolutely; text cells identical" \
  "#   rel X          every number within X relative to the golden value; text cells identical" \
  "#   ... up-to-sign each numeric column may be negated as a whole before comparing" \
  "#   ... rows-unordered the data rows may be compared in any order: both sides are put in a" \
  "#                  canonical row order first, so nothing a row's position implies is compared" \
  "#                  (TASSEL 5's PCA eigenvalues come out of an SVD that does not sort them)" \
  "#   format         headers, row labels, and layout identical; numbers are checked in the" \
  "#                  .full.txt twin (PCA signs are arbitrary, so TASSEL-rounded text can differ)" \
  "#   phenotype      a <Phenotype> file: the tag and types lines identical, the names line" \
  "#                  identical except the taxa column's name (TASSEL 5 always writes Taxa), and" \
  "#                  every data cell identical as text (ADR-0023); with 'rel X' and/or 'abs Y'" \
  "#                  appended, a data or covariate cell may instead be a number within X relative" \
  "#                  to the golden or within Y of it" \
  "# A GLM report is a pair: the rounded file (format) and its .full.txt twin, whose tolerance is" \
  "# the whole-file rel value. Per ADR-0029 the names, Chr, Pos, Obs, minorObs, and *_df columns are" \
  "# exact (a count that differs is far outside any rel 1e-8), and a p below 1e-300 or printed 0.0 is" \
  "# compared as log10 p; the loader applies both, because this file holds one tolerance per file." \
  "" \
  "[tassel5]" \
  "version = $TASSEL5_VERSION" \
  "tag = $TASSEL5_TAG" \
  "artifact = $TASSEL5_URL" \
  "sha256 = $TASSEL5_SHA256" \
  "java = 21+" \
  "" \
  "[inputs]" \
  "artifact = $TUTORIAL_URL" \
  "sha256 = $TUTORIAL_SHA256"
for f in "${INPUTS[@]}"; do
  manifest "tutorial/$f = dataFiles/Tutorial/$f"
done
for f in "${PHENOTYPE_INPUTS[@]}"; do
  manifest "phenotype/$(basename "$f") = dataFiles/$f"
done
manifest "tutorial/vcf_export_edges.hmp.txt = generated by scripts/generate-goldens.sh"
for f in glm_synthetic.hmp.txt glm_synthetic_traits.txt glm_synthetic_factor.txt glm_synthetic_collinear.txt; do
  manifest "association/$f = hand-built, see golden/glm-synthetic/README.md"
done
manifest "association/glm_phenotype_only.txt = hand-built, see golden/glm-phenotype-only/README.md"

echo "Generating goldens with TASSEL $TASSEL5_VERSION (Java $JAVA_MAJOR)"

section hapmap "HapMap export of the tutorial HapMap (reader and writer round trip)"
tassel5 -h "$H" -export golden/hapmap/mdp_genotype -exportType Hapmap
output golden/hapmap/mdp_genotype.hmp.txt exact

section vcf "HapMap export of the tutorial VCF (VCF reader)"
tassel5 -vcf "$V" -export golden/vcf/mdp_genotype_w_depth -exportType Hapmap
output golden/vcf/mdp_genotype_w_depth.hmp.txt exact

section vcf-export "VCF export of the tutorial HapMap, an include-list, and indel, N, and no-allele edges"
tassel5 -h "$H" -export golden/vcf-export/mdp_genotype -exportType VCF
output golden/vcf-export/mdp_genotype.vcf exact
tassel5 -h "$H" -FilterTaxaBuilderPlugin -taxaList W22,B73,CML247,33-16 -includeTaxa true -endPlugin \
  -export "golden/vcf-export/include-list,golden/vcf-export/include-list.filter" -exportType VCF
discard golden/vcf-export/include-list.filter.json.gz
output golden/vcf-export/include-list.vcf exact
tassel5 -h tutorial/vcf_export_edges.hmp.txt -export golden/vcf-export/vcf_export_edges -exportType VCF
output golden/vcf-export/vcf_export_edges.vcf exact
check_vcf_export_edges "$TESTDATA/golden/vcf-export/vcf_export_edges.vcf"

S=golden/summary
section summary "GenotypeSummaryPlugin defaults: overview, allele, site, and taxa tables"
tassel5 -h "$H" -GenotypeSummaryPlugin -endPlugin \
  -export "$S/overview,$S/alleles,$S/sites,$S/taxa" -exportType Table
tassel5_full -h "$H" -GenotypeSummaryPlugin -endPlugin \
  -export "$S/overview.full,$S/alleles.full,$S/sites.full,$S/taxa.full" -exportType Table
for t in overview alleles sites taxa; do output "$S/$t.txt" exact; done
for t in overview alleles sites taxa; do output "$S/$t.full.txt" "abs 1e-12"; done

# filter_sites NAME PLUGIN-ARGS...
filter_sites() {
  local name="$1"
  shift
  tassel5 -h "$H" -FilterSiteBuilderPlugin "$@" -endPlugin \
    -export "golden/filter-sites/$name,golden/filter-sites/$name.filter" -exportType Hapmap
  discard "golden/filter-sites/$name.filter.json.gz"
  output "golden/filter-sites/$name.hmp.txt" exact
}

section filter-sites "FilterSiteBuilderPlugin parameter sets; outputs are HapMap exports"
filter_sites min-count-270 -siteMinCount 270
filter_sites maf-0.05-0.4 -siteMinAlleleFreq 0.05 -siteMaxAlleleFreq 0.4
filter_sites site-range-100-599 -siteRangeFilterType SITES -startSite 100 -endSite 599
filter_sites position-range-5-6 -siteRangeFilterType POSITIONS \
  -startChr 5 -startPos 1000000 -endChr 6 -endPos 50000000
filter_sites remove-minor-snp-states -removeMinorSNPStates true
filter_sites combined -siteMinCount 250 -siteMinAlleleFreq 0.1 -maxHeterozygous 0.02

# filter_taxa NAME PLUGIN-ARGS...
filter_taxa() {
  local name="$1"
  shift
  tassel5 -h "$H" -FilterTaxaBuilderPlugin "$@" -endPlugin \
    -export "golden/filter-taxa/$name,golden/filter-taxa/$name.filter" -exportType Hapmap
  discard "golden/filter-taxa/$name.filter.json.gz"
  output "golden/filter-taxa/$name.hmp.txt" exact
}

section filter-taxa "FilterTaxaBuilderPlugin parameter sets; outputs are HapMap exports"
filter_taxa min-not-missing-0.95 -minNotMissing 0.95
filter_taxa het-0.005-0.02 -minHeterozygous 0.005 -maxHeterozygous 0.02
filter_taxa include-list -taxaList W22,B73,CML247,33-16 -includeTaxa true
filter_taxa exclude-list -taxaList W22,B73,CML247,33-16 -includeTaxa false

section kinship "KinshipPlugin Centered_IBS with defaults; TASSEL square-matrix text (float values)"
tassel5 -h "$H" -KinshipPlugin -method Centered_IBS -endPlugin \
  -export golden/kinship/centered_ibs -exportType SqrMatrix
output golden/kinship/centered_ibs.txt "abs 8e-6"

P=golden/pca
section pca "PrincipalComponentsPlugin defaults (covariance, 5 components); PCs are stored as floats"
tassel5 -h "$H" -PrincipalComponentsPlugin -endPlugin \
  -export "$P/pcs,$P/eigenvalues,$P/eigenvectors"
tassel5_full -h "$H" -PrincipalComponentsPlugin -endPlugin \
  -export "$P/pcs.full,$P/eigenvalues.full,$P/eigenvectors.full"
output "$P/pcs.txt" format
output "$P/eigenvalues.txt" "exact rows-unordered"
output "$P/eigenvectors.txt" format
output "$P/pcs.full.txt" "abs 1e-5 up-to-sign"
output "$P/eigenvalues.full.txt" "rel 1e-9 rows-unordered"
output "$P/eigenvectors.full.txt" "abs 1e-9 up-to-sign"

F=golden/phenotype-format
section phenotype-format "<Phenotype> files loaded (-r) and exported as Phenotype: the format reader and writer"
for f in basic_phenotype_test earht_with_rep; do
  tassel5 -r "phenotype/$f.txt" -export "$F/$f" -exportType Phenotype
  output "$F/$f.txt" phenotype
done

# Legacy <Trait>, <Covariate>, <Factor>, and <Header> files, loaded with -t. TASSEL 5 sorts a
# <Header> file's rows and columns; TASSEL 6 keeps first-appearance order, so rows are unordered.
L=golden/phenotype-legacy
section phenotype-legacy "Legacy <Trait>, <Covariate>, and <Header> files loaded (-t) and exported as Phenotype"
for f in tutorial/mdp_traits tutorial/mdp_population_structure phenotype/version4format_test \
  phenotype/version4format_covariate phenotype/mdp_traits_nomissing; do
  tassel5 -t "$f.txt" -export "$L/$(basename "$f")" -exportType Phenotype
  output "$L/$(basename "$f").txt" "phenotype rows-unordered"
done

# TransformDataPlugin standardizes with a float mean, which moves a standardized value by up to
# ~8e-7. TASSEL 6 computes in double and narrows, so numbers are compared within 1e-6 relative or
# 1e-6 absolute rather than as text; the absolute bound covers values near zero.
X=golden/phenotype-transform
section phenotype-transform "TransformDataPlugin log, power, and standardize, alone and within factor levels"
transform() {
  local name="$1" input="$2"
  shift 2
  tassel5 -t "$input" -TransformDataPlugin "$@" -endPlugin -export "$X/$name" -exportType Phenotype
  output "$X/$name.txt" "phenotype rel 1e-6 abs 1e-6"
}
transform mdp_traits.log-natural tutorial/mdp_traits.txt -log natural
transform mdp_traits.log-base_10 tutorial/mdp_traits.txt -log base_10
transform mdp_traits.power-0.5 tutorial/mdp_traits.txt -power 0.5
transform mdp_traits.standardize tutorial/mdp_traits.txt -standardize true
tassel5 -r phenotype/earht_with_rep.txt -TransformDataPlugin -standardize true -factors rep -endPlugin \
  -export "$X/earht_with_rep.standardize-rep" -exportType Phenotype
output "$X/earht_with_rep.standardize-rep.txt" "phenotype rel 1e-6 abs 1e-6"

# -intersect and -union join the traits and the population structure. They need both inputs as one
# dataset, so each is loaded in its own fork and combined; with plain -t and -q the plugin sees one
# item at a time and waits on a dialog that headless TASSEL 5 never shows. The two share no factors and
# no names, so TASSEL 6's defined order (left, then right) and TASSEL 5's sorted order differ only in
# row order, which `rows-unordered` ignores (ADR-0024).
J=golden/phenotype-join
section phenotype-join "TASSEL 5 -intersect and -union of the tutorial traits and population structure"
for mode in intersect union; do
  tassel5 -fork1 -t tutorial/mdp_traits.txt -fork2 -q tutorial/mdp_population_structure.txt \
    -combine3 -input1 -input2 "-$mode" -export "$J/mdp_traits.q.$mode" -exportType Phenotype \
    -runfork1 -runfork2
  output "$J/mdp_traits.q.$mode.txt" "phenotype rows-unordered"
done

# Colt 1.2.0 cern.jet.stat.Gamma, which LinearModelUtils.Ftest and every TASSEL 5 GLM and MLM p-value call.
# The grid is in scripts/golden/Distributions.java: incompleteBeta over shapes 1e-8..1e7 and x from 1e-300 to 1
# (power series, the a/b swap and its 1 - MACHEP cap, log-space scaling, underflow to 0), the F-test arguments for
# integer and non-integer degrees of freedom, incompleteGamma and its complement over shapes 1e-8..1e5 and x up to
# 1e8, gamma and logGamma from -1e5 to 1e300, and the domain errors, written as `error`. Rows whose value moves by
# more than 2e-14 relative under a 2-ulp nudge of any argument are left out: they are ill-conditioned, and platform
# exp, log, and pow differ enough in the last bit that no port can hold rel 1e-13 on them.
S=golden/stats
section stats "Colt 1.2.0 Gamma (incompleteBeta, incompleteGamma, gamma, logGamma) and Probability (F, t, chi-square tails); commons-math3 3.4.1 quantiles"
# The same program writes the F, Student's t, and chi-square files. distributions.txt is Colt's cdf and survival
# values (F through the incomplete beta calls TASSEL 5 makes) over a grid of degrees of freedom and x with tails
# down to 1e-300, rows outside [1e-300, 1] left out. distribution-quantiles.txt is commons-math3 3.4.1's
# inverseCumulativeProbability, the Brent solve TASSEL 5 uses, for quantiles up to 1e6.
manifest "command = java -cp colt-1.2.0.jar:commons-math3-3.4.1.jar:Distributions.class Distributions $S/special-functions.txt $S/distributions.txt $S/distribution-quantiles.txt"
(cd "$TESTDATA" && "$JAVA" -cp "$COLT:$COMMONS_MATH:$DIST_CLASSES" Distributions \
  "$S/special-functions.txt" "$S/distributions.txt" "$S/distribution-quantiles.txt") ||
  die "Distributions failed"
output "$S/special-functions.txt" "rel 1e-13"
output "$S/distributions.txt" "rel 1e-12"
output "$S/distribution-quantiles.txt" "abs 1e-9"

# --- GLM (FixedEffectLMPlugin), issue #188 ---------------------------------------------------------------
# Every GLM run writes two reports, the site statistics and the effect estimates, each rounded as -export
# writes them (compared with `format`) and at full Double.toString precision through the tassel5_full
# DoubleFormat shadow (compared with rel 1e-8, the ADR-0029 starting value for GLM statistics). The
# genotype, the traits, and the optional covariates are joined through forks, as for the phenotype joins:
# with plain -h, -t, and -q the plugin would see one item at a time.
Q="tutorial/mdp_population_structure.txt"
TR="tutorial/mdp_traits.txt"
H910="tutorial/mdp_genotype_chr9_10.hmp.txt"
GLM_REL="rel 1e-8"

# glm_run NAME [PLUGIN-ARGS...]: uses GLM_PRE (the forks and joins) and GLM_POST (the -runfork flags),
# writes golden/<section>/NAME.{stats,effects}[.full].txt.
glm_run() {
  local name="$1" dir="$GLM_DIR"
  shift
  tassel5 "${GLM_PRE[@]}" -FixedEffectLMPlugin "$@" -endPlugin \
    -export "$dir/$name.stats,$dir/$name.effects" "${GLM_POST[@]}"
  tassel5_full "${GLM_PRE[@]}" -FixedEffectLMPlugin "$@" -endPlugin \
    -export "$dir/$name.stats.full,$dir/$name.effects.full" "${GLM_POST[@]}"
  for r in stats effects; do
    output "$dir/$name.$r.txt" format
    output "$dir/$name.$r.full.txt" "$GLM_REL"
  done
}

GLM_TRAITS_Q=(-fork2 -t "$TR" -fork3 -q "$Q" -excludeLastTrait -combine4 -input1 -input2 -input3 -intersect)
GLM_TRAITS_Q_POST=(-runfork1 -runfork2 -runfork3)

GLM_DIR=golden/glm-tutorial
section glm-tutorial "FixedEffectLMPlugin defaults, marker + Q1 + Q2 (-excludeLastTrait drops Q3), three traits, tutorial HapMap"
GLM_PRE=(-fork1 -h "$H" "${GLM_TRAITS_Q[@]}")
GLM_POST=("${GLM_TRAITS_Q_POST[@]}")
glm_run mdp_traits

GLM_DIR=golden/glm-tutorial-maf
section glm-tutorial-maf "The user manual's GLM recipe: -siteMinAlleleFreq 0.05 on the genotype, then marker + Q"
GLM_PRE=(-fork1 -h "$H" -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin "${GLM_TRAITS_Q[@]}")
glm_run mdp_traits

# The smaller chr 9 and 10 table (414 sites) that TASSEL 5's own DiscreteSitesTest uses keeps the option
# and dosage goldens small; the full tutorial table is used by the two sets above.
GLM_DIR=golden/glm-no-covariates
section glm-no-covariates "Marker-only model (no -q fork), three traits, chr 9 and 10 HapMap"
GLM_PRE=(-fork1 -h "$H910" -fork2 -t "$TR" -combine3 -input1 -input2 -intersect)
GLM_POST=(-runfork1 -runfork2)
glm_run mdp_traits

GLM_DIR=golden/glm-options
section glm-options "One run per option, each on marker + Q (chr 9 and 10 HapMap, three traits)"
GLM_PRE=(-fork1 -h "$H910" "${GLM_TRAITS_Q[@]}")
GLM_POST=("${GLM_TRAITS_Q_POST[@]}")
glm_run min-class-size-5 -minClassSize 5
glm_run biallelic-only -biallelicOnly true
glm_run max-p-0.01 -maxP 0.01
glm_run append-add-dom -appendAddDom true

# NumericalGenotypePlugin writes the reference-probability component, which -genotypeComponent
# ReferenceProbability then reads; the chain works headless, so this set is kept.
GLM_DIR=golden/glm-refprob
section glm-refprob "-genotypeComponent ReferenceProbability on NumericalGenotypePlugin's output (shorter report layouts), chr 9 and 10 HapMap"
GLM_PRE=(-fork1 -h "$H910" -NumericalGenotypePlugin -endPlugin -fork2 -t "$TR" -combine3 -input1 -input2 -intersect)
GLM_POST=(-runfork1 -runfork2)
glm_run mdp_traits -genotypeComponent ReferenceProbability

# Synthetic inputs (testdata/association/, see golden/glm-synthetic/README.md for what each one is for).
A=association
SYN="$A/glm_synthetic.hmp.txt"
GLM_DIR=golden/glm-synthetic
section glm-synthetic "Synthetic edge cases: multiallelic, indel, monomorphic, rare classes, missing calls, collinear covariates, one factor"
GLM_POST=(-runfork1 -runfork2)
GLM_PRE=(-fork1 -h "$SYN" -fork2 -r "$A/glm_synthetic_traits.txt" -combine3 -input1 -input2 -intersect)
glm_run model
glm_run min-class-size-3 -minClassSize 3
glm_run biallelic-only -biallelicOnly true
glm_run append-add-dom -appendAddDom true
GLM_PRE=(-fork1 -h "$SYN" -fork2 -r "$A/glm_synthetic_factor.txt" -combine3 -input1 -input2 -intersect)
glm_run factor-only
GLM_PRE=(-fork1 -h "$SYN" -fork2 -r "$A/glm_synthetic_collinear.txt" -combine3 -input1 -input2 -intersect)
glm_run collinear
cp "$ROOT/scripts/golden/glm-synthetic-README.md" "$GOLDEN/glm-synthetic/README.md"
output golden/glm-synthetic/README.md exact

# Phenotype-only GLM (BLUEs and the taxon ANOVA), issue #193. With no genotype in the input the plugin fits mean +
# factors + covariates + taxa per trait and writes two datasets, the BLUEs (a phenotype, floats) and the ANOVA report.
# The ANOVA is kept rounded (compared with `format`) and at full Double.toString precision (rel 1e-8, ADR-0029), the
# BLUEs once (floats), compared as a phenotype with rows in any order because TASSEL 5 sorts the taxa.
PO_DIR=golden/glm-phenotype-only
GLM_REL="rel 1e-8"

# pheno_only_run NAME LOADER FILE: loads FILE with LOADER (-r or -t) and writes NAME.blues.txt, NAME.anova.txt, and
# NAME.anova.full.txt.
pheno_only_run() {
  local name="$1" loader="$2" file="$3"
  tassel5 "$loader" "$file" -FixedEffectLMPlugin -phenoOnly true -endPlugin \
    -export "$PO_DIR/$name.blues,$PO_DIR/$name.anova"
  tassel5_full "$loader" "$file" -FixedEffectLMPlugin -phenoOnly true -endPlugin \
    -export "$PO_DIR/$name.blues.full,$PO_DIR/$name.anova.full"
  cmp -s "$TESTDATA/$PO_DIR/$name.blues.txt" "$TESTDATA/$PO_DIR/$name.blues.full.txt" ||
    die "$name BLUEs differ between the plain and full-precision runs"
  rm -f "$TESTDATA/$PO_DIR/$name.blues.full.txt"
  output "$PO_DIR/$name.blues.txt" "phenotype rel 1e-6 rows-unordered"
  output "$PO_DIR/$name.anova.txt" format
  output "$PO_DIR/$name.anova.full.txt" "$GLM_REL"
}

section glm-phenotype-only "FixedEffectLMPlugin -phenoOnly true: taxon BLUEs and the ANOVA, replicated, unreplicated, and with a covariate"
pheno_only_run earht_with_rep -r phenotype/earht_with_rep.txt
pheno_only_run mdp_traits -t tutorial/mdp_traits.txt
pheno_only_run synthetic -r association/glm_phenotype_only.txt
cp "$ROOT/scripts/golden/glm-phenotype-only-README.md" "$GOLDEN/glm-phenotype-only/README.md"
output golden/glm-phenotype-only/README.md exact

# TASSEL 5's own expected file for DiscreteSitesTest.testUnreplicatedTaxa (EarHT on the chr 9 and 10 table, no
# covariates): copied as is, and the script checks that this pipeline reproduces it byte for byte.
section glm-t5-expected "TASSEL 5's ExpectedResults/discreteSitesFELM.txt (EarHT, no covariates, chr 9 and 10 HapMap)"
cp "$CACHE/dataFiles/ExpectedResults/discreteSitesFELM.txt" "$GOLDEN/glm-t5-expected/discreteSitesFELM.txt"
manifest "source = dataFiles/ExpectedResults/discreteSitesFELM.txt in the inputs artifact"
output golden/glm-t5-expected/discreteSitesFELM.txt exact
awk -F'\t' 'NR == 1 || $1 == "EarHT"' "$GOLDEN/glm-no-covariates/mdp_traits.stats.txt" |
  cmp -s - "$GOLDEN/glm-t5-expected/discreteSitesFELM.txt" ||
  die "glm-no-covariates EarHT rows differ from TASSEL 5's discreteSitesFELM.txt"

# No replicate golden: GLM refuses replicated taxa. The command and the error are the documented baseline
# for the replicate handling TASSEL 6 adds (#191); nothing is written, so the run is allowed to fail.
section glm-replicates "FixedEffectLMPlugin on replicated taxa fails in TASSEL 5; no output files"
REP_ARGS=(-fork1 -h "$H" -fork2 -r phenotype/earht_with_rep.txt -combine3 -input1 -input2 -intersect
  -FixedEffectLMPlugin -endPlugin -export golden/glm-replicates/never_written -runfork1 -runfork2)
manifest "command = tassel5 ${REP_ARGS[*]}"
(cd "$TESTDATA" && "$JAVA" -Xmx4g -cp "$CLASSPATH_PLAIN" net.maizegenetics.pipeline.TasselPipeline -maxThreads 4 \
  "${REP_ARGS[@]}") >"$CACHE/replicates.log" 2>&1 || true
cat "$CACHE/replicates.log" >>"$LOG"
REP_ERROR="$(grep -m1 'There are more phenotype observations than taxa with genotypes' "$CACHE/replicates.log" |
  sed 's/^.*AbstractPlugin - //')"
[[ -n "$REP_ERROR" ]] || die "GLM on replicated taxa did not fail with the expected error; see $CACHE/replicates.log"
[[ ! -e "$GOLDEN/glm-replicates/never_written1.txt" ]] || die "replicated GLM wrote output"
manifest "error = $REP_ERROR"
rmdir "$GOLDEN/glm-replicates" 2>/dev/null || true

# The cells of the GLM goldens where TASSEL 6 departs from TASSEL 5 on purpose (ADR-0029, #191), as a reviewed list of
# rules; GlmGoldenTest sets exactly these cells aside. Appended as it is so regeneration keeps it.
{ echo; cat "$ROOT/scripts/golden/glm-departures.manifest"; } >>"$MANIFEST"

leftover="$(cd "$TESTDATA" && find golden -type f ! -name MANIFEST | sort |
  while read -r f; do grep -qF "$f = " golden/MANIFEST || echo "$f"; done)"
[[ -z "$leftover" ]] || die "unlisted files in testdata/golden: $leftover"

echo "Wrote $(cd "$TESTDATA" && find golden -type f | wc -l | tr -d ' ') files to testdata/golden ($(du -sh "$GOLDEN" | cut -f1))"
