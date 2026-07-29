#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_mbv_prs_23andme_impact.sh [stage]

Stages:
  prepare   build canonical bases, harmonized target, manifests, and provenance
  smoke     run chromosome 22 C+T smoke tests for both target setups
  ct        run all historical and harmonized Shizhong-style C+T scores
  prsice    run fixed-threshold PRSice scores
  graphpred apply existing LDGM/GraphPred PGS BCFs with bcftools +score
  collect   build the long score registry
  validate  compare no23 reruns with archived scores
  compare   calculate 23andMe-impact metrics and generated report
  all       run all stages in order

Environment overrides: OUT_ROOT, LD_ROOT, THREADS, PLINK, PLINK2, PRSICE,
BCFTOOLS, GWAS_ROOT, ARCHIVE_ROOT.
USAGE
}

STAGE=${1:-all}
case "$STAGE" in
  prepare|smoke|ct|prsice|graphpred|collect|validate|compare|all) ;;
  -h|--help) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCRIPT_DIR="$ROOT/scripts"
OUT_ROOT=${OUT_ROOT:-$ROOT/generated/prs_23andme_impact}
LD_ROOT=${LD_ROOT:-/dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong}
GWAS_ROOT=${GWAS_ROOT:-/home/gpertea/work/ref/GWAS}
ARCHIVE_ROOT=${ARCHIVE_ROOT:-/home/gpertea/work/ref/GWAS/mbv-prs-generated-outputs/mbv-prs}
THREADS=${THREADS:-8}
PLINK=${PLINK:-/opt/sw/bin/plink}
PLINK2=${PLINK2:-/opt/sw/bin/plink2}
PRSICE=${PRSICE:-/opt/sw/bin/PRSice_linux}
BCFTOOLS=${BCFTOOLS:-/opt/sw/bin/bcftools}

LD_PATTERN="$LD_ROOT/EUR_noIndels_chr{CHR}_maf0.01_hg38"
LD_SOURCE_COMBINED="$LD_ROOT/EUR_noIndels_all_maf0.01_hg38"
LD_COMBINED="$LD_ROOT/EUR_noIndels_all_maf0.01_hg38.mbv_qc_canonical"
TARGET_BCF="$ROOT/qc/mbv.qc.chr.bcf"
HIST_TARGET="$ROOT/shizhong-mbv-prs/shan-prs-keri/data/merged_maf01_snps"
TARGET_PREFIX="$OUT_ROOT/target/mbv_qc"

BD_NO23_BCF="$GWAS_ROOT/BD/bip2024_eur_no23andMe.hg38.bcf"
BD_FULL_BASE_SOURCE="$GWAS_ROOT/BD/bip2024_eur.hg38.prsice.tsv.gz"
MDD_NO23_BCF="$GWAS_ROOT/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf"
MDD_FULL_BASE_SOURCE="$GWAS_ROOT/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.prsice.tsv.gz"

BD_NO23_PGS="$GWAS_ROOT/BD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf"
BD_FULL_PGS="$GWAS_ROOT/BD/bip2024_eur.hg38.pgs.b5e-8.bcf"
MDD_NO23_PGS="$GWAS_ROOT/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf"
MDD_FULL_PGS="$GWAS_ROOT/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf"

BD_NO23_BASE="$OUT_ROOT/base/BD.no23andMe.tsv.gz"
BD_FULL_BASE="$OUT_ROOT/base/BD.full.pre-DENTIST.tsv.gz"
MDD_NO23_BASE="$OUT_ROOT/base/MDD.no23andMe.tsv.gz"
MDD_FULL_BASE="$OUT_ROOT/base/MDD.full.tsv.gz"

CT_MANIFEST="$OUT_ROOT/manifests/ct.tsv"
PRSICE_MANIFEST="$OUT_ROOT/manifests/prsice.tsv"
GRAPH_MANIFEST="$OUT_ROOT/manifests/graphpred.tsv"
EXPECTED_SAMPLES="$OUT_ROOT/provenance/expected_samples.txt"
REGISTRY="$OUT_ROOT/results/prs_score_registry.tsv"
HIST_THRESHOLDS='1e-08,1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1'
PRSICE_THRESHOLDS='5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1'

mkdir -p "$OUT_ROOT/logs"
SUCCESS_LOG="$OUT_ROOT/logs/commands_successful.tsv"
[[ -e "$SUCCESS_LOG" ]] || printf 'timestamp\tstage\n' > "$SUCCESS_LOG"

record_success() {
  printf '%s\t%s\n' "$(date '+%F %T')" "$1" >> "$SUCCESS_LOG"
}

require_file() {
  [[ -s "$1" ]] || { printf 'Missing required file: %s\n' "$1" >&2; exit 1; }
}

base_ready() {
  local base=$1 stats="$1.validation.tsv"
  [[ -s "$base" && -s "$stats" ]] || return 1
  awk -F '\t' '$1=="accepted_rows" && $2>=1000000{ok=1} END{exit !ok}' "$stats"
}

prepare() {
  mkdir -p "$OUT_ROOT"/{base,target,manifests,provenance,ct,prsice,graphpred,results,smoke}
  for tool in "$PLINK" "$PLINK2" "$PRSICE" "$BCFTOOLS"; do require_file "$tool"; done
  for file in "$TARGET_BCF" "$HIST_TARGET.pgen" "$HIST_TARGET.pvar" "$HIST_TARGET.psam" \
              "$BD_NO23_BCF" "$BD_FULL_BASE_SOURCE" "$MDD_NO23_BCF" "$MDD_FULL_BASE_SOURCE" \
              "$BD_NO23_PGS" "$BD_FULL_PGS" "$MDD_NO23_PGS" "$MDD_FULL_PGS"; do
    require_file "$file"
  done
  for ext in bed bim fam; do require_file "$LD_SOURCE_COMBINED.$ext"; done
  for chr in {1..22}; do
    prefix=${LD_PATTERN//\{CHR\}/$chr}
    for ext in bed bim fam; do require_file "$prefix.$ext"; done
  done

  ## build one harmonized target from the final-QC dosage BCF
  if [[ ! -s "$TARGET_PREFIX.pgen" || ! -s "$TARGET_PREFIX.pvar" || ! -s "$TARGET_PREFIX.psam" ]]; then
    "$PLINK2" --bcf "$TARGET_BCF" dosage=DS --double-id --make-pgen \
      --out "$TARGET_PREFIX" > "$OUT_ROOT/logs/target-pgen.log" 2>&1
  fi
  if [[ ! -s "$TARGET_PREFIX.bgen" || ! -s "$TARGET_PREFIX.sample" ]]; then
    "$PLINK2" --pfile "$TARGET_PREFIX" --export bgen-1.2 bits=8 ref-first \
      --out "$TARGET_PREFIX" > "$OUT_ROOT/logs/target-bgen.log" 2>&1
  fi

  ## PRSice requires LD-reference IDs to match canonical base and target IDs
  if [[ ! -s "$LD_COMBINED.bed" || ! -s "$LD_COMBINED.bim" || ! -s "$LD_COMBINED.fam" ]]; then
    "$SCRIPT_DIR/build_prsice_ld_panel.sh" \
      --source-prefix "$LD_SOURCE_COMBINED" \
      --ld-pattern "$LD_PATTERN" \
      --target-pvar "$TARGET_PREFIX.pvar" \
      --out-prefix "$LD_COMBINED"
  fi
  for ext in bed bim fam; do require_file "$LD_COMBINED.$ext"; done

  "$BCFTOOLS" query -l "$TARGET_BCF" > "$EXPECTED_SAMPLES"
  awk 'NR>1{print $1}' "$HIST_TARGET.psam" | sort > "$OUT_ROOT/provenance/historical.samples.sorted"
  sort "$EXPECTED_SAMPLES" > "$OUT_ROOT/provenance/harmonized.samples.sorted"
  cmp -s "$OUT_ROOT/provenance/historical.samples.sorted" \
    "$OUT_ROOT/provenance/harmonized.samples.sorted" || {
      printf 'Historical and harmonized targets have different donor IDs\n' >&2
      exit 1
    }
  [[ $(wc -l < "$EXPECTED_SAMPLES") -eq 119 ]] || {
    printf 'Harmonized target does not contain 119 donors\n' >&2
    exit 1
  }

  ## PRSice identifies BGEN samples as FID_IID; use the matching dummy IDs
  {
    printf 'IID DUMMY\n'
    awk 'NR>2{print $1"_"$2,1}' "$TARGET_PREFIX.sample"
  } > "$OUT_ROOT/target/dummy.pheno"

  if ! base_ready "$BD_NO23_BASE"; then
    BCFTOOLS="$BCFTOOLS" "$SCRIPT_DIR/prs_base_from_gwas_bcf.sh" \
      --bcf "$BD_NO23_BCF" --out "$BD_NO23_BASE"
  fi
  if ! base_ready "$MDD_NO23_BASE"; then
    BCFTOOLS="$BCFTOOLS" "$SCRIPT_DIR/prs_base_from_gwas_bcf.sh" \
      --bcf "$MDD_NO23_BCF" --out "$MDD_NO23_BASE"
  fi
  if ! base_ready "$BD_FULL_BASE"; then
    perl "$SCRIPT_DIR/validate_prs_base.pl" --input "$BD_FULL_BASE_SOURCE" \
      --output "$BD_FULL_BASE" --stats "$BD_FULL_BASE.validation.tsv"
  fi
  if ! base_ready "$MDD_FULL_BASE"; then
    perl "$SCRIPT_DIR/validate_prs_base.pl" --input "$MDD_FULL_BASE_SOURCE" \
      --output "$MDD_FULL_BASE" --stats "$MDD_FULL_BASE.validation.tsv"
  fi

  write_manifests
  write_provenance
  record_success prepare
}

append_ct() {
  local run_id=$1 trait=$2 version=$3 setup=$4 base=$5 target=$6 out=$7
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t1000\t0.1\t%s\n' \
    "$run_id" "$trait" "$version" "$setup" "$base" "$target" "$LD_PATTERN" "$out" \
    "$HIST_THRESHOLDS" >> "$CT_MANIFEST"
}

append_prsice() {
  local run_id=$1 trait=$2 version=$3 base=$4 out=$5
  printf '%s\t%s\t%s\tharmonized\t%s\t%s\t%s\t%s\t%s\t250\t0.1\t%s\n' \
    "$run_id" "$trait" "$version" "$base" "$TARGET_PREFIX" "$LD_COMBINED" \
    "$OUT_ROOT/target/dummy.pheno" "$out" "$PRSICE_THRESHOLDS" >> "$PRSICE_MANIFEST"
}

append_graph() {
  local run_id=$1 trait=$2 version=$3 pgs=$4 out=$5
  printf '%s\t%s\t%s\tharmonized\t%s\t%s\t%s\n' \
    "$run_id" "$trait" "$version" "$TARGET_BCF" "$pgs" "$out" >> "$GRAPH_MANIFEST"
}

write_manifests() {
  printf 'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tBASE\tTARGET_PGEN\tLD_PATTERN\tOUT_DIR\tCLUMP_KB\tCLUMP_R2\tTHRESHOLDS\n' > "$CT_MANIFEST"
  for setup in historical harmonized; do
    target=$HIST_TARGET
    [[ "$setup" == harmonized ]] && target=$TARGET_PREFIX
    append_ct "BD_no23_${setup}" BD no23andMe "$setup" "$BD_NO23_BASE" "$target" "$OUT_ROOT/ct/BD/no23andMe/$setup"
    append_ct "BD_full_preDENTIST_${setup}" BD full "$setup" "$BD_FULL_BASE" "$target" "$OUT_ROOT/ct/BD/full/$setup"
    append_ct "MDD_no23_${setup}" MDD no23andMe "$setup" "$MDD_NO23_BASE" "$target" "$OUT_ROOT/ct/MDD/no23andMe/$setup"
    append_ct "MDD_full_${setup}" MDD full "$setup" "$MDD_FULL_BASE" "$target" "$OUT_ROOT/ct/MDD/full/$setup"
  done

  printf 'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tBASE\tTARGET_BGEN\tLD_PREFIX\tPHENO\tOUT_PREFIX\tCLUMP_KB\tCLUMP_R2\tTHRESHOLDS\n' > "$PRSICE_MANIFEST"
  append_prsice BD_no23_prsice BD no23andMe "$BD_NO23_BASE" "$OUT_ROOT/prsice/BD/no23andMe/BD"
  append_prsice BD_full_preDENTIST_prsice BD full "$BD_FULL_BASE" "$OUT_ROOT/prsice/BD/full/BD"
  append_prsice MDD_no23_prsice MDD no23andMe "$MDD_NO23_BASE" "$OUT_ROOT/prsice/MDD/no23andMe/MDD"
  append_prsice MDD_full_prsice MDD full "$MDD_FULL_BASE" "$OUT_ROOT/prsice/MDD/full/MDD"

  printf 'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tTARGET_BCF\tPGS_BCF\tOUT_FILE\n' > "$GRAPH_MANIFEST"
  append_graph BD_no23_graphpred BD no23andMe "$BD_NO23_PGS" "$OUT_ROOT/graphpred/BD.no23andMe.tsv"
  append_graph BD_full_preDENTIST_graphpred BD full "$BD_FULL_PGS" "$OUT_ROOT/graphpred/BD.full.tsv"
  append_graph MDD_no23_graphpred MDD no23andMe "$MDD_NO23_PGS" "$OUT_ROOT/graphpred/MDD.no23andMe.tsv"
  append_graph MDD_full_graphpred MDD full "$MDD_FULL_PGS" "$OUT_ROOT/graphpred/MDD.full.tsv"
}

write_provenance() {
  {
    printf 'bcftools\t'; "$BCFTOOLS" --version | head -n 1
    printf 'plink\t'; "$PLINK" --version 2>&1 | head -n 1
    printf 'plink2\t'; "$PLINK2" --version 2>&1 | head -n 1
    printf 'PRSice\t'; "$PRSICE" --version 2>&1 | head -n 1
    printf 'R\t'; Rscript --version 2>&1
    printf 'perl\t'; perl -e 'print "$^V\n"'
  } > "$OUT_ROOT/provenance/software_versions.tsv"

  sha256sum "$TARGET_BCF" "$HIST_TARGET.pgen" "$HIST_TARGET.pvar" "$HIST_TARGET.psam" \
    "$BD_NO23_BCF" "$BD_FULL_BASE_SOURCE" "$MDD_NO23_BCF" "$MDD_FULL_BASE_SOURCE" \
    "$BD_NO23_PGS" "$BD_FULL_PGS" "$MDD_NO23_PGS" "$MDD_FULL_PGS" \
    > "$OUT_ROOT/provenance/input_checksums.sha256"
}

smoke() {
  require_file "$CT_MANIFEST"
  SMOKE_MANIFEST="$OUT_ROOT/manifests/ct.smoke.tsv"
  head -n 1 "$CT_MANIFEST" > "$SMOKE_MANIFEST"
  awk -F '\t' -v OFS='\t' -v root="$OUT_ROOT/smoke" \
    '$1=="BD_no23_historical"{$8=root"/historical"; print}
     $1=="BD_no23_harmonized"{$8=root"/harmonized"; print}' \
    "$CT_MANIFEST" >> "$SMOKE_MANIFEST"
  PLINK="$PLINK" PLINK2="$PLINK2" "$SCRIPT_DIR/run_plink_ct_manifest.sh" \
    --manifest "$SMOKE_MANIFEST" --chromosomes 22 --threads "$THREADS" --resume
  record_success smoke
}

run_ct() {
  PLINK="$PLINK" PLINK2="$PLINK2" "$SCRIPT_DIR/run_plink_ct_manifest.sh" \
    --manifest "$CT_MANIFEST" --chromosomes 1-22 --threads "$THREADS" --resume
  record_success ct
}

run_prsice() {
  PRSICE="$PRSICE" "$SCRIPT_DIR/run_prsice_manifest.sh" \
    --manifest "$PRSICE_MANIFEST" --threads "$THREADS" --resume
  record_success prsice
}

run_graphpred() {
  BCFTOOLS="$BCFTOOLS" "$SCRIPT_DIR/run_graphpred_manifest.sh" \
    --manifest "$GRAPH_MANIFEST" --resume
  record_success graphpred
}

collect() {
  Rscript "$SCRIPT_DIR/compile_prs_registry.R" \
    --ct-manifest "$CT_MANIFEST" \
    --prsice-manifest "$PRSICE_MANIFEST" \
    --graph-manifest "$GRAPH_MANIFEST" \
    --expected-samples "$EXPECTED_SAMPLES" \
    --output "$REGISTRY"
  record_success collect
}

validate() {
  Rscript "$SCRIPT_DIR/validate_prs_reproduction.R" \
    --registry "$REGISTRY" \
    --archive-root "$ARCHIVE_ROOT" \
    --output "$OUT_ROOT/results/reproduction_validation.tsv"
  record_success validate
}

compare() {
  Rscript "$SCRIPT_DIR/compare_prs_23andme_impact.R" \
    --registry "$REGISTRY" \
    --out-dir "$OUT_ROOT/results" \
    --report "$OUT_ROOT/results/MBv_BD_MDD_23andMe_PRS_impact.md"
  record_success compare
}

if [[ "$STAGE" == all ]]; then
  prepare
  smoke
  run_ct
  run_prsice
  run_graphpred
  collect
  validate
  compare
else
  case "$STAGE" in
    prepare) prepare ;;
    smoke) smoke ;;
    ct) run_ct ;;
    prsice) run_prsice ;;
    graphpred) run_graphpred ;;
    collect) collect ;;
    validate) validate ;;
    compare) compare ;;
  esac
fi

printf '[%s] completed stage: %s\n' "$(date '+%F %T')" "$STAGE"
