#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: run_prs_pipeline.sh --config CONFIG --stage STAGE [--dry-run]

Stages:
gwas target pgs bcftools-score prsice-base prsice prsice-models summarize compare all

The config should define manifest paths:
GWAS_MANIFEST PGS_MANIFEST SCORE_MANIFEST PRSICE_BASE_MANIFEST PRSICE_MANIFEST
PRSICE_MODEL_MANIFEST SCORE_LONG_MANIFEST COMPARE_MANIFEST
EOF
}

CONFIG=
STAGE=all
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG=$2; shift 2 ;;
    --stage) STAGE=$2; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[ -n "$CONFIG" ] || die "missing --config"
source_config "$CONFIG"

run_stage() {
  local stage=$1
  case "$stage" in
    gwas) "$SCRIPT_DIR/prepare_gwas_vcf.sh" --config "$CONFIG" --manifest "${GWAS_MANIFEST:?}" ${DRY_RUN:+--dry-run} ;;
    target) "$SCRIPT_DIR/prepare_target_qc.sh" --config "$CONFIG" ${DRY_RUN:+--dry-run} ;;
    pgs) "$SCRIPT_DIR/compute_pgs_loadings.sh" --config "$CONFIG" --manifest "${PGS_MANIFEST:?}" ${DRY_RUN:+--dry-run} ;;
    bcftools-score) "$SCRIPT_DIR/score_bcftools.sh" --config "$CONFIG" --manifest "${SCORE_MANIFEST:?}" ${DRY_RUN:+--dry-run} ;;
    prsice-base) "$SCRIPT_DIR/build_prsice_base_from_gwas_vcf.sh" --config "$CONFIG" --manifest "${PRSICE_BASE_MANIFEST:?}" ${DRY_RUN:+--dry-run} ;;
    prsice) "$SCRIPT_DIR/run_prsice_scores.sh" --config "$CONFIG" --manifest "${PRSICE_MANIFEST:?}" ${DRY_RUN:+--dry-run} ;;
    prsice-models) Rscript "$SCRIPT_DIR/export_prsice_model_tables.R" "${PRSICE_MODEL_MANIFEST:?}" ;;
    summarize) Rscript "$SCRIPT_DIR/summarize_prs_scores.R" "${SCORE_LONG_MANIFEST:?}" "${SUMMARY_OUT_DIR:?}" ;;
    compare) Rscript "$SCRIPT_DIR/compare_prs_methods.R" "${COMPARE_MANIFEST:?}" "${COMPARE_OUT_DIR:?}" ;;
    *) die "unknown stage: $stage" ;;
  esac
}

if [ "$STAGE" = "all" ]; then
  for s in gwas target pgs bcftools-score prsice-base prsice prsice-models summarize compare; do
    run_stage "$s"
  done
else
  run_stage "$STAGE"
fi
