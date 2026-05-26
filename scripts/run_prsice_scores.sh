#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: run_prsice_scores.sh --config CONFIG --manifest PRSICE.tsv [--dry-run]

Required manifest columns:
label base_file target_prefix target_type ld_ref out_prefix
EOF
}

CONFIG=
MANIFEST=
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG=$2; shift 2 ;;
    --manifest) MANIFEST=$2; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[ -n "$CONFIG" ] || die "missing --config"
[ -n "$MANIFEST" ] || die "missing --manifest"
source_config "$CONFIG"

PRSICE=${PRSICE:-PRSice_linux}
THREADS=${THREADS:-4}
PRSICE_THRESHOLDS=${PRSICE_THRESHOLDS:-5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1}
PRSICE_CLUMP_KB=${PRSICE_CLUMP_KB:-250}
PRSICE_CLUMP_R2=${PRSICE_CLUMP_R2:-0.1}
PRSICE_CLUMP_P=${PRSICE_CLUMP_P:-1}
require_cmd "$PRSICE"
require_fields "$MANIFEST" label base_file target_prefix target_type ld_ref out_prefix

manifest_rows "$MANIFEST" | while IFS=$'\t' read -r label base_file target_prefix target_type ld_ref out_prefix rest; do
  require_file "$base_file"
  require_file "$target_prefix.sample"
  mkdir -p "$(dirname "$out_prefix")"
  pheno_file="$(dirname "$target_prefix")/dummy.pheno"
  require_file "$pheno_file"

  msg "PRSice $label"
  run_cmd "$PRSICE" \
    --base "$base_file" \
    --target "$target_prefix" \
    --type "$target_type" \
    --ld "$ld_ref" \
    --ld-type bed \
    --snp SNP \
    --chr CHR \
    --bp BP \
    --a1 A1 \
    --a2 A2 \
    --stat BETA \
    --pvalue P \
    --beta \
    --ignore-fid \
    --pheno "$pheno_file" \
    --pheno-col DUMMY \
    --no-regress \
    --all-score \
    --score sum \
    --fastscore \
    --print-snp \
    --bar-levels "$PRSICE_THRESHOLDS" \
    --clump-kb "$PRSICE_CLUMP_KB" \
    --clump-r2 "$PRSICE_CLUMP_R2" \
    --clump-p "$PRSICE_CLUMP_P" \
    --thread "$THREADS" \
    --out "$out_prefix"
done
