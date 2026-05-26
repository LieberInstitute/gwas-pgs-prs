#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: compute_pgs_loadings.sh --config CONFIG --manifest PGS.tsv [--dry-run]

Required manifest columns:
label gwas_bcf ldgm_vcfs_file beta_cov max_alpha_hat2 exclude_expr output_bcf
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

BCFTOOLS=${BCFTOOLS:-bcftools}
THREADS=${THREADS:-4}
require_cmd "$BCFTOOLS"
require_fields "$MANIFEST" label gwas_bcf ldgm_vcfs_file beta_cov max_alpha_hat2 exclude_expr output_bcf

manifest_rows "$MANIFEST" | while IFS=$'\t' read -r label gwas_bcf ldgm_vcfs_file beta_cov max_alpha_hat2 exclude_expr output_bcf rest; do
  require_file "$gwas_bcf"
  require_file "$ldgm_vcfs_file"
  mkdir -p "$(dirname "$output_bcf")"

  args=(+pgs -Ob -o "$output_bcf" --ldgm-vcfs-file "$ldgm_vcfs_file")
  [ "$beta_cov" = "." ] || args+=(--beta-cov "$beta_cov")
  [ "$max_alpha_hat2" = "." ] || args+=(--max-alpha-hat2 "$max_alpha_hat2")
  [ "$exclude_expr" = "." ] || args+=(-e "$exclude_expr")

  msg "compute PGS loadings $label"
  run_cmd "$BCFTOOLS" "${args[@]}" "$gwas_bcf"
  run_cmd "$BCFTOOLS" index --threads "$THREADS" -f "$output_bcf"
done
