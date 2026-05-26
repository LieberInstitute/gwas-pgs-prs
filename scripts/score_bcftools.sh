#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: score_bcftools.sh --config CONFIG --manifest SCORE.tsv [--dry-run]

Required manifest columns:
label target_bcf score_bcf output_tsv
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
require_cmd "$BCFTOOLS"
require_fields "$MANIFEST" label target_bcf score_bcf output_tsv

manifest_rows "$MANIFEST" | while IFS=$'\t' read -r label target_bcf score_bcf output_tsv rest; do
  require_file "$target_bcf"
  require_file "$score_bcf"
  mkdir -p "$(dirname "$output_tsv")"
  msg "bcftools score $label"
  run_cmd "$BCFTOOLS" +score --use DS --sample-header --counts -o "$output_tsv" "$target_bcf" "$score_bcf"
done
