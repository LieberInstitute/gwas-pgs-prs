#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: build_prsice_base_from_gwas_vcf.sh --config CONFIG --manifest BASE.tsv [--dry-run]

Required manifest columns:
label gwas_bcf output_base target_bcf

Use target_bcf=. to skip exact target filtering.
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
TMPDIR_ROOT=${TMPDIR_ROOT:-tmp}
require_cmd "$BCFTOOLS"
require_fields "$MANIFEST" label gwas_bcf output_base target_bcf
mkdir -p "$TMPDIR_ROOT"

manifest_rows "$MANIFEST" | while IFS=$'\t' read -r label gwas_bcf output_base target_bcf rest; do
  require_file "$gwas_bcf"
  mkdir -p "$(dirname "$output_base")"
  tmp_tsv="$TMPDIR_ROOT/${label}.prsice_base.tmp.tsv"

  if [ "$target_bcf" != "." ]; then
    require_file "$target_bcf"
    target_ids="$TMPDIR_ROOT/${label}.target.variant_ids"
    msg "write exact target variant IDs $label"
    run_bash "'$BCFTOOLS' query -f '%CHROM:%POS:%REF:%ALT\\n' '$target_bcf' | sort -u > '$target_ids'"
    filter_cmd="awk -F '\\t' 'BEGIN{OFS=\"\\t\"} NR==FNR{keep[\$1]=1; next} FNR==1{print; next} keep[\$8]' '$target_ids' '$tmp_tsv'"
  else
    filter_cmd="cat '$tmp_tsv'"
  fi

  msg "extract PRSice base $label"
  run_bash "printf 'CHR\\tBP\\tSNP\\tA1\\tA2\\tBETA\\tP\\tVARIANT_ID\\n' > '$tmp_tsv' && '$BCFTOOLS' query -f '%CHROM\\t%POS\\t%ID\\t%ALT\\t%REF\\t[%ES]\\t[%LP]\\t%CHROM:%POS:%REF:%ALT\\n' '$gwas_bcf' | awk -F '\\t' 'BEGIN{OFS=\"\\t\"} {chr=\$1; sub(/^chr/,\"\",chr); p=(\$7==\"\" || \$7==\".\" ? \"NA\" : 10^(-\$7)); snp=(\$3==\".\" || \$3==\"\" ? \$8 : \$3); if(\$6!=\"\" && \$6!=\".\" && p!=\"NA\") print chr,\$2,snp,\$4,\$5,\$6,p,\$8}' >> '$tmp_tsv' && $filter_cmd | cut -f1-7 | gzip -c > '$output_base'"
done
