#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: prepare_gwas_vcf.sh --config CONFIG --manifest GWAS.tsv [--dry-run]

Required manifest columns:
label source sample_name source_build source_fasta target_fasta chain columns_file columns_preset output_bcf

Use source_build=GRCh38 and chain=. to skip liftover.
Use columns_file=. when columns_preset is supplied.
Use columns_preset=. when columns_file is supplied.
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
TMPDIR_ROOT=${TMPDIR_ROOT:-tmp}
require_cmd "$BCFTOOLS"
require_fields "$MANIFEST" label source sample_name source_build source_fasta target_fasta chain columns_file columns_preset output_bcf

mkdir -p "$TMPDIR_ROOT"

manifest_rows "$MANIFEST" | while IFS=$'\t' read -r label source sample_name source_build source_fasta target_fasta chain columns_file columns_preset output_bcf rest; do
  [ -n "$label" ] || continue
  require_file "$source"
  out_dir=$(dirname "$output_bcf")
  tmp_bcf="$TMPDIR_ROOT/${label}.munged.bcf"
  mkdir -p "$out_dir"

  munge_args=(+munge -f "$source_fasta" -s "$sample_name" -Ob -o "$tmp_bcf")
  if [ "$columns_file" != "." ]; then
    require_file "$columns_file"
    munge_args=(+munge -C "$columns_file" -f "$source_fasta" -s "$sample_name" -Ob -o "$tmp_bcf")
  elif [ "$columns_preset" != "." ]; then
    munge_args=(+munge -c "$columns_preset" -f "$source_fasta" -s "$sample_name" -Ob -o "$tmp_bcf")
  fi

  msg "munge $label"
  run_cmd "$BCFTOOLS" "${munge_args[@]}" "$source"
  run_cmd "$BCFTOOLS" index --threads "$THREADS" -f "$tmp_bcf"

  if [ "$source_build" = "GRCh38" ] || [ "$chain" = "." ]; then
    msg "sort/index $label without liftover"
    run_bash "'$BCFTOOLS' sort -Ob -o '$output_bcf' '$tmp_bcf' && '$BCFTOOLS' index --threads '$THREADS' -f '$output_bcf'"
  else
    require_file "$source_fasta"
    require_file "$target_fasta"
    require_file "$chain"
    msg "liftover $label to GRCh38"
    run_bash "'$BCFTOOLS' +liftover -Ou '$tmp_bcf' -- -s '$source_fasta' -f '$target_fasta' -c '$chain' | '$BCFTOOLS' sort -Ob -o '$output_bcf' && '$BCFTOOLS' index --threads '$THREADS' -f '$output_bcf'"
  fi
done
