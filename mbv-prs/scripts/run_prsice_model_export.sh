#!/usr/bin/env bash
set -euo pipefail

## resolve paths from the repository root, not the caller cwd
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

PRSICE=${PRSICE:-"$ROOT/tools/prsice-2.3.5/PRSice_linux"}
LDREF=${LDREF:-/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos}
THREADS=${THREADS:-8}
EXPORT_DIR=${EXPORT_DIR:-"$ROOT/prsice_model_export"}
THRESHOLDS=${THRESHOLDS:-5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1}

mkdir -p "$EXPORT_DIR/out" "$EXPORT_DIR/logs"
: > "$EXPORT_DIR/logs/commands_successful.tsv"

run_cmd() {
  local name=$1
  shift
  printf '[%s] %s\n' "$(date '+%F %T')" "$name"
  "$@"
  printf '%s\t%s\n' "$(date '+%F %T')" "$name" >> "$EXPORT_DIR/logs/commands_successful.tsv"
}

run_bash() {
  local name=$1
  local cmd=$2
  printf '[%s] %s\n' "$(date '+%F %T')" "$name"
  bash -o pipefail -c "$cmd"
  printf '%s\t%s\n' "$(date '+%F %T')" "$name" >> "$EXPORT_DIR/logs/commands_successful.tsv"
}

## fail early if required inputs are missing
for path in \
  "$PRSICE" \
  "$LDREF.bed" "$LDREF.bim" "$LDREF.fam" \
  prsice/target/mbv_qc.bgen prsice/target/mbv_qc.sample \
  prsice/target/dummy.pheno \
  prsice/base/BPD.base.tsv.gz prsice/base/MDD.base.tsv.gz
do
  if [ ! -e "$path" ]; then
    printf 'missing required input: %s\n' "$path" >&2
    exit 1
  fi
done

run_bash "record PRSice model-export inputs" \
  "sha256sum '$PRSICE' prsice/base/BPD.base.tsv.gz prsice/base/MDD.base.tsv.gz prsice/target/mbv_qc.bgen prsice/target/mbv_qc.sample > '$EXPORT_DIR/logs/input_hashes.sha256'"

{
  printf 'PRSice:\n'
  "$PRSICE" --version 2>&1
  printf 'PRSice help print-snp:\n'
  "$PRSICE" --help 2>&1 | grep -A6 -- '--print-snp'
} > "$EXPORT_DIR/logs/tool_versions.txt"

run_bash "write PRSice model-export command manifest header" \
  "printf 'disorder\tcommand\n' > '$EXPORT_DIR/logs/prsice_commands.tsv'"

for label in BPD MDD; do
  out_prefix="$EXPORT_DIR/out/$label"
  cmd=(
    "$PRSICE"
    --base "prsice/base/${label}.base.tsv.gz"
    --target prsice/target/mbv_qc
    --type bgen
    --ld "$LDREF"
    --ld-type bed
    --snp SNP
    --chr CHR
    --bp BP
    --a1 A1
    --a2 A2
    --stat BETA
    --pvalue P
    --beta
    --ignore-fid
    --pheno prsice/target/dummy.pheno
    --pheno-col DUMMY
    --no-regress
    --all-score
    --score sum
    --fastscore
    --bar-levels "$THRESHOLDS"
    --clump-kb 250
    --clump-r2 0.1
    --clump-p 1
    --print-snp
    --thread "$THREADS"
    --out "$out_prefix"
  )

  printf '%s\t' "$label" >> "$EXPORT_DIR/logs/prsice_commands.tsv"
  printf '%q ' "${cmd[@]}" >> "$EXPORT_DIR/logs/prsice_commands.tsv"
  printf '\n' >> "$EXPORT_DIR/logs/prsice_commands.tsv"

  run_cmd "PRSice ${label} with print-snp" "${cmd[@]}"
done

run_cmd "export PRSice model tables" \
  Rscript scripts/export_prsice_model_tables.R "$EXPORT_DIR"

printf '[%s] done\n' "$(date '+%F %T')"
