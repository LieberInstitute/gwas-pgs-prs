#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_prsice_manifest.sh --manifest FILE [--threads N] [--resume]

Manifest columns:
RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_BGEN LD_PREFIX PHENO OUT_PREFIX
CLUMP_KB CLUMP_R2 THRESHOLDS
USAGE
}

MANIFEST=
THREADS=${THREADS:-4}
RESUME=0
PRSICE=${PRSICE:-/opt/sw/bin/PRSice_linux}

while (( $# )); do
  case "$1" in
    --manifest) MANIFEST=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -s "$MANIFEST" ]] || { printf 'Missing manifest: %s\n' "$MANIFEST" >&2; exit 1; }
EXPECTED_HEADER=$'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tBASE\tTARGET_BGEN\tLD_PREFIX\tPHENO\tOUT_PREFIX\tCLUMP_KB\tCLUMP_R2\tTHRESHOLDS'
IFS= read -r HEADER < "$MANIFEST"
[[ "$HEADER" == "$EXPECTED_HEADER" ]] || {
  printf 'Unexpected PRSice manifest header\nExpected: %s\nObserved: %s\n' "$EXPECTED_HEADER" "$HEADER" >&2
  exit 1
}

while IFS=$'\t' read -r RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_BGEN LD_PREFIX PHENO OUT_PREFIX CLUMP_KB CLUMP_R2 THRESHOLDS; do
  [[ "$RUN_ID" == "RUN_ID" ]] && continue
  [[ -n "$RUN_ID" ]] || continue

  [[ -s "$BASE" ]] || { printf '[%s] missing base: %s\n' "$RUN_ID" "$BASE" >&2; exit 1; }
  [[ -s "$TARGET_BGEN.bgen" && -s "$TARGET_BGEN.sample" ]] || {
    printf '[%s] incomplete BGEN target: %s\n' "$RUN_ID" "$TARGET_BGEN" >&2
    exit 1
  }
  [[ -s "$LD_PREFIX.bed" && -s "$LD_PREFIX.bim" && -s "$LD_PREFIX.fam" ]] || {
    printf '[%s] incomplete PLINK LD reference: %s\n' "$RUN_ID" "$LD_PREFIX" >&2
    exit 1
  }
  [[ -s "$PHENO" ]] || { printf '[%s] missing phenotype file: %s\n' "$RUN_ID" "$PHENO" >&2; exit 1; }

  mkdir -p "$(dirname "$OUT_PREFIX")"
  if (( RESUME )) && [[ -s "$OUT_PREFIX.all_score" && -s "$OUT_PREFIX.prsice" ]]; then
    printf '[%s] reused completed PRSice output\n' "$RUN_ID"
    continue
  fi

  printf '[%s] starting PRSice run\n' "$RUN_ID"
  "$PRSICE" \
    --base "$BASE" \
    --target "$TARGET_BGEN" \
    --type bgen \
    --ld "$LD_PREFIX" \
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
    --pheno "$PHENO" \
    --pheno-col DUMMY \
    --no-regress \
    --all-score \
    --score sum \
    --fastscore \
    --bar-levels "$THRESHOLDS" \
    --clump-kb "$CLUMP_KB" \
    --clump-r2 "$CLUMP_R2" \
    --clump-p 1 \
    --thread "$THREADS" \
    --out "$OUT_PREFIX" > "$OUT_PREFIX.stdout.log" 2>&1
  printf 'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\n%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "$TRAIT" "$GWAS_VERSION" "$SETUP" > "$OUT_PREFIX.run-identity.tsv"
  printf '[%s] completed PRSice run\n' "$RUN_ID"
done < "$MANIFEST"
