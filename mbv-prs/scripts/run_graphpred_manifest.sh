#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_graphpred_manifest.sh --manifest FILE [--regions LIST] [--resume]

Manifest columns:
RUN_ID TRAIT GWAS_VERSION SETUP TARGET_BCF PGS_BCF OUT_FILE
USAGE
}

MANIFEST=
REGIONS=
RESUME=0
BCFTOOLS=${BCFTOOLS:-/opt/sw/bin/bcftools}

while (( $# )); do
  case "$1" in
    --manifest) MANIFEST=$2; shift 2 ;;
    --regions) REGIONS=$2; shift 2 ;;
    --resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -s "$MANIFEST" ]] || { printf 'Missing manifest: %s\n' "$MANIFEST" >&2; exit 1; }
EXPECTED_HEADER=$'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tTARGET_BCF\tPGS_BCF\tOUT_FILE'
IFS= read -r HEADER < "$MANIFEST"
[[ "$HEADER" == "$EXPECTED_HEADER" ]] || {
  printf 'Unexpected GraphPred manifest header\nExpected: %s\nObserved: %s\n' "$EXPECTED_HEADER" "$HEADER" >&2
  exit 1
}

while IFS=$'\t' read -r RUN_ID TRAIT GWAS_VERSION SETUP TARGET_BCF PGS_BCF OUT_FILE; do
  [[ "$RUN_ID" == "RUN_ID" ]] && continue
  [[ -n "$RUN_ID" ]] || continue
  [[ -s "$TARGET_BCF" && -s "$PGS_BCF" ]] || {
    printf '[%s] missing target or PGS BCF\n' "$RUN_ID" >&2
    exit 1
  }

  mkdir -p "$(dirname "$OUT_FILE")"
  if (( RESUME )) && [[ -s "$OUT_FILE" ]]; then
    printf '[%s] reused completed GraphPred output\n' "$RUN_ID"
    continue
  fi

  REGION_ARGS=()
  [[ -n "$REGIONS" ]] && REGION_ARGS=(-r "$REGIONS")
  printf '[%s] starting GraphPred score\n' "$RUN_ID"
  "$BCFTOOLS" +score \
    --use DS \
    --sample-header SAMPLE \
    --counts \
    "${REGION_ARGS[@]}" \
    -o "$OUT_FILE" \
    "$TARGET_BCF" "$PGS_BCF" 2> "$OUT_FILE.log"
  printf 'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\n%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "$TRAIT" "$GWAS_VERSION" "$SETUP" > "$OUT_FILE.run-identity.tsv"
  printf '[%s] completed GraphPred score\n' "$RUN_ID"
done < "$MANIFEST"
