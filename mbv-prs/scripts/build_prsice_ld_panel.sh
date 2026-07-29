#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: build_prsice_ld_panel.sh --source-prefix PREFIX --ld-pattern PREFIX_{CHR}
       --target-pvar FILE --out-prefix PREFIX

Create a target-overlapping copy of a PLINK LD panel whose variant IDs match
canonical target IDs. Genotypes and samples are unchanged.
USAGE
}

SOURCE=
LD_PATTERN=
TARGET_PVAR=
OUT=
PLINK=${PLINK:-/opt/sw/bin/plink}

while (( $# )); do
  case "$1" in
    --source-prefix) SOURCE=$2; shift 2 ;;
    --ld-pattern) LD_PATTERN=$2; shift 2 ;;
    --target-pvar) TARGET_PVAR=$2; shift 2 ;;
    --out-prefix) OUT=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SOURCE" && -n "$LD_PATTERN" && -n "$TARGET_PVAR" && -n "$OUT" ]] || {
  usage >&2
  exit 2
}
for ext in bed bim fam; do
  [[ -s "$SOURCE.$ext" ]] || { printf 'Missing source %s.%s\n' "$SOURCE" "$ext" >&2; exit 1; }
done
[[ -s "$TARGET_PVAR" ]] || { printf 'Missing target PVAR: %s\n' "$TARGET_PVAR" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK="$(dirname "$OUT")/prsice-canonical-provenance"
mkdir -p "$WORK"
UPDATE="$WORK/ld_to_target_ids.tsv"

perl "$SCRIPT_DIR/map_plink_ld_to_target.pl" \
  --ld-pattern "$LD_PATTERN" \
  --target-pvar "$TARGET_PVAR" \
  --output "$UPDATE" \
  --stats "$WORK/ld_to_target_ids.stats.tsv" \
  --work-dir "$WORK/target-split"
## PLINK applies --update-name before --extract, so extract by the new IDs
cut -f2 "$UPDATE" > "$WORK/extract.ids"

"$PLINK" --bfile "$SOURCE" \
  --extract "$WORK/extract.ids" \
  --update-name "$UPDATE" 2 1 \
  --make-bed \
  --out "$OUT" > "$WORK/build.stdout.log" 2>&1

mapped=$(wc -l < "$UPDATE")
variants=$(wc -l < "$OUT.bim")
samples=$(wc -l < "$OUT.fam")
[[ "$mapped" -eq "$variants" ]] || {
  printf 'Canonical LD count mismatch: map=%s panel=%s\n' "$mapped" "$variants" >&2
  exit 1
}
awk '$2 !~ /^chr([1-9]|1[0-9]|2[0-2]):[0-9]+:[ACGT]+:[ACGT]+$/{bad++} END{exit bad>0}' \
  "$OUT.bim" || { printf 'Noncanonical IDs remain in %s.bim\n' "$OUT" >&2; exit 1; }

printf 'metric\tvalue\nvariants\t%s\nsamples\t%s\n' "$variants" "$samples" \
  > "$WORK/canonical-panel-counts.tsv"
sha256sum "$SOURCE.bed" "$SOURCE.bim" "$SOURCE.fam" "$TARGET_PVAR" \
  > "$WORK/canonical-panel-input-checksums.sha256"
sha256sum "$OUT.bed" "$OUT.bim" "$OUT.fam" \
  > "$WORK/canonical-panel-output-checksums.sha256"
printf 'Wrote %s with %s variants and %s samples\n' "$OUT" "$variants" "$samples"
