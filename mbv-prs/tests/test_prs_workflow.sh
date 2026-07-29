#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FIXTURES="$ROOT/tests/fixtures/prs"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

perl "$ROOT/scripts/validate_prs_base.pl" \
  --input "$FIXTURES/base.valid.tsv" \
  --output "$TMP/base.valid.tsv.gz" \
  --stats "$TMP/base.validation.tsv"

rows=$(gzip -cd "$TMP/base.valid.tsv.gz" | wc -l)
[[ "$rows" -eq 4 ]] || { printf 'Unexpected validated row count: %s\n' "$rows" >&2; exit 1; }

for invalid in duplicate nonauto allele_mismatch; do
  if perl "$ROOT/scripts/validate_prs_base.pl" \
    --input "$FIXTURES/base.$invalid.tsv" \
    --output "$TMP/base.$invalid.tsv.gz" \
    --strict > "$TMP/$invalid.stdout" 2> "$TMP/$invalid.stderr"; then
    printf 'Validator accepted invalid fixture: %s\n' "$invalid" >&2
    exit 1
  fi
done

## non-strict validation must remove every occurrence of a duplicate ID
perl "$ROOT/scripts/validate_prs_base.pl" \
  --input "$FIXTURES/base.duplicate.tsv" \
  --output "$TMP/base.duplicate.filtered.tsv.gz" \
  --stats "$TMP/base.duplicate.validation.tsv"

duplicate_ids=$(awk -F '\t' '$1=="duplicate_ids"{print $2}' "$TMP/base.duplicate.validation.tsv")
duplicate_rows=$(awk -F '\t' '$1=="rejected_duplicate_rows"{print $2}' "$TMP/base.duplicate.validation.tsv")
accepted_rows=$(awk -F '\t' '$1=="accepted_rows"{print $2}' "$TMP/base.duplicate.validation.tsv")
[[ "$duplicate_ids" -eq 1 && "$duplicate_rows" -eq 2 && "$accepted_rows" -eq 0 ]] || {
  printf 'Unexpected duplicate counts: ids=%s rows=%s accepted=%s\n' \
    "$duplicate_ids" "$duplicate_rows" "$accepted_rows" >&2
  exit 1
}

perl "$ROOT/scripts/prepare_ct_inputs.pl" \
  --base "$TMP/base.valid.tsv.gz" \
  --target-pvar "$FIXTURES/target.pvar" \
  --ld-pattern "$FIXTURES/ld.chr{CHR}" \
  --out-dir "$TMP/mapping" \
  --chromosomes 1

matched=$(awk -F '\t' '$1=="matched"{print $2}' "$TMP/mapping/chr1.mapping-stats.tsv")
target_mismatch=$(awk -F '\t' '$1=="target_allele_mismatch"{print $2}' "$TMP/mapping/chr1.mapping-stats.tsv")
missing_ld=$(awk -F '\t' '$1=="missing_ld"{print $2}' "$TMP/mapping/chr1.mapping-stats.tsv")
[[ "$matched" -eq 1 && "$target_mismatch" -eq 1 && "$missing_ld" -eq 1 ]] || {
  printf 'Unexpected mapping counts: matched=%s target_mismatch=%s missing_ld=%s\n' \
    "$matched" "$target_mismatch" "$missing_ld" >&2
  exit 1
}

selected=$(perl "$ROOT/scripts/extract_ct_scores.pl" \
  --clumped "$FIXTURES/chr1.clumped" \
  --lookup "$TMP/mapping/chr1.lookup.tsv" \
  --output "$TMP/selected.tsv")
[[ "$selected" -eq 1 ]] || { printf 'Unexpected selected count: %s\n' "$selected" >&2; exit 1; }
grep -q $'^target_100\tG\t0.10\t0.01$' "$TMP/selected.tsv"

printf 'PRS workflow fixtures passed\n'
