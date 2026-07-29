#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: merge_plink_chromosomes.sh --pattern PREFIX_{CHR} --out PREFIX [options]

Options:
  --chromosomes LIST  comma/range list, default 1-22
  --plink FILE        PLINK 1 executable

The script verifies identical FAM files, records input/output checksums and
variant counts, and creates one PLINK BED prefix suitable for PRSice --ld.
USAGE
}

PATTERN=
OUT=
CHROMOSOMES=1-22
PLINK=${PLINK:-/opt/sw/bin/plink}

while (( $# )); do
  case "$1" in
    --pattern) PATTERN=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --chromosomes) CHROMOSOMES=$2; shift 2 ;;
    --plink) PLINK=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$PATTERN" && -n "$OUT" ]] || { usage >&2; exit 2; }
[[ "$PATTERN" == *'{CHR}'* ]] || { printf 'Pattern must contain {CHR}\n' >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
WORK="$(dirname "$OUT")/provenance"
mkdir -p "$WORK"

mapfile -t CHRS < <(perl -e '
  my %seen;
  for my $token (split /,/, $ARGV[0]) {
    my @x = $token =~ /^(\d+)-(\d+)$/ ? ($1 <= $2 ? $1..$2 : reverse $2..$1)
          : $token =~ /^\d+$/ ? ($token) : die "Invalid chromosome: $token\n";
    for my $chr (@x) {
      die "Non-autosomal chromosome: $chr\n" if $chr < 1 || $chr > 22;
      print "$chr\n" unless $seen{$chr}++;
    }
  }' "$CHROMOSOMES")

: > "$WORK/input-files.txt"
: > "$WORK/merge-list.txt"
: > "$WORK/variant-counts.tsv"
printf 'chromosome\tvariants\n' > "$WORK/variant-counts.tsv"
first_prefix=
first_fam=
total=0

for chr in "${CHRS[@]}"; do
  prefix=${PATTERN//\{CHR\}/$chr}
  for ext in bed bim fam; do
    [[ -s "$prefix.$ext" ]] || { printf 'Missing %s.%s\n' "$prefix" "$ext" >&2; exit 1; }
    printf '%s\n' "$prefix.$ext" >> "$WORK/input-files.txt"
  done
  count=$(wc -l < "$prefix.bim")
  printf '%s\t%s\n' "$chr" "$count" >> "$WORK/variant-counts.tsv"
  total=$((total + count))

  if [[ -z "$first_prefix" ]]; then
    first_prefix=$prefix
    first_fam=$prefix.fam
  else
    cmp -s "$first_fam" "$prefix.fam" || {
      printf 'FAM differs from chromosome 1: %s.fam\n' "$prefix" >&2
      exit 1
    }
    printf '%s\n' "$prefix" >> "$WORK/merge-list.txt"
  fi
done

: > "$WORK/input-checksums.sha256"
while IFS= read -r file; do
  sha256sum "$file" >> "$WORK/input-checksums.sha256"
done < "$WORK/input-files.txt"
"$PLINK" --bfile "$first_prefix" --merge-list "$WORK/merge-list.txt" \
  --make-bed --out "$OUT" > "$WORK/merge.stdout.log" 2>&1

merged=$(wc -l < "$OUT.bim")
samples=$(wc -l < "$OUT.fam")
printf 'total_input\t%s\nmerged\t%s\nsamples\t%s\n' "$total" "$merged" "$samples" \
  >> "$WORK/variant-counts.tsv"
[[ "$merged" -eq "$total" ]] || {
  printf 'Merged variant count differs: %s != %s\n' "$merged" "$total" >&2
  exit 1
}
sha256sum "$OUT.bed" "$OUT.bim" "$OUT.fam" > "$WORK/output-checksums.sha256"
printf 'Wrote %s with %s variants and %s samples\n' "$OUT" "$merged" "$samples"
