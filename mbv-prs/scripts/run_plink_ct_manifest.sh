#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_plink_ct_manifest.sh --manifest FILE [options]

Options:
  --chromosomes LIST  comma/range list, default 1-22
  --threads N         PLINK2 scoring threads, default 4
  --resume            reuse completed chromosome and score outputs

Manifest columns:
RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_PGEN LD_PATTERN OUT_DIR
CLUMP_KB CLUMP_R2 THRESHOLDS

LD_PATTERN must contain {CHR}. THRESHOLDS is a comma-separated list.
USAGE
}

MANIFEST=
CHROMOSOMES=1-22
THREADS=${THREADS:-4}
RESUME=0
PLINK=${PLINK:-/opt/sw/bin/plink}
PLINK2=${PLINK2:-/opt/sw/bin/plink2}

while (( $# )); do
  case "$1" in
    --manifest) MANIFEST=$2; shift 2 ;;
    --chromosomes) CHROMOSOMES=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -s "$MANIFEST" ]] || { printf 'Missing manifest: %s\n' "$MANIFEST" >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

EXPECTED_HEADER=$'RUN_ID\tTRAIT\tGWAS_VERSION\tSETUP\tBASE\tTARGET_PGEN\tLD_PATTERN\tOUT_DIR\tCLUMP_KB\tCLUMP_R2\tTHRESHOLDS'
IFS= read -r HEADER < "$MANIFEST"
[[ "$HEADER" == "$EXPECTED_HEADER" ]] || {
  printf 'Unexpected C+T manifest header\nExpected: %s\nObserved: %s\n' "$EXPECTED_HEADER" "$HEADER" >&2
  exit 1
}

expand_chromosomes() {
  perl -e '
    use strict; use warnings;
    my %seen;
    for my $token (split /,/, $ARGV[0]) {
      my @values = $token =~ /^(\d+)-(\d+)$/ ? ($1 <= $2 ? $1..$2 : reverse $2..$1)
                 : $token =~ /^\d+$/ ? ($token) : die "Invalid chromosome: $token\n";
      for my $chr (@values) {
        die "Non-autosomal chromosome: $chr\n" if $chr < 1 || $chr > 22;
        print "$chr\n" unless $seen{$chr}++;
      }
    }' "$1"
}
mapfile -t CHRS < <(expand_chromosomes "$CHROMOSOMES")

run_one() {
  local run_id=$1 trait=$2 version=$3 setup=$4 base=$5 target=$6
  local ld_pattern=$7 out_dir=$8 clump_kb=$9 clump_r2=${10} thresholds=${11}

  [[ -s "$base" ]] || { printf '[%s] missing base: %s\n' "$run_id" "$base" >&2; return 1; }
  [[ -s "$target.pgen" && -s "$target.pvar" && -s "$target.psam" ]] || {
    printf '[%s] incomplete target PGEN prefix: %s\n' "$run_id" "$target" >&2
    return 1
  }
  [[ "$ld_pattern" == *'{CHR}'* ]] || {
    printf '[%s] LD_PATTERN lacks {CHR}: %s\n' "$run_id" "$ld_pattern" >&2
    return 1
  }

  mkdir -p "$out_dir" "$out_dir/clump" "$out_dir/selected" "$out_dir/score" "$out_dir/logs"
  printf 'chromosome\tstatus\tselected_variants\n' > "$out_dir/chromosome-status.tsv"

  ## map canonical base variants to both LD-reference and target IDs
  perl "$SCRIPT_DIR/prepare_ct_inputs.pl" \
    --base "$base" \
    --target-pvar "$target.pvar" \
    --ld-pattern "$ld_pattern" \
    --out-dir "$out_dir/mapping" \
    --chromosomes "$CHROMOSOMES"

  for chr in "${CHRS[@]}"; do
    local ref=${ld_pattern//\{CHR\}/$chr}
    local prefix="$out_dir/clump/chr$chr"
    local selected="$out_dir/selected/chr$chr.tsv"

    if (( RESUME )) && [[ -s "$prefix.clumped" && -s "$selected" ]]; then
      local selected_n
      selected_n=$(awk 'END{print NR>0 ? NR-1 : 0}' "$selected")
      printf '%s\treused\t%s\n' "$chr" "$selected_n" >> "$out_dir/chromosome-status.tsv"
      continue
    fi

    if ! "$PLINK" \
      --bfile "$ref" \
      --clump "$out_dir/mapping/chr$chr.clump.tsv" \
      --clump-p1 1 \
      --clump-p2 1 \
      --clump-r2 "$clump_r2" \
      --clump-kb "$clump_kb" \
      --out "$prefix" > "$out_dir/logs/chr$chr.stdout.log" 2>&1; then
      printf '%s\tfailed\t0\n' "$chr" >> "$out_dir/chromosome-status.tsv"
      printf '[%s] PLINK clumping failed on chromosome %s\n' "$run_id" "$chr" >&2
      return 1
    fi

    local selected_n
    if ! selected_n=$(perl "$SCRIPT_DIR/extract_ct_scores.pl" \
      --clumped "$prefix.clumped" \
      --lookup "$out_dir/mapping/chr$chr.lookup.tsv" \
      --output "$selected"); then
      printf '%s\tfailed_extract\t0\n' "$chr" >> "$out_dir/chromosome-status.tsv"
      return 1
    fi
    printf '%s\tcompleted\t%s\n' "$chr" "$selected_n" >> "$out_dir/chromosome-status.tsv"
  done

  ## combine selected variants and retain legacy q-score labels exactly
  printf 'SNP\tA1\tBETA\n' > "$out_dir/score/profile.tsv"
  printf 'SNP\tP\n' > "$out_dir/score/profile_p.tsv"
  for chr in "${CHRS[@]}"; do
    awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3}' "$out_dir/selected/chr$chr.tsv" >> "$out_dir/score/profile.tsv"
    awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$4}' "$out_dir/selected/chr$chr.tsv" >> "$out_dir/score/profile_p.tsv"
  done

  : > "$out_dir/score/ranges.tsv"
  IFS=',' read -r -a LEVELS <<< "$thresholds"
  for level in "${LEVELS[@]}"; do
    printf '%s\t0\t%s\n' "$level" "$level" >> "$out_dir/score/ranges.tsv"
  done

  if (( RESUME )); then
    local complete=1
    for level in "${LEVELS[@]}"; do
      [[ -s "$out_dir/score/score.$level.sscore" ]] || complete=0
    done
    if (( complete )); then
      printf '[%s] reused completed scores\n' "$run_id"
      return 0
    fi
  fi

  "$PLINK2" \
    --pfile "$target" \
    --score "$out_dir/score/profile.tsv" 1 2 3 header ignore-dup-ids \
    --q-score-range "$out_dir/score/ranges.tsv" "$out_dir/score/profile_p.tsv" header \
    --threads "$THREADS" \
    --out "$out_dir/score/score" > "$out_dir/logs/scoring.stdout.log" 2>&1

  printf '%s\t%s\t%s\t%s\n' "$run_id" "$trait" "$version" "$setup" > "$out_dir/run-identity.tsv"
}

while IFS=$'\t' read -r RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_PGEN LD_PATTERN OUT_DIR CLUMP_KB CLUMP_R2 THRESHOLDS; do
  [[ "$RUN_ID" == "RUN_ID" ]] && continue
  [[ -n "$RUN_ID" ]] || continue
  printf '[%s] starting C+T run\n' "$RUN_ID"
  run_one "$RUN_ID" "$TRAIT" "$GWAS_VERSION" "$SETUP" "$BASE" "$TARGET_PGEN" \
    "$LD_PATTERN" "$OUT_DIR" "$CLUMP_KB" "$CLUMP_R2" "$THRESHOLDS"
  printf '[%s] completed C+T run\n' "$RUN_ID"
done < "$MANIFEST"
