#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: prs_base_from_gwas_bcf.sh --bcf FILE --out FILE.gz [--sample NAME]

Create and validate a canonical CHR/BP/SNP/A1/A2/BETA/P base table from a
single-sample GWAS-VCF BCF. ES is interpreted as the ALT-allele beta and LP as
-log10(P). Only PASS autosomal biallelic SNPs are retained.
USAGE
}

BCF_FILE=
OUT=
SAMPLE=
BCFTOOLS=${BCFTOOLS:-/opt/sw/bin/bcftools}
BGZIP=${BGZIP:-bgzip}

while (( $# )); do
  case "$1" in
    --bcf) BCF_FILE=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --sample) SAMPLE=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$BCF_FILE" && -n "$OUT" ]] || { usage >&2; exit 2; }
[[ -s "$BCF_FILE" ]] || { printf 'Missing BCF: %s\n' "$BCF_FILE" >&2; exit 1; }

mapfile -t SAMPLES < <("$BCFTOOLS" query -l "$BCF_FILE")
if [[ -z "$SAMPLE" ]]; then
  [[ ${#SAMPLES[@]} -eq 1 ]] || {
    printf 'Expected one GWAS sample, found %d; use --sample\n' "${#SAMPLES[@]}" >&2
    exit 1
  }
  SAMPLE=${SAMPLES[0]}
fi

mkdir -p "$(dirname "$OUT")"
RAW=$(mktemp "$(dirname "$OUT")/.prs-base.raw.XXXXXX.tsv.gz")
trap 'rm -f "$RAW"' EXIT
REGIONS=$(printf 'chr%s,' {1..22})
REGIONS=${REGIONS%,}

## convert LP without allowing floating-point underflow to produce P=0
{
  printf 'CHR\tBP\tSNP\tA1\tA2\tBETA\tP\n'
  "$BCFTOOLS" query \
    -r "$REGIONS" \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%ES\t%LP]\n' \
    -s "$SAMPLE" \
    -i '(FILTER="PASS" || FILTER=".") && N_ALT=1 && TYPE="snp"' \
    "$BCF_FILE" |
    awk -F '\t' 'BEGIN{OFS="\t"}
      NF==6 && $5!="." && $6!="." {
        lp=$6+0
        p=(lp>=300 ? 1e-300 : exp(-lp*log(10)))
        print $1,$2,$1":"$2":"$3":"$4,$4,$3,$5,p
      }'
} | "$BGZIP" -c > "$RAW"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
perl "$SCRIPT_DIR/validate_prs_base.pl" \
  --input "$RAW" \
  --output "$OUT" \
  --stats "$OUT.validation.tsv"

printf 'Wrote %s\n' "$OUT"
