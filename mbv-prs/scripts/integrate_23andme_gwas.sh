#!/usr/bin/env bash
set -euo pipefail

usage() {
    sed -n '2,95p' "$0" | sed -n 's/^## //p'
}

## Integrate the public no-23andMe and 23andMe-only European GWAS components.
##
## Usage:
##   scripts/integrate_23andme_gwas.sh --trait BD|MDD --out-prefix PATH \
##     --confirm-no-sample-overlap [options]
##
## Required:
##   --trait BD|MDD                 study profile
##   --out-prefix PATH              output path prefix
##   --confirm-no-sample-overlap    assert that the two components are independent
##
## Important options:
##   --public FILE                  override public no-23andMe summary statistics
##   --23me-assoc FILE              override 23andMe association file
##   --23me-annotation FILE         override 23andMe all_snp_info file
##   --annotation-position-offset N add N to supplied annotation positions
##   --min-ne-fraction FLOAT        override paper effective-N threshold
##   --allow-bd-v7.2-annotations    reproduce a noncanonical BD v7.2 mapping run
##   --exclude-list FILE            canonical GRCh37 CHROM POS list to exclude
##   --threads INT                  worker threads [4]
##   --resume                       reuse existing nonempty stage outputs
##   --skip-sha256                  omit input SHA-256 hashes from the report
##   --help                         show this help

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

require_file() {
    [[ -s "$1" ]] || die "missing or empty file: $1"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

stream_file() {
    case "$1" in
        *.gz) gzip -cd -- "$1" ;;
        *) command cat -- "$1" ;;
    esac
}

check_output() {
    local output=$1
    if [[ -e "$output" && "$resume" -eq 0 ]]; then
        die "output exists; use --resume only after checking it: $output"
    fi
}

trait=
out_prefix=
public=
assoc=
annotation=
annotation_was_overridden=0
annotation_position_offset=
min_ne_fraction=
confirm_independent=0
allow_bd_annotation_mismatch=0
exclude_list=
threads=4
resume=0
skip_sha256=0

while (($#)); do
    case "$1" in
        --trait) trait=${2:-}; shift 2 ;;
        --out-prefix) out_prefix=${2:-}; shift 2 ;;
        --public) public=${2:-}; shift 2 ;;
        --23me-assoc) assoc=${2:-}; shift 2 ;;
        --23me-annotation) annotation=${2:-}; annotation_was_overridden=1; shift 2 ;;
        --annotation-position-offset) annotation_position_offset=${2:-}; shift 2 ;;
        --min-ne-fraction) min_ne_fraction=${2:-}; shift 2 ;;
        --confirm-no-sample-overlap) confirm_independent=1; shift ;;
        --allow-bd-v7.2-annotations) allow_bd_annotation_mismatch=1; shift ;;
        --exclude-list) exclude_list=${2:-}; shift 2 ;;
        --threads) threads=${2:-}; shift 2 ;;
        --resume) resume=1; shift ;;
        --skip-sha256) skip_sha256=1; shift ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

[[ "$trait" == "BD" || "$trait" == "MDD" ]] || die "--trait must be BD or MDD"
[[ -n "$out_prefix" ]] || die "--out-prefix is required"
((confirm_independent == 1)) || die "--confirm-no-sample-overlap is required"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"

root=/home/gpertea/work/ref/GWAS
delivery="$root/23andMe_MDD_BD"
ref37=/dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta
ref38=/dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
chain=/dbdata/cdb/ref/hg19ToHg38.over.chain.gz
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
columns="$script_dir/gwas_meta_colheaders.tsv"
prepare_awk="$script_dir/prepare_23andme_gwas.awk"
extra_headers="$script_dir/gwas_meta_extra_headers.txt"

if [[ "$trait" == "BD" ]]; then
    public=${public:-$root/BD/bip2024_eur_no23andMe.gz}
    assoc=${assoc:-$delivery/Bipolar-Disorder-O_Connell-2025/OConnell_2025_bipolar_european-7.0/bipolar.dat.gz}
    bd_v70_annotation=$delivery/7.0-Annotations/v7.0_europe/all_snp_info.txt.gz
    annotation=${annotation:-$bd_v70_annotation}
    min_ne_fraction=${min_ne_fraction:-0.75}
    final_sample=BD_2024_FULL_EUR
    expected_max_ne=440934
    if [[ "$annotation" == "$bd_v70_annotation" ]]; then
        annotation_position_offset=${annotation_position_offset:--1}
        annotation_release_note="BD used matching v7.0 association and annotation releases."
    elif [[ "$annotation" == "$delivery/7.2-Annotations/v7.2_europe/all_snp_info.txt.gz" ]]; then
        annotation_position_offset=${annotation_position_offset:-0}
        ((allow_bd_annotation_mismatch == 1)) || die \
            "BD association is v7.0 but annotation is v7.2; use the v7.0 annotation or pass --allow-bd-v7.2-annotations"
        annotation_release_note="BD used noncanonical v7.2 annotations with the v7.0 association file."
    else
        ((annotation_was_overridden == 1)) || die "internal error resolving BD annotation"
        [[ -n "$annotation_position_offset" ]] || die \
            "a custom BD annotation requires --annotation-position-offset"
        annotation_release_note="BD used a custom annotation with an operator-specified position offset."
    fi
else
    public=${public:-$root/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz}
    assoc=${assoc:-$delivery/MDD-Adams-2025/Adams_2025_mdd_european-7.2/mdd.dat.gz}
    annotation=${annotation:-$delivery/7.2-Annotations/v7.2_europe/all_snp_info.txt.gz}
    min_ne_fraction=${min_ne_fraction:-0.80}
    final_sample=MDD_2025_FULL_EUR
    expected_max_ne=1577206
    annotation_position_offset=${annotation_position_offset:-0}
    annotation_release_note="MDD used matching v7.2 association and annotation releases."
fi

[[ "$annotation_position_offset" =~ ^-?[0-9]+$ ]] || \
    die "--annotation-position-offset must be an integer"
awk -v value="$min_ne_fraction" 'BEGIN{exit !(value > 0 && value <= 1)}' || \
    die "--min-ne-fraction must be in (0,1]"

for command_name in awk bcftools bgzip gzip paste sha256sum tabix; do
    require_command "$command_name"
done
for input in "$public" "$assoc" "$annotation" "$ref37" "$ref38" "$chain" \
    "$columns" "$prepare_awk" "$extra_headers"; do
    require_file "$input"
done
[[ -z "$exclude_list" ]] || require_file "$exclude_list"

export BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-/opt/sw/bcf-plugins}
available_plugins=$(bcftools plugin -l)
for plugin_name in munge metal liftover; do
    grep -qx "$plugin_name" <<< "$available_plugins" || \
        die "bcftools +$plugin_name is unavailable in $BCFTOOLS_PLUGINS"
done

mkdir -p "$(dirname "$out_prefix")"
prepared="$out_prefix.23andme.grch37.ssf.tsv.gz"
public_bcf="$out_prefix.public-no23.grch37.bcf"
component_bcf="$out_prefix.23andme.grch37.bcf"
unfiltered_bcf="$out_prefix.meta-unfiltered.grch37.bcf"
filtered_bcf="$out_prefix.full.grch37.bcf"
final_bcf="$out_prefix.full.grch38.bcf"
meta_tsv="$out_prefix.full.grch38.meta.tsv.gz"
prsice_tsv="$out_prefix.full.grch38.prsice.tsv.gz"
report="$out_prefix.integration-report.md"

for output in "$prepared" "$public_bcf" "$component_bcf" "$unfiltered_bcf" \
    "$filtered_bcf" "$final_bcf" "$meta_tsv" "$prsice_tsv" "$report"; do
    check_output "$output"
done

if [[ ! -s "$prepared" ]]; then
    log "preparing QC-passing 23andMe autosomal SNVs"
    if [[ "$assoc" == *.gz ]]; then gzip -t "$assoc"; fi
    if [[ "$annotation" == *.gz ]]; then gzip -t "$annotation"; fi
    prepared_tmp="$prepared.partial.$$"
    [[ ! -e "$prepared_tmp" ]] || die "partial output exists: $prepared_tmp"
    paste <(stream_file "$annotation") <(stream_file "$assoc") | \
        awk -v position_offset="$annotation_position_offset" -f "$prepare_awk" | \
        bgzip -@ "$threads" -c > "$prepared_tmp"
    mv "$prepared_tmp" "$prepared"
else
    log "reusing $prepared"
fi

if [[ ! -s "$public_bcf" ]]; then
    log "normalizing public no-23andMe component on GRCh37"
    bcftools +munge --no-version -Ou -C "$columns" -f "$ref37" \
        -s "$final_sample" "$public" | \
        bcftools view -Ou -f PASS,. -m2 -M2 -v snps | \
        bcftools norm -Ou -f "$ref37" -c e -d exact | \
        bcftools sort -Ob -o "$public_bcf" --write-index
else
    log "reusing $public_bcf"
    [[ -s "$public_bcf.csi" ]] || bcftools index -f "$public_bcf"
fi

if [[ ! -s "$component_bcf" ]]; then
    log "normalizing 23andMe component on GRCh37"
    bcftools +munge --no-version -Ou -C "$columns" -f "$ref37" \
        -s "$final_sample" "$prepared" | \
        bcftools view -Ou -f PASS,. -m2 -M2 -v snps | \
        bcftools norm -Ou -f "$ref37" -c e -d exact | \
        bcftools sort -Ob -o "$component_bcf" --write-index
else
    log "reusing $component_bcf"
    [[ -s "$component_bcf.csi" ]] || bcftools index -f "$component_bcf"
fi

[[ "$(bcftools query -l "$public_bcf")" == "$final_sample" ]] || \
    die "public BCF must contain exactly sample $final_sample"
[[ "$(bcftools query -l "$component_bcf")" == "$final_sample" ]] || \
    die "23andMe BCF must contain exactly sample $final_sample"

if [[ ! -s "$unfiltered_bcf" ]]; then
    log "running standard-error inverse-variance fixed-effect meta-analysis"
    bcftools +metal --no-version --het --esd --threads "$threads" \
        -Ob -o "$unfiltered_bcf" --write-index "$public_bcf" "$component_bcf"
else
    log "reusing $unfiltered_bcf"
    [[ -s "$unfiltered_bcf.csi" ]] || bcftools index -f "$unfiltered_bcf"
fi
[[ "$(bcftools query -l "$unfiltered_bcf")" == "$final_sample" ]] || \
    die "meta-analysis BCF must contain exactly sample $final_sample"

max_ne=$(bcftools query -f '[%NE\n]' "$unfiltered_bcf" | \
    awk '$1 != "." && $1+0 > max {max=$1+0} END{if(max>0) printf "%.10g",max; else exit 1}')
ne_threshold=$(awk -v max="$max_ne" -v fraction="$min_ne_fraction" \
    'BEGIN{printf "%.10g", max*fraction}')

if [[ ! -s "$filtered_bcf" ]]; then
    log "applying effective-N threshold: NE >= $ne_threshold ($min_ne_fraction of $max_ne)"
    filter_expression="FMT/NE >= $ne_threshold"
    if [[ -n "$exclude_list" ]]; then
        bcftools view -Ou -i "$filter_expression" "$unfiltered_bcf" | \
            bcftools view -T "^$exclude_list" -Ob -o "$filtered_bcf" --write-index
    else
        bcftools view -Ob -i "$filter_expression" -o "$filtered_bcf" \
            --write-index "$unfiltered_bcf"
    fi
else
    log "reusing $filtered_bcf"
    [[ -s "$filtered_bcf.csi" ]] || bcftools index -f "$filtered_bcf"
fi

if [[ ! -s "$final_bcf" ]]; then
    log "lifting the combined result once, from GRCh37 to GRCh38"
    bcftools +liftover --no-version -Ou "$filtered_bcf" -- \
        -s "$ref37" -f "$ref38" -c "$chain" | \
        bcftools annotate -Ou -h "$extra_headers" | \
        bcftools sort -Ob -o "$final_bcf" --write-index
else
    log "reusing $final_bcf"
    [[ -s "$final_bcf.csi" ]] || bcftools index -f "$final_bcf"
fi

if [[ ! -s "$meta_tsv" ]]; then
    log "exporting canonical meta-analysis table"
    {
        printf '#CHROM\tPOS\tRSID\tVARIANT\tREF\tALT\tEA\tNEA\tBETA\tSE\tP\tOR\tNS\tNC\tNE\n'
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NS\t%NC\t%NE]\n' \
            "$final_bcf" | awk 'BEGIN{FS=OFS="\t"} {
                p=exp(-2.302585092994046*$8); if(p==0) p=1e-300;
                variant=$1 ":" $2 ":" $4 ":" $5;
                print $1,$2,$3,variant,$4,$5,$5,$4,$6,$7,sprintf("%.12g",p),
                    sprintf("%.12g",exp($6)),$9,$10,$11
            }'
    } | bgzip -@ "$threads" -c > "$meta_tsv"
    tabix -f -s 1 -b 2 -e 2 "$meta_tsv"
else
    log "reusing $meta_tsv"
fi

if [[ ! -s "$prsice_tsv" ]]; then
    log "exporting PRSice base table"
    {
        printf 'CHR\tBP\tSNP\tA1\tA2\tBETA\tP\n'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%ES\t%LP]\n' "$final_bcf" | \
            awk 'BEGIN{FS=OFS="\t"} {
                p=exp(-2.302585092994046*$6); if(p==0) p=1e-300;
                print $1,$2,$1 ":" $2 ":" $3 ":" $4,$4,$3,$5,sprintf("%.12g",p)
            }'
    } | bgzip -@ "$threads" -c > "$prsice_tsv"
else
    log "reusing $prsice_tsv"
fi

records=$(bcftools index -n "$final_bcf")
gws=$(bcftools query -i 'FMT/LP>=7.30103' -f '%CHROM\n' "$final_bcf" | wc -l)
ne_ratio=$(awk -v observed="$max_ne" -v expected="$expected_max_ne" \
    'BEGIN{printf "%.5f",observed/expected}')
awk -v ratio="$ne_ratio" 'BEGIN{exit !(ratio < 0.95 || ratio > 1.05)}' && \
    log "WARNING: maximum NE differs from paper scale by more than 5% (ratio $ne_ratio)"

if ((skip_sha256 == 0)); then
    log "hashing source inputs"
    public_sha=$(sha256sum "$public" | awk '{print $1}')
    assoc_sha=$(sha256sum "$assoc" | awk '{print $1}')
    annotation_sha=$(sha256sum "$annotation" | awk '{print $1}')
else
    public_sha=not_calculated
    assoc_sha=not_calculated
    annotation_sha=not_calculated
fi

cat > "$report" <<EOF
# GWAS integration report: $trait

Created: $(date --iso-8601=seconds)

## Inputs

| Component | File | SHA-256 |
|---|---|---|
| Public EUR, no 23andMe | \`$public\` | \`$public_sha\` |
| 23andMe EUR association | \`$assoc\` | \`$assoc_sha\` |
| 23andMe annotation | \`$annotation\` | \`$annotation_sha\` |

The run asserted that the public and 23andMe components have no overlapping
participants. No genomic-control multiplier was applied to the delivered
23andMe beta or standard error.

## Method

- Effect model: standard-error inverse-variance fixed effect
- Variant scope: QC-passing, biallelic autosomal SNVs
- Effect scale: log odds
- Annotation-position offset applied before normalization: $annotation_position_offset
- GRCh37 effective-N inclusion fraction: $min_ne_fraction
- Observed maximum NE: $max_ne
- Expected paper-scale maximum NE: $expected_max_ne
- Observed/expected maximum NE: $ne_ratio
- Applied minimum NE: $ne_threshold
- Optional post-meta exclusion list: ${exclude_list:-none}

## Outputs

| Output | Path |
|---|---|
| Final GRCh37 BCF | \`$filtered_bcf\` |
| Final GRCh38 BCF | \`$final_bcf\` |
| Canonical GRCh38 table | \`$meta_tsv\` |
| PRSice GRCh38 base table | \`$prsice_tsv\` |

The final GRCh38 BCF contains $records variants, including $gws variants with
P <= 5e-8. These are variant counts, not independent signals or loci.

## Known Limits

- The supplied 23andMe annotations encode indels symbolically as D/I. This run
  excludes them because exact REF/ALT sequence alleles cannot be inferred from
  D/I alone. The public BD and MDD BCFs contain only 2 and 1 indels,
  respectively, so this does not remove a substantial shared variant class.
- $annotation_release_note
- The BD paper's post-meta DENTIST QC is not reproduced unless its rejected
  GRCh37 variants are supplied with \`--exclude-list\`.
- Heterogeneity across the two aggregate components is calculated, but the
  original cohort-level heterogeneity statistics cannot be reconstructed.
- \`FORMAT/SI\` is defined but missing. Public HRC INFO and 23andMe rsq are
  component-specific and are not combined into a fabricated quality score.
EOF

log "integration complete: $final_bcf"
