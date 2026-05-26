#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
## shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat <<'EOF'
usage: prepare_target_qc.sh --config CONFIG [--dry-run]

Required config variables:
TARGET_VCF REF_FASTA WORK_DIR TARGET_PREFIX

Optional config variables:
BCFTOOLS PLINK2 THREADS AUTO_REGIONS TARGET_FILTER DEMOGRAPHICS PHENO_ID_COL PHENO_GROUP_COL HWE_CONTROL_VALUE
EOF
}

CONFIG=
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG=$2; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[ -n "$CONFIG" ] || die "missing --config"
source_config "$CONFIG"

BCFTOOLS=${BCFTOOLS:-bcftools}
PLINK2=${PLINK2:-plink2}
THREADS=${THREADS:-4}
AUTO_REGIONS=${AUTO_REGIONS:-$(printf "chr%s," {1..22} | sed 's/,$//')}
TARGET_FILTER=${TARGET_FILTER:-'N_ALT=1 && TYPE="snp" && INFO/R2>=0.8 && INFO/MAF>=0.01'}
TARGET_PREFIX=${TARGET_PREFIX:?}
WORK_DIR=${WORK_DIR:?}
TARGET_VCF=${TARGET_VCF:?}
REF_FASTA=${REF_FASTA:?}

require_cmd "$BCFTOOLS"
require_cmd "$PLINK2"
require_file "$TARGET_VCF"
require_file "$REF_FASTA"

qc_dir="$WORK_DIR/target/qc"
prsice_target_dir="$WORK_DIR/target/prsice"
mkdir -p "$qc_dir" "$prsice_target_dir"

norm_bcf="$qc_dir/${TARGET_PREFIX}.autosomal.snps.norm.bcf"
dsok_bcf="$qc_dir/${TARGET_PREFIX}.autosomal.snps.norm.dsok.bcf"
qc_pfx="$qc_dir/${TARGET_PREFIX}.qc"

msg "normalize target variants"
run_bash "'$BCFTOOLS' view --threads '$THREADS' -Ou -r '$AUTO_REGIONS' -f PASS,. -i '$TARGET_FILTER' '$TARGET_VCF' | '$BCFTOOLS' norm --threads '$THREADS' -Ou -f '$REF_FASTA' -c e -d exact | '$BCFTOOLS' annotate --threads '$THREADS' -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' -o '$norm_bcf' --write-index"

msg "check target DS range"
run_bash "'$BCFTOOLS' view -H -i 'MAX(FMT/DS)>2 || MIN(FMT/DS)<0' '$norm_bcf' | cut -f1-5 > '$qc_dir/ds_out_of_range_sites.tsv'"

if [ "${DRY_RUN:-0}" = "1" ] || [ ! -s "$qc_dir/ds_out_of_range_sites.tsv" ]; then
  run_bash "ln -sf '$(basename "$norm_bcf")' '$dsok_bcf' && ln -sf '$(basename "$norm_bcf").csi' '$dsok_bcf.csi'"
else
  run_bash "awk 'BEGIN{OFS=\"\\t\"}{print \$1\":\"\$2\":\"\$4\":\"\$5}' '$qc_dir/ds_out_of_range_sites.tsv' > '$qc_dir/ds_out_of_range.ids' && '$BCFTOOLS' view -Ob -e 'ID=@$qc_dir/ds_out_of_range.ids' -o '$dsok_bcf' --write-index '$norm_bcf'"
fi

msg "import target to PLINK2"
run_cmd "$PLINK2" --bcf "$dsok_bcf" dosage=DS --double-id --make-pgen --out "$qc_dir/${TARGET_PREFIX}.ds"
run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.ds" --rm-dup force-first list --make-pgen --out "$qc_dir/${TARGET_PREFIX}.dedup"
run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.dedup" --geno 0.02 --maf 0.01 --make-pgen --out "$qc_dir/${TARGET_PREFIX}.varqc"
run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.varqc" --mind 0.02 --make-pgen --out "$qc_dir/${TARGET_PREFIX}.qc1"

if [ -n "${DEMOGRAPHICS:-}" ] && [ -n "${PHENO_ID_COL:-}" ] && [ -n "${PHENO_GROUP_COL:-}" ] && [ -n "${HWE_CONTROL_VALUE:-}" ]; then
  require_file "$DEMOGRAPHICS"
  msg "apply controls-only HWE filter"
  run_bash "awk -F '\\t' -v id='$PHENO_ID_COL' -v grp='$PHENO_GROUP_COL' -v val='$HWE_CONTROL_VALUE' 'NR==1{for(i=1;i<=NF;i++){if(\$i==id) idc=i; if(\$i==grp) grpc=i} next} \$grpc==val{print \$idc, \$idc}' '$DEMOGRAPHICS' > '$qc_dir/controls.keep'"
  run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.qc1" --keep "$qc_dir/controls.keep" --hwe 1e-6 0 keep-fewhet --write-snplist --out "$qc_dir/controls_hwe_pass"
  run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.qc1" --extract "$qc_dir/controls_hwe_pass.snplist" --make-pgen --out "$qc_pfx"
else
  run_cmd "$PLINK2" --pfile "$qc_dir/${TARGET_PREFIX}.qc1" --make-pgen --out "$qc_pfx"
fi

msg "write final target exports"
run_cmd "$PLINK2" --pfile "$qc_pfx" --export bcf vcf-dosage=DS --out "$qc_pfx"
run_cmd "$BCFTOOLS" index --threads "$THREADS" -f "$qc_pfx.bcf"
run_bash "awk '!/^#/{print \$3}' '$qc_pfx.pvar' > '$qc_dir/${TARGET_PREFIX}.qc.ids'"
run_bash "'$BCFTOOLS' view --threads '$THREADS' -Ob -i 'ID=@$qc_dir/${TARGET_PREFIX}.qc.ids' -o '$qc_dir/${TARGET_PREFIX}.qc.chr.bcf' --write-index '$dsok_bcf'"
run_cmd "$PLINK2" --pfile "$qc_pfx" --export bgen-1.2 bits=8 ref-first --out "$prsice_target_dir/${TARGET_PREFIX}_qc"
run_bash "printf 'IID DUMMY\\n' > '$prsice_target_dir/dummy.pheno' && awk 'NR>2{print \$1\"_\"\$2, 1}' '$prsice_target_dir/${TARGET_PREFIX}_qc.sample' >> '$prsice_target_dir/dummy.pheno'"
