#!/usr/bin/env bash
set -euo pipefail

## local tools
BCF=/opt/sw/bin/bcftools
PLINK2=plink2
PLINK=plink
PRSICE=/opt/sw/bin/PRSice_linux

## local inputs
GT=genotypes/merged_R.8_MAF.01.RSann.vcf.gz
DEMOG=MBv_demographics_n119.tab
REF=/dbdata/cdb/genotyping/ref/hg38c.fa
LDREF=/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos

## summary and score inputs
BPD_GWAS=/home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf
MDD_GWAS=/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf
SCZD_GWAS=/home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf
BPD_PGS=/home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf
MDD_PGS=/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf
SCZD_PGS=/home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf

THREADS=${THREADS:-8}
AUTO_REGIONS=$(printf "chr%s," {1..22} | sed 's/,$//')

mkdir -p qc scores_bcftools prsice/base prsice/target prsice/ref prsice/out logs scripts
: > logs/commands_successful.tsv

run_cmd() {
  local name=$1
  shift
  printf '[%s] %s\n' "$(date '+%F %T')" "$name"
  "$@"
  printf '%s\t%s\n' "$(date '+%F %T')" "$name" >> logs/commands_successful.tsv
}

run_bash() {
  local name=$1
  local cmd=$2
  printf '[%s] %s\n' "$(date '+%F %T')" "$name"
  bash -o pipefail -c "$cmd"
  printf '%s\t%s\n' "$(date '+%F %T')" "$name" >> logs/commands_successful.tsv
}

## record versions and input identity
{
  echo "bcftools:"
  "$BCF" --version | head -n 2
  echo "plink2:"
  "$PLINK2" --version
  echo "plink:"
  "$PLINK" --version
  echo "PRSice:"
  "$PRSICE" --help 2>&1 | head -n 1
} > logs/tool_versions.txt

run_bash "write input checksums" \
  "sha256sum '$GT' '$GT.csi' '$DEMOG' '$BPD_GWAS' '$MDD_GWAS' '$SCZD_GWAS' '$BPD_PGS' '$MDD_PGS' '$SCZD_PGS' > logs/input_checksums.sha256"

## confirm sample IDs are identical between demographics and VCF
run_bash "write sample lists" \
  "awk -F '\t' 'NR>1{print \$1}' '$DEMOG' | sort > qc/demog.samples && '$BCF' query -l '$GT' | sort > qc/vcf.samples"
run_bash "compare sample lists" \
  "comm -3 qc/demog.samples qc/vcf.samples > qc/sample_mismatches.txt && test ! -s qc/sample_mismatches.txt"

## basic cohort counts
run_bash "write demographics summary" \
  "awk -F '\t' 'NR>1{n++; dx[\$5]++; race[\$3]++; sex[\$4]++} END{print \"metric\tlevel\tcount\"; print \"n\tall\t\" n; for(k in dx) print \"PrimaryDx\t\" k \"\t\" dx[k]; for(k in race) print \"Race\t\" k \"\t\" race[k]; for(k in sex) print \"Sex\t\" k \"\t\" sex[k]}' '$DEMOG' > qc/demographics_summary.tsv"

run_bash "write input variant count" \
  "printf 'stage\tvariants\n' > qc/variant_counts.tsv && printf 'input\t%s\n' \"\$('$BCF' index -n '$GT')\" >> qc/variant_counts.tsv"

## keep autosomal, biallelic, well-imputed common SNPs and normalize IDs
run_bash "create normalized autosomal target BCF" \
  "'$BCF' view --threads '$THREADS' -Ou -r '$AUTO_REGIONS' -f PASS,. -i 'N_ALT=1 && TYPE=\"snp\" && INFO/R2>=0.8 && INFO/MAF>=0.01' '$GT' | '$BCF' norm --threads '$THREADS' -Ou -f '$REF' -c e -d exact | '$BCF' annotate --threads '$THREADS' -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' -o qc/mbv.autosomal.snps.r8.maf01.norm.bcf --write-index"

## check dosage range before PLINK import
run_bash "write DS out of range sites" \
  "'$BCF' view -H -i 'MAX(FMT/DS)>2 || MIN(FMT/DS)<0' qc/mbv.autosomal.snps.r8.maf01.norm.bcf | cut -f1-5 > qc/ds_out_of_range_sites.tsv"

if [ -s qc/ds_out_of_range_sites.tsv ]; then
  run_bash "exclude DS out of range sites" \
    "cut -f1,2,4,5 qc/ds_out_of_range_sites.tsv | awk 'BEGIN{OFS=\"\t\"}{print \$1\":\"\$2\":\"\$3\":\"\$4}' > qc/ds_out_of_range.ids && '$BCF' view -Ob -e 'ID=@qc/ds_out_of_range.ids' -o qc/mbv.autosomal.snps.r8.maf01.norm.dsok.bcf --write-index qc/mbv.autosomal.snps.r8.maf01.norm.bcf"
else
  ln -sf mbv.autosomal.snps.r8.maf01.norm.bcf qc/mbv.autosomal.snps.r8.maf01.norm.dsok.bcf
  ln -sf mbv.autosomal.snps.r8.maf01.norm.bcf.csi qc/mbv.autosomal.snps.r8.maf01.norm.dsok.bcf.csi
fi

## convert dosages for QC and downstream export
run_cmd "import target to PLINK2 PGEN" \
  "$PLINK2" --bcf qc/mbv.autosomal.snps.r8.maf01.norm.dsok.bcf dosage=DS --double-id --make-pgen --out qc/mbv.ds

run_cmd "deduplicate target variants" \
  "$PLINK2" --pfile qc/mbv.ds --rm-dup force-first list --make-pgen --out qc/mbv.dedup

run_cmd "apply variant missingness and MAF QC" \
  "$PLINK2" --pfile qc/mbv.dedup --geno 0.02 --maf 0.01 --make-pgen --out qc/mbv.varqc

run_cmd "apply sample missingness QC after variant QC" \
  "$PLINK2" --pfile qc/mbv.varqc --mind 0.02 --make-pgen --out qc/mbv.qc1

## controls-only HWE variant filter, applied back to all samples
run_bash "write control keep file" \
  "awk -F '\t' 'NR>1 && \$5==\"Control\"{print \$1, \$1}' '$DEMOG' > qc/controls.keep"

run_cmd "write controls HWE pass list" \
  "$PLINK2" --pfile qc/mbv.qc1 --keep qc/controls.keep --hwe 1e-6 0 keep-fewhet --write-snplist --out qc/controls_hwe_pass

run_cmd "apply controls HWE pass list" \
  "$PLINK2" --pfile qc/mbv.qc1 --extract qc/controls_hwe_pass.snplist --make-pgen --out qc/mbv.qc

run_bash "append final variant count" \
  "grep -v '^final_qc' qc/variant_counts.tsv > qc/variant_counts.tmp && mv qc/variant_counts.tmp qc/variant_counts.tsv && printf 'final_qc\t%s\n' \"\$(awk '!/^#/{n++} END{print n}' qc/mbv.qc.pvar)\" >> qc/variant_counts.tsv"

## final QC metrics
run_cmd "write final sample missingness" \
  "$PLINK2" --pfile qc/mbv.qc --missing sample-only --out qc/mbv.final

run_cmd "write final variant missingness" \
  "$PLINK2" --pfile qc/mbv.qc --missing variant-only --out qc/mbv.final

run_cmd "write final allele frequencies" \
  "$PLINK2" --pfile qc/mbv.qc --freq --out qc/mbv.final

run_cmd "LD prune final target" \
  "$PLINK2" --pfile qc/mbv.qc --indep-pairwise 200 50 0.1 --out qc/mbv.prune

run_cmd "compute target PCs" \
  "$PLINK2" --pfile qc/mbv.qc --extract qc/mbv.prune.prune.in --pca 10 --out qc/mbv.pca

run_cmd "compute KING relatedness table" \
  "$PLINK2" --pfile qc/mbv.qc --extract qc/mbv.prune.prune.in --make-king-table --king-table-filter 0.0884 --out qc/mbv.king

## export final target BCF and BGEN for scoring methods
run_cmd "export QC target BCF with DS" \
  "$PLINK2" --pfile qc/mbv.qc --export bcf vcf-dosage=DS --out qc/mbv.qc
run_cmd "index QC target BCF" \
  "$BCF" index --threads "$THREADS" -f qc/mbv.qc.bcf

## preserve chr-prefixed VCF representation for exact GWAS-VCF scoring
run_bash "write final QC variant IDs" \
  "awk '!/^#/{print \$3}' qc/mbv.qc.pvar > qc/mbv.qc.ids"
run_bash "export chr-prefixed QC target BCF" \
  "'$BCF' view --threads '$THREADS' -Ob -i 'ID=@qc/mbv.qc.ids' -o qc/mbv.qc.chr.bcf --write-index qc/mbv.autosomal.snps.r8.maf01.norm.dsok.bcf"

run_cmd "export QC target BGEN for PRSice" \
  "$PLINK2" --pfile qc/mbv.qc --export bgen-1.2 bits=8 ref-first --out prsice/target/mbv_qc
run_bash "write dummy PRSice phenotype" \
  "printf 'IID DUMMY\n' > prsice/target/dummy.pheno && awk 'NR>2{print \$1\"_\"\$2, 1}' prsice/target/mbv_qc.sample >> prsice/target/dummy.pheno"

## score GraphPred loadings with bcftools
run_cmd "bcftools score BPD" \
  "$BCF" +score --use DS --sample-header --counts -o scores_bcftools/BPD.graphpred.tsv qc/mbv.qc.chr.bcf "$BPD_PGS"
run_cmd "bcftools score MDD" \
  "$BCF" +score --use DS --sample-header --counts -o scores_bcftools/MDD.graphpred.tsv qc/mbv.qc.chr.bcf "$MDD_PGS"
run_cmd "bcftools score SCZD" \
  "$BCF" +score --use DS --sample-header --counts -o scores_bcftools/SCZD.graphpred.tsv qc/mbv.qc.chr.bcf "$SCZD_PGS"

run_bash "write scoring overlap counts" \
  "printf 'disorder\tvariants\n' > qc/scoring_overlap.tsv; printf 'BPD\t%s\n' \"\$('$BCF' isec -n=2 -w1 qc/mbv.qc.chr.bcf '$BPD_PGS' | wc -l)\" >> qc/scoring_overlap.tsv; printf 'MDD\t%s\n' \"\$('$BCF' isec -n=2 -w1 qc/mbv.qc.chr.bcf '$MDD_PGS' | wc -l)\" >> qc/scoring_overlap.tsv; printf 'SCZD\t%s\n' \"\$('$BCF' isec -n=2 -w1 qc/mbv.qc.chr.bcf '$SCZD_PGS' | wc -l)\" >> qc/scoring_overlap.tsv"

## build PRSice base files from GWAS-VCF summary statistics
make_base() {
  local label=$1
  local gwas=$2
  local out=prsice/base/${label}.base.tsv
  printf 'CHR\tBP\tSNP\tA1\tA2\tBETA\tP\n' > "$out"
  "$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%ES\t%LP]\n' \
    -i 'N_ALT=1 && TYPE="snp" && FILTER!="IFFY" && FILTER!="REF_MISMATCH"' "$gwas" | \
    awk 'BEGIN{OFS="\t"} $5!="." && $6!="." {chr=$1; sub(/^chr/,"",chr); snp=$1":"$2":"$3":"$4; if(!seen[snp]++){p=10^(-$6); print chr,$2,snp,$4,$3,$5,p}}' >> "$out"
  gzip -f "$out"
}

make_base BPD "$BPD_GWAS"
printf '%s\t%s\n' "$(date '+%F %T')" "make PRSice BPD base" >> logs/commands_successful.tsv
make_base MDD "$MDD_GWAS"
printf '%s\t%s\n' "$(date '+%F %T')" "make PRSice MDD base" >> logs/commands_successful.tsv
make_base SCZD "$SCZD_GWAS"
printf '%s\t%s\n' "$(date '+%F %T')" "make PRSice SCZD base" >> logs/commands_successful.tsv

## run PRSice fixed-threshold scores only
THRESHOLDS=5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1
for label in BPD MDD SCZD; do
  run_cmd "PRSice ${label}" \
    "$PRSICE" \
      --base "prsice/base/${label}.base.tsv.gz" \
      --target prsice/target/mbv_qc \
      --type bgen \
      --ld "$LDREF" \
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
      --pheno prsice/target/dummy.pheno \
      --pheno-col DUMMY \
      --no-regress \
      --all-score \
      --score sum \
      --fastscore \
      --bar-levels "$THRESHOLDS" \
      --clump-kb 250 \
      --clump-r2 0.1 \
      --clump-p 1 \
      --thread "$THREADS" \
      --out "prsice/out/${label}"
done

## summarize outputs with base R to avoid package dependencies
Rscript scripts/summarize_mbv_prs.R
printf '%s\t%s\n' "$(date '+%F %T')" "summarize PRS outputs" >> logs/commands_successful.tsv

printf '[%s] done\n' "$(date '+%F %T')"
