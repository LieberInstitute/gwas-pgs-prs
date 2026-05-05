# Tutorial: LDGM Setup, GWAS GRCh38 Conversion, and PGS Computation

This document replaces the audit-style work log `work_26-04-28_19-48_1817973.md` with a step-by-step tutorial. It preserves the exact commands used for the 2024 BIP, 2025 MDD, and 2022 SCZ EUR analyses, and it records file provenance more explicitly.

## 1. Software Context

This workflow assumes the locally built BCFtools and freescore plugins are installed here:

```bash
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
```

The installed BCFtools build used for the run was:

```text
bcftools 1.23.1
Using htslib 1.23.1
```

The relevant plugins are:

```text
munge
liftover
pgs
```

## 2. File Provenance

### LDGM matrices

The LDGM data archive used in this workflow was:

```text
/scratch/gpertea/osrc/bcftools/ldgms/ldgms.GRCh38.zip
```

Its upstream origin was the Broad SCORE software page:

```text
https://software.broadinstitute.org/software/score/
https://software.broadinstitute.org/software/score/ldgms.GRCh38.zip
```

The Broad page describes this file as:

```text
ldgms.GRCh38.zip - LDGM-VCF precision matrix files, approximately 3 GB, updated on 2022-11-29
```

The archive was already local during this implementation, so it was not network-downloaded during this run. To reproduce acquisition from scratch, use:

```bash
mkdir -p /scratch/gpertea/osrc/bcftools/ldgms
curl -L --fail --show-error \
  -o /scratch/gpertea/osrc/bcftools/ldgms/ldgms.GRCh38.zip \
  https://software.broadinstitute.org/software/score/ldgms.GRCh38.zip
```

### GWAS summary statistics

BIP 2024 EUR no23andMe was downloaded during this run:

```text
Study: Bipolar Disorder 2024
Publication DOI: 10.1038/s41586-024-08468-9
Dataset page from freescore examples: http://figshare.com/articles/dataset/bip2024/27216117
File URL used: https://ndownloader.figshare.com/files/49760772
Local file: ~/work/ref/GWAS/BPD/bip2024_eur_no23andMe.gz
```

MDD 2025 EUR no23andMe was already local:

```text
Study: Major Depressive Disorder 2025
Publication DOI in file header: 10.1016/j.cell.2024.12.002
Preprint DOI in file header: 10.1101/2024.04.29.24306535
Header sumstats URL: https://www.med.unc.edu/pgc/download-results/
Freescore examples file URL: http://figshare.com/ndownloader/files/51487019
Local file: ~/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz
```

SCZ 2022 European autosome was already local:

```text
Study: Schizophrenia 2022
Publication DOI from freescore examples: 10.1038/s41586-022-04434-5
Dataset page from freescore examples: http://figshare.com/articles/dataset/scz2022/19426775
Freescore examples file URL: http://figshare.com/ndownloader/files/34517828
Local file: ~/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz
```

### Reference files

These reference files were already local and were not downloaded during this run:

```text
GRCh37 source FASTA:
~/work/ref/GRCh37/human_g1k_v37.fasta

GRCh38 destination FASTA:
~/work/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

UCSC liftOver chain:
~/work/ref/hg19ToHg38.over.chain.gz
```

The chain file is the standard UCSC hg19-to-hg38 chain also referenced by `bcftools +liftover` help:

```text
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

The exact original download commands for the two FASTA files were not part of this work session or prior work log. Their filenames identify the intended assemblies/reference sets used by the freescore examples workflow.

## 3. Install LDGM Data

Unpack the Broad SCORE LDGM archive into the requested target location:

```bash
mkdir -p "$HOME/work/ref/lgdms/GRCh38" && unzip -o /scratch/gpertea/osrc/bcftools/ldgms/ldgms.GRCh38.zip -d "$HOME/work/ref/lgdms/GRCh38" && printf '%s\n' "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf" > "$HOME/work/ref/lgdms/GRCh38/ldgm-vcfs.EUR.list" && printf '%s\n' "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.AFR.bcf" "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.AMR.bcf" "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EAS.bcf" "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf" "$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.SAS.bcf" > "$HOME/work/ref/lgdms/GRCh38/ldgm-vcfs.ALL.list"
```

Expected files:

```text
~/work/ref/lgdms/GRCh38/1kg_ldgm.AFR.bcf
~/work/ref/lgdms/GRCh38/1kg_ldgm.AMR.bcf
~/work/ref/lgdms/GRCh38/1kg_ldgm.EAS.bcf
~/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf
~/work/ref/lgdms/GRCh38/1kg_ldgm.SAS.bcf
```

Each has a matching `.bcf.csi` index.

## 4. Download BIP GWAS Data

The BIP file was downloaded with:

```bash
set -euo pipefail
mkdir -p "$HOME/work/ref/GWAS/BPD"
if [ ! -s "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.gz" ]; then
  curl -L --fail --show-error -o "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.gz" https://ndownloader.figshare.com/files/49760772
fi
ls -lh "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.gz"
```

The downloaded file size was:

```text
332M /home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.gz
```

## 5. Create Column Header Mapping

Create a shared mapping file for `bcftools +munge`:

```bash
cat > "$HOME/work/ref/GWAS/colheaders.tsv" <<'EOF'
SNP	SNP
ID	SNP
CHR	CHR
CHROM	CHR
#CHROM	CHR
BP	BP
POS	BP
A1	A1
EA	A1
A2	A2
NEA	A2
OR	OR
BETA	BETA
SE	SE
P	P
PVAL	P
INFO	INFO
IMPINFO	INFO
NCAS	N_CAS
NCON	N_CON
Nca	N_CAS
Nco	N_CON
NEFF	NEFF
NEFFDIV2	NEFFDIV2
Neff_half	NEFFDIV2
HETI	HET_I2
HetISqt	HET_I2
HETPVAL	HET_P
HetPVa	HET_P
Direction	DIRE
HRC_FRQ_A1	FRQ
EOF
```

This mapping is shared across the three psychiatric GWAS inputs. It is not disorder-specific. `FCAS` and `FCON` are intentionally not mapped because they are case/control-specific allele frequencies, not one overall effect allele frequency.

## 6. Convert GWAS Files to GRCh38 GWAS-VCF BCF

### BIP 2024 EUR

Input columns:

```text
SNP CHR BP A1 A2 INFO OR SE P ngt Direction HetISqt HetDf HetPVa Nca Nco Neff_half HRC_FRQ_A1
```

Command:

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
COL="$HOME/work/ref/GWAS/colheaders.tsv"
SRC="$HOME/work/ref/GRCh37/human_g1k_v37.fasta"
DST="$HOME/work/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN="$HOME/work/ref/hg19ToHg38.over.chain.gz"
cd "$HOME/work/ref/GWAS/BPD"
$BCF +munge --no-version -Ou -C "$COL" -f "$SRC" -s BIP_2024.EUR bip2024_eur_no23andMe.gz | \
  $BCF +liftover --no-version -Ou -- -s "$SRC" -f "$DST" -c "$CHAIN" | \
  $BCF sort -o bip2024_eur_no23andMe.hg38.bcf -Ob --write-index
```

Liftover summary:

```text
Lines total/swapped/reference added/rejected: 6939126/14873/6/362
```

### MDD 2025 EUR

Input columns:

```text
#CHROM POS ID EA NEA BETA SE PVAL FCAS FCON IMPINFO NEFF NCAS NCON HETI HETDF HETPVAL
```

Command:

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
COL="$HOME/work/ref/GWAS/colheaders.tsv"
SRC="$HOME/work/ref/GRCh37/human_g1k_v37.fasta"
DST="$HOME/work/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN="$HOME/work/ref/hg19ToHg38.over.chain.gz"
cd "$HOME/work/ref/GWAS/MDD"
$BCF +munge --no-version -Ou -C "$COL" -f "$SRC" -s MDD_2025 pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz | \
  $BCF +liftover --no-version -Ou -- -s "$SRC" -f "$DST" -c "$CHAIN" | \
  $BCF sort -o pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf -Ob --write-index
```

Liftover summary:

```text
Lines total/swapped/reference added/rejected: 7363302/15504/4/624
```

### SCZ 2022 European autosome

Input columns:

```text
CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL NCAS NCON NEFF
```

The terminal `NEFF` column was rewritten to `NEFFDIV2`, following the freescore examples page.

Command:

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
COL="$HOME/work/ref/GWAS/colheaders.tsv"
SRC="$HOME/work/ref/GRCh37/human_g1k_v37.fasta"
DST="$HOME/work/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN="$HOME/work/ref/hg19ToHg38.over.chain.gz"
cd "$HOME/work/ref/GWAS/SCZD"
gzip -cd PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz | sed 's/NEFF$/NEFFDIV2/' | \
  $BCF +munge --no-version -Ou -C "$COL" -f "$SRC" -s SCZ_2022.EUR | \
  $BCF +liftover --no-version -Ou -- -s "$SRC" -f "$DST" -c "$CHAIN" | \
  $BCF sort -o PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf -Ob --write-index
```

Liftover summary:

```text
Lines total/swapped/reference added/rejected: 7659767/16099/120/1280
```

## 7. Compute PGS with EUR LDGM

All three GWAS datasets are EUR-only, so the EUR LDGM was used:

```bash
LDGM="$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf"
```

### BIP 2024 EUR PGS

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
LDGM="$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf"
cd "$HOME/work/ref/GWAS/BPD"
$BCF +pgs --no-version --beta-cov 5e-8 --max-alpha-hat2 0.001 --exclude 'FILTER="IFFY"' bip2024_eur_no23andMe.hg38.bcf "$LDGM" --output bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf --output-type b --log bip2024_eur_no23andMe.hg38.pgs.b5e-8.log --write-index
```

### MDD 2025 EUR PGS

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
LDGM="$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf"
cd "$HOME/work/ref/GWAS/MDD"
$BCF +pgs --no-version --beta-cov 2e-8 --max-alpha-hat2 0.0005 --exclude 'FILTER="IFFY"' pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf "$LDGM" --output pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf --output-type b --log pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.log --write-index
```

### SCZ 2022 EUR PGS

```bash
set -euo pipefail
export BCFTOOLS_PLUGINS=/opt/sw/bcf-plugins
BCF=/opt/sw/bin/bcftools
LDGM="$HOME/work/ref/lgdms/GRCh38/1kg_ldgm.EUR.bcf"
cd "$HOME/work/ref/GWAS/SCZD"
$BCF +pgs --no-version --beta-cov 2e-7 --max-alpha-hat2 0.002 --exclude 'FILTER="IFFY"' PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf "$LDGM" --output PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf --output-type b --log PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.log --write-index
```

## 8. QC Commands

Check LDGM files and EUR record count:

```bash
set -euo pipefail
cd "$HOME/work/ref/lgdms/GRCh38"
ls -lh 1kg_ldgm.{AFR,AMR,EAS,EUR,SAS}.bcf 1kg_ldgm.{AFR,AMR,EAS,EUR,SAS}.bcf.csi ldgm-vcfs.*.list
/opt/sw/bin/bcftools index -n 1kg_ldgm.EUR.bcf
```

Check converted GWAS BCFs:

```bash
BCF=/opt/sw/bin/bcftools
for f in "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf" "$HOME/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf" "$HOME/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf"; do
  echo "== $f =="
  test -s "$f" && test -s "$f.csi" && echo indexed=yes
  printf 'sample: '; $BCF query -l "$f"
  printf 'records: '; $BCF index -n "$f" 2>/dev/null
  printf 'first: '; $BCF view -H "$f" 2>/dev/null | sed -n '1p'
  printf 'REF_MISMATCH count: '; $BCF view -H -f REF_MISMATCH "$f" 2>/dev/null | wc -l
done
```

Check PGS BCFs:

```bash
BCF=/opt/sw/bin/bcftools
for f in "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf" "$HOME/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf" "$HOME/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf"; do
  echo "== $f =="
  test -s "$f" && test -s "$f.csi" && echo indexed=yes
  printf 'sample: '; $BCF query -l "$f"
  printf 'records: '; $BCF index -n "$f" 2>/dev/null
  printf 'first: '; $BCF view -H "$f" 2>/dev/null | sed -n '1p'
done
```

Extract PGS summaries from logs:

```bash
for f in "$HOME/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.log" "$HOME/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.log" "$HOME/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.log"; do
  echo "== $f =="
  perl -pe 's/\e\[[0-9;]*[A-Za-z]//g; s/\r/\n/g' "$f" | rg '^(=== PARAMETERS ===|alpha:|betaCov:|heritability|max alpha|random seed|sigmasq|=== SUMMARY ===|.*non_missing=|Tr\(|Advised options:)'
done
```

## 9. QC Results

LDGM EUR record count:

```text
1kg_ldgm.EUR.bcf records=7944261
```

Converted GWAS BCF results:

```text
BIP_2024.EUR records=6938764 REF_MISMATCH=0
MDD_2025 records=7362678 REF_MISMATCH=0
SCZ_2022.EUR records=7658487 REF_MISMATCH=0
```

PGS output record counts:

```text
BIP_2024.EUR_pgs_a0.5_b5e-08 records=5660257
MDD_2025_pgs_a0.5_b2e-08 records=5737561
SCZ_2022.EUR_pgs_a0.5_b2e-07 records=6076466
```

PGS log summaries:

```text
BIP_2024.EUR 1kg_ldgm.EUR.bcf non_missing=5660257 missing=1223752 lambda_GC=1.3901 sigmasqInf=4.19e-08
Tr(S_inf)=0.137235 Tr(S_non_inf)=0.163591 selected_effects=1169.12
Advised options: --alpha-param -0.5000 --beta-cov 8.641e-08 --max-alpha-hat2 0.0005094

MDD_2025 1kg_ldgm.EUR.bcf non_missing=5737561 missing=1146448 lambda_GC=1.6815 sigmasqInf=1.002e-08
Tr(S_inf)=0.0526465 Tr(S_non_inf)=0.0886379 selected_effects=1451
Advised options: --alpha-param -0.5000 --beta-cov 2.183e-08 --max-alpha-hat2 0.0001119

SCZ_2022.EUR 1kg_ldgm.EUR.bcf non_missing=6076466 missing=807543 lambda_GC=1.5708 sigmasqInf=7.684e-08
Tr(S_inf)=0.56841 Tr(S_non_inf)=0.119985 selected_effects=495.75
Advised options: --alpha-param -0.5000 --beta-cov 1.643e-07 --max-alpha-hat2 0.001492
```

## 10. Main Output Files

Converted GRCh38 GWAS BCFs:

```text
~/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf
~/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf
~/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf
```

PGS BCFs:

```text
~/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf
~/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf
~/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf
```

Each BCF has a matching `.csi` index.

## 11. Notes

- The old work log was accurate about the local LDGM archive used, but incomplete about upstream LDGM provenance.
- This tutorial records the LDGM origin as the Broad SCORE software page and its `ldgms.GRCh38.zip` download.
- The freescore examples page was used to select the GWAS download URLs and PGS parameters for BIP, MDD, and SCZ.
- The MDD and SCZ files were already local in `~/work/ref/GWAS`; their original local download commands were not part of this session.
