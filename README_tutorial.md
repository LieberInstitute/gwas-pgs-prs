# PRS Calculation Tutorial

This tutorial shows how to calculate PRS for a new GWAS using the BPD GWAS and MBv genotypes as example data.

The workflow standardizes the GWAS first, then uses the same GRCh38 GWAS BCF for both `bcftools` scoring and PRSice-2 scoring.

## 1. Review the example config

Example files are in:

```bash
cd /dbdata/cdb/gwas-pgs-prs
ls examples/mbv_bpd
```

The main config is:

```bash
examples/mbv_bpd/config.env
```

It defines tool paths, the work directory, target genotype file, reference FASTA, PRSice thresholds, and manifest paths.

For a new project, copy `examples/mbv_bpd/` to a new example directory and edit the paths.

The BPD/MDD/SCZD GWAS conversion performed previously used these local reference files:

```text
GRCh37 source FASTA: /dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta
GRCh38 liftover FASTA: /dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
UCSC liftover chain: /dbdata/cdb/ref/hg19ToHg38.over.chain.gz
```

The MBv target genotype QC step uses a separate GRCh38 FASTA:

```text
/dbdata/cdb/genotyping/ref/hg38c.fa
```

## 2. Prepare the GWAS column map

`bcftools +munge` needs to know how source GWAS columns map to canonical GWAS fields. This repository provides:

```bash
colheaders.tsv
```

For a new GWAS, inspect the header first:

```bash
zcat GWAS/BPD/bip2024_eur_no23andMe.gz | head -n 1
```

If the source uses different column names, update a project-specific copy of `colheaders.tsv`.

## 3. Convert the GWAS to GWAS-VCF/BCF

Edit:

```bash
examples/mbv_bpd/gwas_manifest.tsv
```

Set these fields correctly:

- `source`: source GWAS summary statistics.
- `sample_name`: synthetic sample name for the GWAS BCF.
- `source_build`: usually `GRCh37` or `GRCh38`.
- `source_fasta`: `/dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta` for the BPD/MDD/SCZD examples.
- `target_fasta`: `/dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` for GWAS liftover.
- `chain`: `/dbdata/cdb/ref/hg19ToHg38.over.chain.gz` for GRCh37/hg19 to GRCh38.
- `columns_file`: column map, usually `colheaders.tsv`.
- `output_bcf`: standardized GRCh38 GWAS BCF.

The BPD example row is:

```text
label	source	sample_name	source_build	source_fasta	target_fasta	chain	columns_file	columns_preset	output_bcf
BPD	/dbdata/cdb/gwas-pgs-prs/GWAS/BPD/bip2024_eur_no23andMe.gz	BIP_2024.EUR	GRCh37	/dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta	/dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna	/dbdata/cdb/ref/hg19ToHg38.over.chain.gz	/dbdata/cdb/gwas-pgs-prs/colheaders.tsv	.	/dbdata/cdb/gwas-pgs-prs/work/mbv_bpd/gwas/BPD.hg38.bcf
```

The prior psychiatric GWAS conversions used these synthetic sample names:

```text
BPD: BIP_2024.EUR
MDD: MDD_2025
SCZD: SCZ_2022.EUR
```

For SCZD, the source file's terminal `NEFF` header was rewritten to `NEFFDIV2` before `bcftools +munge`, matching the earlier workflow:

```bash
gzip -cd PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz | sed 's/NEFF$/NEFFDIV2/' | \
  /opt/sw/bin/bcftools +munge --no-version -Ou -C colheaders.tsv \
    -f /dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta \
    -s SCZ_2022.EUR
```

Run:

```bash
scripts/prepare_gwas_vcf.sh \
  --config examples/mbv_bpd/config.env \
  --manifest examples/mbv_bpd/gwas_manifest.tsv
```

This runs `bcftools +munge`, then `bcftools +liftover` when the source build is not already GRCh38.

Expected output:

```text
work/mbv_bpd/gwas/BPD.hg38.bcf
work/mbv_bpd/gwas/BPD.hg38.bcf.csi
```

Check the standardized fields:

```bash
/opt/sw/bin/bcftools view -h work/mbv_bpd/gwas/BPD.hg38.bcf | grep 'ID=ES\|ID=LP'
```

## 4. Compute LD-aware PGS loadings

Edit:

```bash
examples/mbv_bpd/pgs_manifest.tsv
```

For the BPD example, the LDGM list is:

```text
/dbdata/cdb/ref/lgdms/GRCh38/ldgm-vcfs.EUR.list
```

Run:

```bash
scripts/compute_pgs_loadings.sh \
  --config examples/mbv_bpd/config.env \
  --manifest examples/mbv_bpd/pgs_manifest.tsv
```

Expected output:

```text
work/mbv_bpd/pgs/BPD.pgs.bcf
work/mbv_bpd/pgs/BPD.pgs.bcf.csi
```

This file contains the variant weights used by `bcftools +score`.

## 5. Prepare target genotypes

The MBv example target VCF is:

```text
mbv-prs/genotypes/merged_R.8_MAF.01.RSann.vcf.gz
```

Run:

```bash
scripts/prepare_target_qc.sh --config examples/mbv_bpd/config.env
```

Expected key outputs:

```text
work/mbv_bpd/target/qc/mbv.qc.chr.bcf
work/mbv_bpd/target/prsice/mbv_qc.bgen
work/mbv_bpd/target/prsice/mbv_qc.sample
```

The `.qc.chr.bcf` file preserves chr-prefixed contigs for exact matching to GWAS-VCF/PGS BCF files.

## 6. Score samples with bcftools

Run:

```bash
scripts/score_bcftools.sh \
  --config examples/mbv_bpd/config.env \
  --manifest examples/mbv_bpd/score_manifest.tsv
```

Expected output:

```text
work/mbv_bpd/scores_bcftools/BPD.graphpred.tsv
```

This is a sample-level score table. It is not a SNP-level model table.

## 7. Build the PRSice base file

Run:

```bash
scripts/build_prsice_base_from_gwas_vcf.sh \
  --config examples/mbv_bpd/config.env \
  --manifest examples/mbv_bpd/prsice_base_manifest.tsv
```

Expected output:

```text
work/mbv_bpd/prsice/base/BPD.base.tsv.gz
```

The base file columns are:

```text
CHR BP SNP A1 A2 BETA P
```

The values are extracted from the standardized GWAS BCF:

- `A1=ALT`
- `A2=REF`
- `BETA=ES`
- `P=10^-LP`

## 8. Run PRSice-2

Run:

```bash
scripts/run_prsice_scores.sh \
  --config examples/mbv_bpd/config.env \
  --manifest examples/mbv_bpd/prsice_manifest.tsv
```

Expected outputs:

```text
work/mbv_bpd/prsice/out/BPD.all_score
work/mbv_bpd/prsice/out/BPD.prsice
work/mbv_bpd/prsice/out/BPD.snp
```

The wrapper uses fixed thresholds from `PRSICE_THRESHOLDS` in `config.env` and enables `--print-snp`.

## 9. Export PRSice model tables

Run:

```bash
Rscript scripts/export_prsice_model_tables.R \
  examples/mbv_bpd/prsice_model_manifest.tsv
```

Expected outputs:

```text
work/mbv_bpd/prsice_model_export/BPD.postclump_model.tsv.gz
work/mbv_bpd/prsice_model_export/BPD.threshold_model.tsv.gz
work/mbv_bpd/prsice_model_export/model_manifest.tsv
work/mbv_bpd/prsice_model_export/threshold_count_check.tsv
```

These are SNP-level model tables. They contain PRSice-retained SNPs and weights. They do not contain donor/sample PRS values.

## 10. Summarize and compare scores

Edit:

```bash
examples/mbv_bpd/score_long_manifest.tsv
```

Add one row per score column to compare. Then run:

```bash
Rscript scripts/summarize_prs_scores.R \
  examples/mbv_bpd/score_long_manifest.tsv \
  work/mbv_bpd/summary

Rscript scripts/compare_prs_methods.R \
  examples/mbv_bpd/score_long_manifest.tsv \
  work/mbv_bpd/comparison
```

Expected outputs include:

```text
work/mbv_bpd/summary/scores_long.tsv
work/mbv_bpd/summary/scores_raw.tsv
work/mbv_bpd/summary/scores_z.tsv
work/mbv_bpd/comparison/all_pairwise_correlations.tsv
work/mbv_bpd/comparison/combined_scores_raw.tsv
work/mbv_bpd/comparison/combined_scores_z.tsv
```

## 11. Run by stage

The stage wrapper can run one part at a time:

```bash
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage gwas
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage target
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage pgs
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage bcftools-score
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage prsice-base
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage prsice
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage prsice-models
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage summarize
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage compare
```

Use:

```bash
scripts/run_prs_pipeline.sh --config examples/mbv_bpd/config.env --stage all
```

only after all manifests have been checked.

## Notes

- GWAS conversion is the most dataset-specific step.
- Once a GWAS is in standardized GRCh38 BCF format, downstream scoring should be reproducible and mostly manifest-driven.
- Sample-level PRS outputs are separate from SNP-level model tables.
- For PRSice, the SNP-level model is the retained clumped SNP set plus the base-file weights used at each threshold.
- For the bcftools GraphPred path, the SNP-level model is the `bcftools +pgs` output BCF.
