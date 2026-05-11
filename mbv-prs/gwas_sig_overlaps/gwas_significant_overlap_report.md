# GWAS Significant Variant Overlap Report

Threshold: `P <= 5e-08`, equivalent to `LP >= 7.30103`.

Primary overlap unit: exact `CHROM:POS:REF:ALT` variant ID.
Secondary overlap unit: `CHROM:POS` position ID.

## Exact Variant Counts

        region count
      BPD_only  1246
      MDD_only 11809
     SCZD_only 17040
  BPD_MDD_only    94
 BPD_SCZD_only   872
 MDD_SCZD_only  1740
  BPD_MDD_SCZD   805
     BPD_total  3017
     MDD_total 14448
    SCZD_total 20457
   union_total 33606

## Position Counts

        region count
      BPD_only  1246
      MDD_only 11809
     SCZD_only 17037
  BPD_MDD_only    94
 BPD_SCZD_only   872
 MDD_SCZD_only  1740
  BPD_MDD_SCZD   805
     BPD_total  3017
     MDD_total 14448
    SCZD_total 20454
   union_total 33603

## Pairwise Overlaps

                        universe     pair intersection union    jaccard
 exact_variant_CHROM_POS_REF_ALT  BPD_MDD          899 16566 0.05426778
 exact_variant_CHROM_POS_REF_ALT BPD_SCZD         1677 21797 0.07693719
 exact_variant_CHROM_POS_REF_ALT MDD_SCZD         2545 32360 0.07864648
              position_CHROM_POS  BPD_MDD          899 16566 0.05426778
              position_CHROM_POS BPD_SCZD         1677 21794 0.07694778
              position_CHROM_POS MDD_SCZD         2545 32357 0.07865377

## Interpretation Notes

The BPD significant-variant count is much lower than MDD and SCZD in these local inputs:

- BPD: 3,017 significant variants from 6,938,758 finite-LP rows; max LP 16.783.
- MDD: 14,448 significant variants from 7,362,674 finite-LP rows; max LP 26.6519.
- SCZD: 20,457 significant variants from 7,658,367 finite-LP rows; max LP 39.184.

Possible explanations:

- The BPD input is `bip2024_eur_no23andMe.hg38.bcf`, so it excludes 23andMe and is EUR-only. That can materially reduce discovery power versus a fuller meta-analysis.
- These are raw significant variant rows, not independent loci. Large LD blocks can produce many significant neighboring variants for MDD or SCZD, inflating row counts relative to locus counts.
- Disorder GWAS files differ in effective sample size, case/control balance, phenotype definition, imputation, and study QC.
- Liftover/normalization and GWAS-VCF conversion can remove or mark records differently across files. The similar total row counts make a gross parsing failure unlikely, but the exact release and p-value field should still be checked against source documentation.
- If comparing to published locus counts, clump or fine-map first; this table counts every variant with `LP >= 7.30103`.

## Files

- `gwas_significant_variants_long.csv`
- `exact_variant_membership.csv`
- `position_membership.csv`
- `gwas_significant_overlap_counts.csv`
- `gwas_significant_pairwise_overlaps.csv`
- `gwas_significant_exact_variant_venn.svg`
- `gwas_significant_position_venn.svg`
