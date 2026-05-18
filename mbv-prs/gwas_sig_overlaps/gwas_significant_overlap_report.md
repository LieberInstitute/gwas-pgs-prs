# GWAS Significant Variant Overlap Report

Threshold: `P <= 5e-08`, equivalent to `LP >= 7.30103`.

Primary overlap unit: exact canonical `CHROM:POS:REF:ALT` variant ID.
Secondary overlap unit: `CHROM:POS` position ID.
The exact and position diagrams are nearly identical because only three SCZD-only positions collapse when REF/ALT is ignored.

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

## Files

- `gwas_significant_variants_long.csv`
- `exact_variant_membership.csv`
- `position_membership.csv`
- `gwas_significant_overlap_counts.csv`
- `gwas_significant_pairwise_overlaps.csv`
- `gwas_significant_exact_variant_venn.svg`
- `gwas_significant_position_venn.svg`
