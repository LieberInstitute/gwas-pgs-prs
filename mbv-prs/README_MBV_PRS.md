# MBv PRS Results

## What Was Run

The approved plan is in `plan_26-05-10_20-37_master.md`.

The reproducible driver is:

```bash
./scripts/run_mbv_prs.sh
```

The actual run required two fixes now included in the script:

- Apply variant missingness before sample missingness. The input VCF is a union across imputation batches/arrays, so `--mind` on the unfiltered union removes every sample.
- Use a chr-prefixed BCF rebuilt from the normalized VCF for `bcftools +score`. PLINK2 BCF export changes contig names from `chr1` to `1`, breaking exact GWAS-VCF matching.

## Main Outputs

- `qc/target_qc_report.tsv`: sample, variant, relatedness, PCA, and overlap summary.
- `scores_bcftools/mbv_graphpred_scores.tsv`: GraphPred/LDGM scores from `bcftools +score`.
- `scores_bcftools/mbv_graphpred_scores_z.tsv`: z-scored GraphPred scores.
- `prsice/out/mbv_prsice_scores.tsv`: PRSice clumping-and-thresholding scores.
- `prsice/out/mbv_prsice_scores_z.tsv`: z-scored PRSice scores.
- `mbv_prs_merged_scores.tsv`: merged raw score table for all methods and diagnoses.
- `logs/commands_successful.tsv`: successful commands and non-trivial actions.

## QC Summary

Final target QC used autosomal biallelic SNPs, `R2 >= 0.8`, `MAF >= 0.01`, exact REF checks, DS range checks, variant missingness, sample missingness after variant filtering, control-only HWE, LD-pruned PCA, and KING relatedness.

Current summary:

```text
input samples: 119
final samples: 119
input variants: 14494877
final variants: 6639276
DS out-of-range sites: 0
KING related pairs >= 0.0884: 0
PCA samples: 119
GraphPred overlaps: BPD 4996043, MDD 4832784, SCZD 5128978
```

## Method Notes

`bcftools +pgs` produced LD-aware GraphPred loadings from the GWAS BCFs and EUR LDGM. `bcftools +score` only applies those loadings to MBv dosages.

PRSice used the converted GWAS BCFs as base data, not the GraphPred loading BCFs. It used the external EUR PLINK LD panel:

```text
/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos
```

PRSice cannot use the LDGM BCF files as its LD reference because PRSice clumping expects genotype-like PLINK/BGEN data, not sparse LDGM precision matrices.
