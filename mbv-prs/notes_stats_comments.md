# Notes: PRS and GWAS Overlap Follow-up Comments

This file summarizes observations and discussion points added after the main workbook/worklog files were written.

## PRS Method Correlation Check

Correlation/verification was performed between the GraphPred/+score PRS and PRSice PRS outputs for all three disorders.

Output file:

- `qc/prs_method_correlations.tsv`

Best observed Pearson correlations between GraphPred/+score and PRSice threshold series:

| disorder | best PRSice threshold | Pearson r |
|---|---:|---:|
| BPD | 0.1 | 0.828 |
| MDD | 0.5 | 0.798 |
| SCZD | 0.05 | 0.825 |

Interpretation:

- The methods are directionally consistent but not identical.
- This is expected because GraphPred/+score used precomputed LD-aware PGS loadings, while PRSice used clumping plus p-value thresholding.
- The comparison should be treated as a sanity check, not as evidence that the two methods estimate the same PRS model.

## +score Versus PRSice Threshold Series

The bcftools `+score` run used here produced one PRS series per disorder from the supplied score/weight file.

PRSice produced multiple PRS series because it explicitly ran p-value thresholding over several `Pt` values.

Practical implication:

- `+score` does not automatically emulate PRSice-style multiple p-value thresholds from one score file.
- To produce multiple thresholded `+score` outputs, prepare multiple score files, one per threshold, or use a score table organized for threshold-specific scoring if supported by the chosen plugin/options.
- With the current precomputed GraphPred/LD-aware loadings, one score per disorder is the expected behavior.

## GWAS Significant Variant Overlap Report

GWAS significant-variant overlap outputs were generated under:

- `gwas_sig_overlaps/`

Main files:

- `gwas_sig_overlaps/gwas_significant_overlap_report.md`
- `gwas_sig_overlaps/gwas_significant_overlap_counts.csv`
- `gwas_sig_overlaps/gwas_significant_pairwise_overlaps.csv`
- `gwas_sig_overlaps/gwas_significant_variants_long.csv`
- `gwas_sig_overlaps/exact_variant_membership.csv`
- `gwas_sig_overlaps/position_membership.csv`
- `gwas_sig_overlaps/gwas_significant_exact_variant_venn.svg`
- `gwas_sig_overlaps/gwas_significant_position_venn.svg`

Threshold used:

- `P <= 5e-8`
- Equivalent GWAS-VCF field threshold: `LP >= 7.30103`

Primary overlap unit:

- Exact variant ID: `CHROM:POS:REF:ALT`

Secondary overlap unit:

- Position ID: `CHROM:POS`

Exact-variant counts:

| region | count |
|---|---:|
| BPD only | 1,246 |
| MDD only | 11,809 |
| SCZD only | 17,040 |
| BPD and MDD only | 94 |
| BPD and SCZD only | 872 |
| MDD and SCZD only | 1,740 |
| BPD and MDD and SCZD | 805 |
| BPD total | 3,017 |
| MDD total | 14,448 |
| SCZD total | 20,457 |
| union total | 33,606 |

Pairwise exact-variant overlaps:

| pair | intersection | union | Jaccard |
|---|---:|---:|---:|
| BPD-MDD | 899 | 16,566 | 0.0543 |
| BPD-SCZD | 1,677 | 21,797 | 0.0769 |
| MDD-SCZD | 2,545 | 32,360 | 0.0786 |

The Venn diagrams are schematic and are not area-proportional.

## BPD Low Significant-Variant Count

BPD had fewer genome-wide significant variant rows than MDD or SCZD:

| disorder | finite LP rows | significant variants | max LP |
|---|---:|---:|---:|
| BPD | 6,938,758 | 3,017 | 16.783 |
| MDD | 7,362,674 | 14,448 | 26.6519 |
| SCZD | 7,658,367 | 20,457 | 39.184 |

This does not look like an obvious parsing failure because the total finite-LP row counts are in the same range across disorders.

Likely or plausible explanations:

- The BPD input file is `bip2024_eur_no23andMe.hg38.bcf`, but MDD also uses a `no23andMe` file. Therefore, 23andMe exclusion alone cannot explain why BPD is much lower than MDD.
- BPD is an outlier relative to MDD and SCZD in this specific raw significant-variant count.
- The more likely explanation is lower BPD discovery power or different BPD source-study composition, effective sample size, phenotype definition, case/control balance, or QC, not simply 23andMe exclusion.
- These counts are raw significant variant rows, not independent loci.
- MDD and SCZD may contain larger LD blocks where many correlated variants cross the threshold, inflating raw row counts.
- The three GWAS differ in effective sample size, case/control balance, phenotype definition, imputation, source cohorts, and source QC.
- Published numbers may refer to independent loci, all-ancestry analyses, full meta-analyses, or releases that include 23andMe.
- Liftover, normalization, allele harmonization, GWAS-VCF conversion, and filtering can affect the count differently across source files.

## Next Audit Suggestions

Recommended follow-up checks before treating the BPD count as biological rather than file/release-specific:

1. Confirm that `bip2024_eur_no23andMe.hg38.bcf` is the intended BPD GWAS release for this analysis.
2. Compare counts against the original BPD GWAS paper or release notes, specifically checking effective sample size, included cohorts, ancestry subset, and whether the local file is a reduced-release or public-only subset.
3. Check whether the source documentation reports independent loci rather than raw significant variants.
4. Run LD clumping on significant variants for all three disorders and compare independent-locus counts.
5. Count significant variants after restricting to `FILTER=PASS` only, and compare with the current all-filter count.
6. Inspect records marked `IFFY` or `REF_MISMATCH` to see whether they disproportionately affect one disorder.
7. Verify that `LP` is the correct p-value-derived field for all three GWAS-VCF files.
8. Compare exact `CHROM:POS:REF:ALT` overlap with rsID overlap where rsIDs are reliable.
9. Confirm allele normalization against the target genotype reference build and target variant IDs.
10. Preserve both raw-variant and LD-clumped/locus-level summaries in production reports.

## Production Workflow Notes

For production PRS workflows, keep these checks explicit:

- Record the exact GWAS file name, checksum, genome build, ancestry subset, and cohort exclusions.
- Record whether 23andMe or other restricted cohorts are included.
- Report raw significant-variant counts separately from LD-clumped independent loci.
- Record the threshold definition, p-value field, and whether all filters or only `PASS` variants were counted.
- Keep exact-variant and position-level overlap tables, because allele normalization can affect exact overlap.
- Treat method-to-method PRS correlations as QC diagnostics, not as validation that two methods are equivalent.

## Diagnosis Association Follow-up

Outputs were added under:

- `prs_diagnosis_assoc/`

Main files:

- `prs_diagnosis_assoc/prs_diagnosis_association_report.md`
- `prs_diagnosis_assoc/prs_association_all_scores.csv`
- `prs_diagnosis_assoc/prs_best_thresholds_exploratory.csv`
- `prs_diagnosis_assoc/prs_group_summary.csv`
- `prs_diagnosis_assoc/prs_overlap_interpretation.csv`
- `prs_diagnosis_assoc/prs_auc_by_threshold.svg`
- `prs_diagnosis_assoc/prs_selected_distributions.svg`

Demographics:

- 119 samples total.
- 39 MDD, 40 Bipolar/BPD, 40 Controls.
- No SCZ diagnosis group is present, so SCZD PRS was evaluated only as a cross-disorder score.

Models:

- Unadjusted logistic model: `case_status ~ PRS_z`.
- Sensitivity model: `case_status ~ PRS_z + Sex + PC1 + PC2 + PC3 + PC4 + PC5`.
- PRS values were standardized over all 119 samples.

Key exploratory observations:

- MDD vs Control showed weak separation:
  - SCZD PRSice Pt=0.1: AUC 0.566, OR 1.335 per 1 SD, 95% CI 0.843-2.115, p=0.218, within-contrast FDR 0.959, Cohen d 0.279.
  - MDD PRSice Pt=5e-8: AUC 0.553, OR 1.105, 95% CI 0.747-1.634, p=0.617, FDR 0.959, Cohen d 0.111.
  - MDD GraphPred: AUC 0.547, OR 1.329, 95% CI 0.829-2.132, p=0.237, FDR 0.959, Cohen d 0.267.
- BPD vs Control showed the clearest separation:
  - MDD PRSice Pt=0.0001: AUC 0.692, OR 2.029, 95% CI 1.193-3.448, p=0.00897, FDR 0.0573, Cohen d 0.640.
  - MDD GraphPred: AUC 0.687, OR 2.101, 95% CI 1.246-3.543, p=0.00536, FDR 0.0573, Cohen d 0.688.
  - SCZD PRSice Pt=0.05: AUC 0.676, OR 2.067, 95% CI 1.232-3.468, p=0.00594, FDR 0.0573, Cohen d 0.678.
  - SCZD GraphPred: AUC 0.668, OR 1.932, 95% CI 1.174-3.181, p=0.00959, FDR 0.0573, Cohen d 0.624.
  - SCZD PRSice Pt=0.1: AUC 0.667, OR 1.944, 95% CI 1.174-3.217, p=0.00975, FDR 0.0573, Cohen d 0.629.
  - BPD PRSice Pt=5e-8: AUC 0.634, OR 1.952, 95% CI 1.069-3.564, p=0.0295, FDR 0.0679, Cohen d 0.521.
- BPD vs MDD showed exploratory separation:
  - BPD PRSice Pt=1e-6: AUC 0.679, OR 1.988, 95% CI 1.184-3.339, p=0.00936, FDR 0.0705, Cohen d 0.646.
  - MDD PRSice Pt=0.1: AUC 0.675, OR 1.974, 95% CI 1.192-3.267, p=0.00818, FDR 0.0705, Cohen d 0.645.
  - MDD PRSice Pt=0.5: AUC 0.674, OR 1.958, 95% CI 1.176-3.258, p=0.00972, FDR 0.0705, Cohen d 0.625.
  - MDD PRSice Pt=1: AUC 0.663, OR 1.894, 95% CI 1.141-3.142, p=0.0135, FDR 0.0705, Cohen d 0.592.
  - MDD PRSice Pt=0.0001: AUC 0.658, OR 1.907, 95% CI 1.131-3.217, p=0.0155, FDR 0.0705, Cohen d 0.582.
- MDD+BPD vs Control showed modest separation:
  - MDD GraphPred: AUC 0.618, OR 1.633, 95% CI 1.078-2.472, p=0.0206, FDR 0.303, Cohen d 0.468.
  - SCZD PRSice Pt=0.1: AUC 0.617, OR 1.616, 95% CI 1.068-2.443, p=0.0230, FDR 0.303, Cohen d 0.457.
  - MDD PRSice Pt=0.0001: AUC 0.614, OR 1.436, 95% CI 0.958-2.153, p=0.0796, FDR 0.421, Cohen d 0.346.
  - SCZD PRSice Pt=0.05: AUC 0.612, OR 1.570, 95% CI 1.044-2.362, p=0.0303, FDR 0.303, Cohen d 0.434.
  - SCZD GraphPred: AUC 0.598, OR 1.413, 95% CI 0.955-2.092, p=0.0841, FDR 0.421, Cohen d 0.341.
- Output coverage checks after implementation:
  - `prs_association_all_scores.csv`: 240 rows.
  - `prs_best_thresholds_exploratory.csv`: 24 rows.
  - `prs_group_summary.csv`: 120 rows.
  - Contrasts: MDD vs Control, BPD vs Control, MDD+BPD vs Control, BPD vs MDD.
  - Models: unadjusted and adjusted for Sex plus PC1-PC5.

Threshold-selection comments:

- The internally best PRSice thresholds should not be used as confirmatory downstream thresholds without external validation.
- Choosing the threshold that best separates cases and controls in these same 119 samples overfits the MBv sample.
- For downstream analyses, prefer one of these strategies:
  1. Use a pre-specified score or threshold from literature or an external validation cohort.
  2. Use the GraphPred/LD-aware score as the main score and PRSice thresholds as sensitivity analyses.
  3. If using PRSice threshold optimization, select the threshold in an independent validation set, then test in MBv or another held-out set.
  4. If no independent sample exists, report all thresholds and label the internally best threshold exploratory.
  5. Use PRSice permutation/empirical p-values for optimized-threshold association testing, while still noting that optimized R2/AUC remains overfit.

Relation to GWAS significant-variant overlap:

- MDD-BPD exact significant-variant overlap was 899 variants.
- Jaccard overlap was 0.0543.
- Because BPD had fewer significant variants, 899 shared variants represent 29.8% of BPD significant variants but only 6.2% of MDD significant variants.
- This overlap helps contextualize cross-disorder PRS behavior but does not directly predict PRS separation, because PRSice thresholds include many sub-significant variants and GraphPred uses LD-aware genome-wide weights.
