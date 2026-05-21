# PRS Method Comparison

## Purpose

This directory compares three PRS score sources for the same 119 MBv samples:

- `Shizhong_C+T`: colleague-provided clumping-and-thresholding scores from `PRS_Bipolar.csv`, `PRS_MDD.csv`, and `PRS_SCZ.csv`.
- `PRSice_C+T`: our PRSice-2 clumping-and-thresholding scores from `../prsice/out/mbv_prsice_scores.tsv`.
- `GraphPred_bcftools`: our LDGM/GraphPred score loadings applied with `bcftools +score`, from `../scores_bcftools/mbv_graphpred_scores.tsv`.

The goal is not to declare one score "correct." The goal is to quantify how similarly the methods rank the same samples, identify which p-value thresholds are directly comparable, and prepare tables that can be shown or reused in reports.

## Method Context

Shizhong and PRSice are both LD-aware C+T approaches. They create multiple scores per disorder by using LD clumping to select approximately independent variants and then applying p-value thresholds. They are comparable in concept, but not identical:

- Shizhong uses custom PLINK 1.9 clumping plus PLINK2 scoring.
- PRSice uses PRSice-2 with our validated base files and configured LD reference.
- Shizhong used an 18-threshold grid.
- PRSice used a 9-threshold grid.
- Shizhong and PRSice share only 8 exact thresholds.

GraphPred is different in how it uses LD. C+T methods use LD for variant selection; GraphPred/LDGM uses LD structure upstream to estimate adjusted genome-wide score weights. `bcftools +score` then applies those weights. GraphPred gives one score per disorder rather than a threshold grid, so it should be compared with C+T scores by correlation across all thresholds, not by threshold matching.

## What `compare_prs_methods.R` Did

The script `../compare_prs_methods.R` performed these steps:

1. Read inputs:
   - `../MBv_demographics_n119.tab`
   - `../scores_bcftools/mbv_graphpred_scores.tsv`
   - `../prsice/out/mbv_prsice_scores.tsv`
   - `../PRS_Bipolar.csv`
   - `../PRS_MDD.csv`
   - `../PRS_SCZ.csv`
2. Checked that every score source contains the same 119 sample IDs as the demographics file.
3. Reordered all score tables to the demographics sample order.
4. Standardized method labels:
   - `Shizhong_C+T`
   - `PRSice_C+T`
   - `GraphPred_bcftools`
5. Standardized disorder labels:
   - `BPD`
   - `MDD`
   - `SCZD`
6. Parsed p-value thresholds for Shizhong and PRSice.
7. Z-scored every score column over all 119 samples.
8. Built combined raw and z-scored score matrices.
9. Computed:
   - PRSice vs Shizhong correlations at exact shared thresholds only.
   - GraphPred vs every PRSice and Shizhong threshold.
   - all within-disorder pairwise score correlations.
   - diagnosis group summaries using z-scored scores.

The script uses base R only.

## Generated Files

- `score_inventory.tsv`
  - One row per method and disorder.
  - Documents input file, score columns, threshold count, and sample count.

- `threshold_availability.tsv`
  - One row per disorder and p-value threshold.
  - Flags whether the threshold exists in Shizhong, PRSice, or both.
  - Use this table to explain why only some C+T thresholds are compared directly.

- `combined_scores_raw.tsv`
  - One row per sample.
  - Demographics plus all raw Shizhong, PRSice, and GraphPred scores.

- `combined_scores_z.tsv`
  - Same shape as `combined_scores_raw.tsv`.
  - Every score column is standardized across all 119 samples.
  - Best table for comparing score scales visually or summarizing by diagnosis.

- `matched_threshold_correlations.tsv`
  - PRSice vs Shizhong only at exact shared thresholds.
  - Columns include Pearson, Spearman, and z-score RMSE.
  - This is the main table for direct C+T method comparison.

- `graphpred_threshold_correlations.tsv`
  - GraphPred vs every PRSice and Shizhong threshold.
  - Includes a `best_abs_pearson_for_method` flag for the threshold most correlated with GraphPred within each disorder/method.

- `all_pairwise_correlations.tsv`
  - Complete within-disorder correlation grid for audit.
  - Useful if a reviewer asks about a threshold or method pair not shown in the summary.

- `diagnosis_group_summary.tsv`
  - Z-score mean, SD, median, and IQR by method, disorder, threshold, and `PrimaryDx`.
  - Descriptive only; not a confirmatory association test.

## Threshold Matching

PRSice and Shizhong share these 8 thresholds:

```text
1e-06
1e-04
0.001
0.01
0.05
0.1
0.5
1
```

PRSice-only threshold:

```text
5e-08
```

Shizhong-only thresholds:

```text
1e-08
1e-07
1e-05
0.2
0.3
0.4
0.6
0.7
0.8
0.9
```

Conclusion: compare PRSice vs Shizhong directly only at the 8 shared thresholds. Do not force `5e-08` to match `1e-08` or `1e-07`.

## Main Results

All comparisons used `n = 119` samples.

### PRSice vs Shizhong

PRSice and Shizhong agree strongly at shared thresholds, consistent with both being C+T methods.

Best Pearson correlations by disorder:

```text
BPD   p <= 0.1    r = 0.9296
MDD   p <= 1      r = 0.9307
SCZD  p <= 1e-06  r = 0.9242
```

At `p <= 1e-06`, the first shared threshold:

```text
BPD   r = 0.9002, Spearman = 0.9119
MDD   r = 0.8957, Spearman = 0.8788
SCZD  r = 0.9242, Spearman = 0.9080
```

Interpretation: the two C+T workflows produce highly similar sample rankings, but they are not identical. Differences are expected because Shizhong and PRSice used different software, preprocessing, LD clumping windows, score construction details, and threshold grids.

### GraphPred vs C+T

GraphPred is one LD-model-based score per disorder, not a threshold grid. Its strongest C+T correlations were:

```text
BPD vs PRSice     p <= 0.1    r = 0.8289
BPD vs Shizhong   p <= 0.1    r = 0.8198
MDD vs PRSice     p <= 0.5    r = 0.7973
MDD vs Shizhong   p <= 0.1    r = 0.7812
SCZD vs PRSice    p <= 0.05   r = 0.8255
SCZD vs Shizhong  p <= 0.001  r = 0.7831
```

Interpretation: GraphPred tracks the broad C+T PRS signal well, but less tightly than PRSice and Shizhong track each other. This is expected because the methods use LD differently. PRSice and Shizhong use LD for clumping/variant selection before scoring retained marginal GWAS effects. GraphPred/LDGM uses LD structure to estimate adjusted genome-wide loadings before `bcftools +score` applies them to target dosages.

### Diagnosis Group Patterns

Using z-scored values, Bipolar samples tend to have higher average PRS than controls for several selected comparisons, including GraphPred and broad C+T thresholds. For example, at `p <= 0.1`:

- BPD PRSice BPD score mean z:
  - Control: `-0.1096`
  - MDD: `-0.2059`
  - Bipolar: `0.3104`
- BPD Shizhong BPD score mean z:
  - Control: `-0.0681`
  - MDD: `-0.2347`
  - Bipolar: `0.2969`

GraphPred shows a similar descriptive pattern for Bipolar samples:

- GraphPred BPD mean z:
  - Control: `-0.1523`
  - MDD: `-0.1408`
  - Bipolar: `0.2896`

These group summaries are descriptive. They should not be presented as confirmatory association results without a separate statistical analysis and multiple-testing strategy.

## Detailed Conclusions

1. The inputs are harmonized.
   - All three score sources cover the same 119 MBv samples.
   - No samples were lost during comparison.

2. PRSice and Shizhong are highly concordant C+T implementations.
   - Pearson correlations are generally near or above 0.90 at shared thresholds.
   - This supports that both pipelines are capturing the same broad PRS signal.

3. PRSice and Shizhong are not interchangeable.
   - The threshold grids differ.
   - Clumping settings and implementation details differ.
   - Direct comparisons should be limited to exact shared thresholds.

4. GraphPred is related but methodologically distinct.
   - Its correlations with C+T scores are high but lower than PRSice-vs-Shizhong.
   - This supports using GraphPred as a complementary LD-model-based score rather than treating it as another p-thresholded C+T score.

5. Threshold choice matters.
   - Best-correlated thresholds differ by disorder and method.
   - This reinforces that internal threshold selection in `n = 119` should be treated as exploratory.

6. For presentation, use:
   - `matched_threshold_correlations.tsv` for PRSice vs Shizhong.
   - `graphpred_threshold_correlations.tsv` for GraphPred vs thresholded methods.
   - `threshold_availability.tsv` to explain threshold mismatch.
   - `diagnosis_group_summary.tsv` for descriptive group plots or tables.

## Recommended Wording

The Shizhong and PRSice pipelines produce strongly concordant LD-aware clumping-and-thresholding PRS across the same MBv samples, especially at shared p-value thresholds. GraphPred/bcftools scores are also substantially correlated with both C+T approaches, but less tightly, consistent with a distinct LD-model-based genome-wide weighting strategy. Therefore, the three score sources are mutually reassuring but should not be treated as identical; PRSice and Shizhong should be compared only at exact shared thresholds, and GraphPred should be reported as a complementary LDGM/GraphPred-weighted score.
