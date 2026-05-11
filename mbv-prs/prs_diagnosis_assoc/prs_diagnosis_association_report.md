# PRS Diagnosis Association Report

Samples: 119 total; 39 MDD, 40 Bipolar/BPD, 40 Controls. No SCZ diagnosis group is present, so SCZD PRS is evaluated only as a cross-disorder score.

Models: unadjusted logistic regression and sensitivity logistic regression adjusted for Sex plus PC1-PC5. Each PRS was standardized over all 119 samples. PRSice best thresholds are internal exploratory flags, not confirmatory choices.

## Best Internal PRSice Thresholds

             contrast                model score_disorder threshold n_case
       BPD_vs_Control adjusted_sex_PC1_PC5            BPD     5e-08     40
       BPD_vs_Control adjusted_sex_PC1_PC5            MDD    0.0001     40
       BPD_vs_Control adjusted_sex_PC1_PC5           SCZD      0.05     40
       BPD_vs_Control           unadjusted            BPD     5e-08     40
       BPD_vs_Control           unadjusted            MDD    0.0001     40
       BPD_vs_Control           unadjusted           SCZD      0.05     40
           BPD_vs_MDD adjusted_sex_PC1_PC5            BPD     1e-06     40
           BPD_vs_MDD adjusted_sex_PC1_PC5            MDD       0.1     40
           BPD_vs_MDD adjusted_sex_PC1_PC5           SCZD     0.001     40
           BPD_vs_MDD           unadjusted            BPD     1e-06     40
           BPD_vs_MDD           unadjusted            MDD       0.1     40
           BPD_vs_MDD           unadjusted           SCZD     0.001     40
 CaseUnion_vs_Control adjusted_sex_PC1_PC5            BPD     5e-08     79
 CaseUnion_vs_Control adjusted_sex_PC1_PC5            MDD    0.0001     79
 CaseUnion_vs_Control adjusted_sex_PC1_PC5           SCZD       0.1     79
 CaseUnion_vs_Control           unadjusted            BPD     5e-08     79
 CaseUnion_vs_Control           unadjusted            MDD    0.0001     79
 CaseUnion_vs_Control           unadjusted           SCZD       0.1     79
       MDD_vs_Control adjusted_sex_PC1_PC5            BPD     5e-08     39
       MDD_vs_Control adjusted_sex_PC1_PC5            MDD     5e-08     39
       MDD_vs_Control adjusted_sex_PC1_PC5           SCZD       0.1     39
       MDD_vs_Control           unadjusted            BPD     5e-08     39
       MDD_vs_Control           unadjusted            MDD     5e-08     39
       MDD_vs_Control           unadjusted           SCZD       0.1     39
 n_control       AUC        OR    OR_low  OR_high           p
        40 0.6343750 2.0834840 1.0992305 3.949040 0.024448809
        40 0.6918750 2.1698531 1.2315151 3.823146 0.007348932
        40 0.6762500 2.3455768 1.3219357 4.161875 0.003568816
        40 0.6343750 1.9518476 1.0690175 3.563748 0.029459680
        40 0.6918750 2.0285817 1.1933699 3.448339 0.008973238
        40 0.6762500 2.0672221 1.2321663 3.468206 0.005944587
        39 0.6788462 2.0853001 1.2105820 3.592055 0.008078649
        39 0.6750000 2.0469607 1.1990124 3.494583 0.008662227
        39 0.6339744 1.9711363 1.1405342 3.406630 0.015053939
        39 0.6788462 1.9882550 1.1839872 3.338852 0.009361055
        39 0.6750000 1.9735796 1.1923772 3.266597 0.008184028
        39 0.6339744 1.9062714 1.1197825 3.245158 0.017463565
        40 0.5829114 1.2944785 0.8692711 1.927678 0.203935191
        40 0.6142405 1.4650678 0.9597907 2.236346 0.076758482
        40 0.6170886 1.7453970 1.1211135 2.717308 0.013655722
        40 0.5829114 1.2336506 0.8464902 1.797887 0.274516237
        40 0.6142405 1.4363571 0.9581539 2.153226 0.079592858
        40 0.6170886 1.6156878 1.0684065 2.443309 0.022991915
        40 0.5301282 1.0366547 0.6811238 1.577765 0.866590433
        40 0.5532051 1.1859837 0.7666725 1.834626 0.443484303
        40 0.5660256 1.4411436 0.8848852 2.347078 0.141958099
        40 0.5301282 0.9885425 0.6624061 1.475252 0.955010341
        40 0.5532051 1.1049180 0.7472895 1.633696 0.617048999
        40 0.5660256 1.3351040 0.8427050 2.115215 0.218309498
 q_fdr_within_contrast_model     cohen_d
                  0.06963326  0.52116687
                  0.05187527  0.64034926
                  0.05187527  0.67794832
                  0.06791291  0.52116687
                  0.05730266  0.64034926
                  0.05730266  0.67794832
                  0.06101724  0.64570464
                  0.06101724  0.64496196
                  0.06101724  0.57142246
                  0.07052869  0.64570464
                  0.07052869  0.64496196
                  0.07052869  0.57142246
                  0.40702743  0.21287448
                  0.38379241  0.34640959
                  0.25520998  0.45731044
                  0.45752706  0.21287448
                  0.42051419  0.34640959
                  0.30253710  0.45731044
                  0.92848975 -0.01253443
                  0.92848975  0.11144138
                  0.92848975  0.27858993
                  0.95903944 -0.01253443
                  0.95903944  0.11144138
                  0.95903944  0.27858993

## Top Unadjusted Scores by AUC

             contrast                             score          method
       BPD_vs_Control              PRSice_MDD_Pt_0.0001          PRSice
       BPD_vs_Control      MDD_MDD_2025_pgs_a0.5_b2e-08 GraphPred_score
       BPD_vs_Control               PRSice_SCZD_Pt_0.05          PRSice
       BPD_vs_Control SCZD_SCZ_2022.EUR_pgs_a0.5_b2e-07 GraphPred_score
       BPD_vs_Control                PRSice_SCZD_Pt_0.1          PRSice
           BPD_vs_MDD               PRSice_BPD_Pt_1e-06          PRSice
           BPD_vs_MDD                 PRSice_MDD_Pt_0.1          PRSice
           BPD_vs_MDD                 PRSice_MDD_Pt_0.5          PRSice
           BPD_vs_MDD                   PRSice_MDD_Pt_1          PRSice
           BPD_vs_MDD              PRSice_MDD_Pt_0.0001          PRSice
 CaseUnion_vs_Control      MDD_MDD_2025_pgs_a0.5_b2e-08 GraphPred_score
 CaseUnion_vs_Control                PRSice_SCZD_Pt_0.1          PRSice
 CaseUnion_vs_Control              PRSice_MDD_Pt_0.0001          PRSice
 CaseUnion_vs_Control               PRSice_SCZD_Pt_0.05          PRSice
 CaseUnion_vs_Control SCZD_SCZ_2022.EUR_pgs_a0.5_b2e-07 GraphPred_score
       MDD_vs_Control                PRSice_SCZD_Pt_0.1          PRSice
       MDD_vs_Control               PRSice_SCZD_Pt_0.01          PRSice
       MDD_vs_Control               PRSice_MDD_Pt_5e-08          PRSice
       MDD_vs_Control      MDD_MDD_2025_pgs_a0.5_b2e-08 GraphPred_score
       MDD_vs_Control               PRSice_SCZD_Pt_0.05          PRSice
 score_disorder threshold       AUC       OR           p   cohen_d
            MDD    0.0001 0.6918750 2.028582 0.008973238 0.6403493
            MDD GraphPred 0.6868750 2.100981 0.005357285 0.6879584
           SCZD      0.05 0.6762500 2.067222 0.005944587 0.6779483
           SCZD GraphPred 0.6675000 1.932489 0.009587442 0.6237062
           SCZD       0.1 0.6668750 1.943616 0.009750428 0.6288069
            BPD     1e-06 0.6788462 1.988255 0.009361055 0.6457046
            MDD       0.1 0.6750000 1.973580 0.008184028 0.6449620
            MDD       0.5 0.6737179 1.957812 0.009716835 0.6249652
            MDD         1 0.6628205 1.893555 0.013491367 0.5921649
            MDD    0.0001 0.6583333 1.907194 0.015483158 0.5820144
            MDD GraphPred 0.6180380 1.632625 0.020575940 0.4677832
           SCZD       0.1 0.6170886 1.615688 0.022991915 0.4573104
            MDD    0.0001 0.6142405 1.436357 0.079592858 0.3464096
           SCZD      0.05 0.6117089 1.570249 0.030253710 0.4336162
           SCZD GraphPred 0.5984177 1.413028 0.084102837 0.3407501
           SCZD       0.1 0.5660256 1.335104 0.218309498 0.2785899
           SCZD      0.01 0.5538462 1.259554 0.298088092 0.2342472
            MDD     5e-08 0.5532051 1.104918 0.617048999 0.1114414
            MDD GraphPred 0.5474359 1.329348 0.237289040 0.2674444
           SCZD      0.05 0.5455128 1.230149 0.376627050 0.1982004

## GWAS Significant Variant Overlap Context

                                   metric        value
                 BPD_significant_variants 3.017000e+03
                 MDD_significant_variants 1.444800e+04
      BPD_MDD_shared_significant_variants 8.990000e+02
                          BPD_MDD_jaccard 5.426778e-02
 fraction_BPD_significant_shared_with_MDD 2.979781e-01
 fraction_MDD_significant_shared_with_BPD 6.222315e-02

MDD-BPD genome-wide significant exact-variant overlap is modest by Jaccard, but asymmetric by disorder because BPD has fewer significant variants. This overlap is not equivalent to PRS overlap because PRSice thresholds include sub-significant variants and GraphPred uses LD-aware genome-wide weights.

## Threshold Selection Guidance

Do not select a PRSice threshold for confirmatory downstream analysis solely because it separates MBv cases and controls best. This is overfit at n=119. Prefer a pre-specified score, an external validation cohort, cross-validation if sample size allows, or PRSice permutation/empirical p-values for association. Report all thresholds here and label internally best thresholds as exploratory.

## Output Files

- `prs_association_all_scores.csv`
- `prs_best_thresholds_exploratory.csv`
- `prs_group_summary.csv`
- `prs_overlap_interpretation.csv`
- `prs_auc_by_threshold.svg`
- `prs_selected_distributions.svg`
