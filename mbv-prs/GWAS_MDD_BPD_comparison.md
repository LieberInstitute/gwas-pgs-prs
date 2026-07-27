# GWAS MDD vs BPD cohort-size comparison

## Question

Compare the BPD and MDD GWAS datasets used for EUR PRS work, focusing on:

- cohort sizes with and without 23andMe
- proportion of discovery power contributed by 23andMe
- whether missing 23andMe data plausibly explains the lower BPD significant-variant count

## Main EUR comparison

`N_eff/2` is the half effective sample size used in the papers/tables.

| Study | EUR dataset | Cases | Controls | N_eff/2 | 23andMe cases | 23andMe controls | 23andMe % cases | 23andMe % controls | 23andMe % N_eff/2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BPD | no 23andMe | 59,287 | 781,022 | 81,649 | NA | NA | NA | NA | NA |
| BPD | with 23andMe | 131,969 | 2,322,416 | 220,467 | 72,682 | 1,541,394 | 55.1% | 66.4% | 63.0% |
| MDD | no 23andMe | 412,305 | 1,588,397 | 576,327 | NA | NA | NA | NA | NA |
| MDD | with 23andMe | 525,197 | 3,362,335 | 788,603 | 112,892 | 1,773,938 | 21.5% | 52.8% | 26.9% |

## Power gain from 23andMe

| Study | EUR N_eff/2 no 23andMe | EUR N_eff/2 with 23andMe | Fold gain from 23andMe |
|---|---:|---:|---:|
| BPD | 81,649 | 220,467 | 2.70x |
| MDD | 576,327 | 788,603 | 1.37x |

Interpretation:

- 23andMe contributes much more proportional EUR discovery power to BPD than to MDD.
- Even without 23andMe, MDD is much larger than BPD: 412,305 vs 59,287 EUR cases, and 576,327 vs 81,649 EUR N_eff/2.
- The BPD no-23andMe EUR file therefore has about 7.0x fewer cases and 7.1x lower N_eff/2 than the MDD no-23andMe EUR file.

## Significant discovery counts from papers

These are independent signals/loci from the papers, not raw significant variant rows from local GWAS-VCF files.

| Study | Analysis | Loci or regions | Independent GWS signals/SNPs |
|---|---|---:|---:|
| BPD | EUR clinical + community, no self-report/23andMe | 88 | 94 |
| BPD | EUR clinical + community + self-report/23andMe | 229 | 261 |
| MDD | EUR full analysis with 23andMe | 570 | 622 |
| MDD | all-ancestry full analysis with 23andMe | 635 | 697 |

Note: I did not find a published MDD no-23andMe independent-locus count in the checked sources. The public no-23andMe MDD summary-stat release provides cohort sizes and summary statistics, but the paper reports the main EUR locus count for the full EUR analysis.

## All-ancestry check

This is not the primary PRS input, but it confirms the same pattern.

| Study | Dataset | Cases | Controls | N_eff/2 | 23andMe cases | 23andMe controls | 23andMe % cases | 23andMe % N_eff/2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BPD | no 23andMe | 67,948 | 867,710 | 95,861 | NA | NA | NA | NA |
| BPD | with 23andMe | 158,036 | 2,796,499 | 267,860 | 90,088 | 1,928,789 | 57.0% | 64.2% |
| MDD | no 23andMe | 547,355 | 2,071,918 | 739,180 | NA | NA | NA | NA |
| MDD | with 23andMe | 688,808 | 4,364,225 | 1,004,459 | 141,453 | 2,292,307 | 20.5% | 26.4% |

## Local raw significant row context

From `notes_stats_comments.md`, using `P <= 5e-8` / `LP >= 7.30103`:

| Disorder | Finite LP rows | Significant raw variant rows | Max LP |
|---|---:|---:|---:|
| BPD | 6,938,758 | 3,017 | 16.783 |
| MDD | 7,362,674 | 14,448 | 26.6519 |
| SCZD | 7,658,367 | 20,457 | 39.184 |

These are raw variant rows, not LD-clumped independent loci. They should not be compared directly with paper locus counts.

## Data provenance and calculations

- BPD EUR with and without 23andMe comes directly from BPD Supplementary Table 2.
- BPD 23andMe EUR contribution is the Supplementary Table 2 self-report block: 72,682 cases and 1,541,394 controls.
- BPD public no-23andMe EUR cohort metadata from Figshare matches the no-23andMe EUR counts: 59,287 cases and 781,022 controls.
- MDD full EUR counts come from the MDD paper Table 1: 525,197 cases, 3,362,335 controls, N_eff/2 = 788,603.
- MDD public no-23andMe EUR counts come from the Figshare cohort metadata file `pgc-mdd2025_no23andMe_eur_v3.49.24.11.txt`: 412,305 cases, 1,588,397 controls.
- MDD no-23andMe EUR N_eff/2 = 576,327 comes from the MDD paper/ScienceDirect text for PGS training data excluding 23andMe.
- MDD 23andMe contribution above is calculated as full EUR minus public no-23andMe EUR: 112,892 cases, 1,773,938 controls, and 212,276 N_eff/2.

## Bottom line

For EUR analyses, missing 23andMe data is much more consequential for BPD than for MDD. Adding 23andMe would increase BPD EUR N_eff/2 by about 2.70x, but MDD EUR N_eff/2 by about 1.37x. However, even after adding 23andMe, BPD EUR N_eff/2 remains much smaller than MDD EUR N_eff/2: 220,467 vs 788,603.

This supports the low BPD significant-variant count as primarily a discovery-power issue in the public no-23andMe BPD EUR release, not a parsing failure.

## Sources

- BPD paper: O'Connell et al. Genomics yields biological and phenotypic insights into bipolar disorder. Nature 639, 968-975 (2025). https://doi.org/10.1038/s41586-024-08468-9
- BPD Supplementary Tables 1-35: linked from the Nature paper, Supplementary Tables download. https://www.nature.com/articles/s41586-024-08468-9
- BPD public summary statistics: https://figshare.com/articles/dataset/bip2024/27216117
- MDD paper: Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium. Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies. Cell 188, 640-652.e9 (2025). https://doi.org/10.1016/j.cell.2024.12.002
- MDD public no-23andMe summary statistics: https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255
- MDD ScienceDirect article page, including PGS training effective sample sizes excluding 23andMe: https://www.sciencedirect.com/science/article/pii/S0092867424014156

## BIP 2021 vs BIP 2024 EUR comparison

This table compares the older BIP 2021 GWAS with the 2024 EUR BPD analyses. The 2021 values come from Supplementary Table 1 total row and Supplementary Table 2 locus rows. The 2024 values come from Supplementary Table S1 and Supplementary Table S2.

| BPD GWAS | EUR dataset | Cases | Controls | N_eff/2 | Significant loci |
|---|---|---:|---:|---:|---:|
| BIP 2021 | PGC3 EUR | 41,917 | 371,549 | 50,981 | 64 |
| BIP 2024 | EUR clinical + community, no self-report/23andMe | 59,287 | 781,022 | 81,649 | 88 |
| BIP 2024 | EUR clinical + community + self-report/23andMe | 131,969 | 2,322,416 | 220,467 | 229 |

Checked interpretation:

- BIP 2021 did not have more cases or larger EUR N_eff/2 than the 2024 no-23andMe EUR dataset used here.
- BIP 2021 had 17,370 fewer EUR cases, 409,473 fewer controls, and lower N_eff/2 than BIP 2024 no-23andMe.
- The `229` loci value is correct for the BIP 2024 EUR analysis with self-report/23andMe.
- The larger BIP 2024 number sometimes cited is `298` loci, which is for the full multi-ancestry analysis with self-report/23andMe. The corresponding independent GWS signal count is `337`.
- For BIP 2021, use the reported Supplementary Table 1 `N eff half = 50,981`. The aggregate formula from total cases and controls gives a larger value, but the paper reports the summed cohort-level effective sample size.

Additional BIP 2021 source:

- BIP 2021 paper: Mullins et al. Genome-wide association study of more than 40,000 bipolar disorder cases provides new insights into the underlying biology. Nature Genetics 53, 817-829 (2021). https://doi.org/10.1038/s41588-021-00857-4
- BIP 2021 Supplementary Tables: linked from the Nature Genetics paper, including Supplementary Table 1 and Supplementary Table 2. https://www.nature.com/articles/s41588-021-00857-4
