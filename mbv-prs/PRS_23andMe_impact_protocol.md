# PRS Sensitivity Protocol for Added GWAS Cohorts

## Purpose

This protocol recalculates polygenic scores from two versions of the same GWAS:

- a public summary-statistics release excluding a restricted cohort; and
- an integrated release that adds that independent cohort by fixed-effect
  inverse-variance meta-analysis.

The comparison measures how added discovery sample size changes scores and donor
ranks. It does not select a p-value threshold or test phenotype association in
the scoring cohort.

Reusable scripts are tracked in `scripts/`. Cohort manifests, genotype targets,
scores, figures, logs, and reports belong under an ignored `generated/`
directory. Large GWAS and LD-reference files remain in external reference
storage.

## Canonical Base Contract

All C+T and PRSice runs consume a seven-column tab-delimited base table:

```text
CHR BP SNP A1 A2 BETA P
```

`SNP` must be `chr:pos:REF:ALT`, `A1` must be ALT, `A2` must be REF, and `BETA`
must be the ALT-allele log-odds effect. `scripts/validate_prs_base.pl` rejects:

- non-autosomal rows;
- invalid or identical alleles;
- a canonical ID that disagrees with `A2/A1`;
- duplicate canonical IDs;
- non-finite beta or p values; and
- p values outside `0 < P <= 1`.

`scripts/prs_base_from_gwas_bcf.sh` creates this table from a single-sample
GWAS-VCF BCF, using `FORMAT/ES` and `FORMAT/LP`. It retains PASS or unfiltered
autosomal biallelic SNPs and interprets ES as the ALT-allele beta.

## LD Reference

The primary LD reference is Shizhong's original 489-person European founder
panel: 8,670,255 hg38 SNPs after no-indel and MAF >= 0.01 filtering. The newer
633-person panel is not used for primary C+T because it includes 108
non-founders. All 489 historical samples are present in the expanded panel; the
difference follows from the expanded family-based 1000 Genomes release and
different variant filters. See the [IGSR expanded high-coverage release](https://www.internationalgenome.org/announcements/3202-samples-at-high-coverage-from-NYGC/)
and [IGSR related-sample guidance](https://www.internationalgenome.org/faq/which-datasets-include-related-individuals/).

Store the 22 PLINK BED/BIM/FAM triplets outside Git. First build one combined
prefix with:

```bash
scripts/merge_plink_chromosomes.sh \
  --pattern '/dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong/EUR_noIndels_chr{CHR}_maf0.01_hg38' \
  --out /dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong/EUR_noIndels_all_maf0.01_hg38
```

The merge script requires identical FAM files, checks that the merged BIM count
equals the sum of all chromosome BIM counts, and records SHA-256 checksums.

PRSice joins base, target, and LD-reference variants by ID. If the original LD
panel uses rsIDs while the base and target use canonical IDs, create a
target-overlapping canonical-ID copy without changing LD genotypes or samples:

```bash
scripts/build_prsice_ld_panel.sh \
  --source-prefix /dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong/EUR_noIndels_all_maf0.01_hg38 \
  --ld-pattern '/dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong/EUR_noIndels_chr{CHR}_maf0.01_hg38' \
  --target-pvar generated/prs_23andme_impact/target/mbv_qc.pvar \
  --out-prefix /dbdata/cdb/ref/1000g_hg38/1KG_EUR_489_shizhong/EUR_noIndels_all_maf0.01_hg38.mbv_qc_canonical
```

The canonical copy is target-specific reference infrastructure and remains
outside Git. Allele mismatches and ambiguous target mappings are excluded and
counted.

## Methods

### Shizhong-Style C+T

`scripts/run_plink_ct_manifest.sh` reproduces the legacy two-tool method:

1. Match base, LD-reference, and target variants by chromosome, position, and
   exact unordered allele pair.
2. Use LD-reference IDs for PLINK 1.9 clumping.
3. Map selected variants back to target IDs and use base `A1/BETA` for scoring.
4. Use PLINK2 `--q-score-range` and retain `SCORE1_AVG`.

Primary settings are `--clump-kb 1000`, `--clump-r2 0.1`, `p1=1`, and `p2=1`.
The historical grid is:

```text
1e-08,1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.05,0.1,0.2,
0.3,0.4,0.5,0.6,0.7,0.8,0.9,1
```

Mapping statistics distinguish absent positions, allele mismatches, and
ambiguous duplicate mappings for every chromosome. Nonmatching alleles are not
silently flipped or scored.

### PRSice

`scripts/run_prsice_manifest.sh` uses the harmonized target, the combined
489-person LD prefix, `250 kb`, `r2=0.1`, `--beta`, `--score sum`, and:

```text
5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1
```

PRSice sums and Shizhong `SCORE1_AVG` are different raw scales. Correlation and
within-series rank comparisons remain valid; raw values should not be compared
as if they were calibrated to the same units.

### LDGM/GraphPred

`scripts/run_graphpred_manifest.sh` applies existing full and no-23andMe PGS
BCFs to the same final-QC dosage BCF using:

```bash
bcftools +score --use DS --sample-header SAMPLE --counts
```

This stage does not rerun `bcftools +pgs`. It records the plugin-reported matched
variant count with each score.

## Manifest Contracts

C+T:

```text
RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_PGEN LD_PATTERN OUT_DIR CLUMP_KB CLUMP_R2 THRESHOLDS
```

PRSice:

```text
RUN_ID TRAIT GWAS_VERSION SETUP BASE TARGET_BGEN LD_PREFIX PHENO OUT_PREFIX CLUMP_KB CLUMP_R2 THRESHOLDS
```

GraphPred:

```text
RUN_ID TRAIT GWAS_VERSION SETUP TARGET_BCF PGS_BCF OUT_FILE
```

Fields are tab-separated. `LD_PATTERN` must contain `{CHR}`.

## MBv Execution

The cohort-specific driver writes its manifests and all work products under
`generated/prs_23andme_impact/`:

```bash
scripts/run_mbv_prs_23andme_impact.sh all
```

Stages can be resumed separately: `prepare`, `smoke`, `ct`, `prsice`,
`graphpred`, `collect`, `validate`, and `compare`. The smoke stage runs
chromosome 22 for both the historical and harmonized targets before production
C+T starts.

The long registry contains:

```text
SAMPLE TRAIT GWAS_VERSION METHOD SETUP THRESHOLD RAW_SCORE SCORED_VARIANTS
```

Every series must contain exactly the same 119 unique donors. The comparison
script reports Pearson and Spearman correlations, z-score RMSE, donor rank
changes, top/bottom-decile overlap, scored-variant counts, exact-threshold C+T
method correlations, and GraphPred-to-threshold correlations.

## Required Validation

- Run shell, Perl, and R syntax checks and `tests/test_prs_workflow.sh`.
- Require all chromosome status files to contain no failures.
- Require historical no-23andMe Shizhong scores to correlate at least 0.999
  with archived scores at every threshold.
- Require no-23andMe GraphPred scores to reproduce archived scores within
  `1e-6` maximum absolute error.
- Retain software versions, input checksums, base validation counts, mapping
  counts, clumped/scored counts, and sample-list checks.
- Confirm `generated/` and external reference data are ignored or outside Git.

## BD Qualification

The reconstructed full BD base and PGS inputs are usable for sensitivity
analysis but are internally identified as pre-DENTIST. The v7.0 23andMe
association file was paired with supplied v7.2 European annotations, and the
paper's final HRC-based post-meta-analysis DENTIST filter has not been rerun.
Do not describe this reconstruction as the exact final paper summary statistics.
See `BD-DENTIST-issue-v7.0-vs-v7.2.md`.
