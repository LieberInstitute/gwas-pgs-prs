# Why Use bcftools GraphPred Scores Here?

## Short answer

For this MBv project, `bcftools +score` is best understood as a fast, exact VCF/BCF scoring engine. It applies already-computed effect weights to target genotype dosages. The GraphPred part is upstream: `bcftools +pgs` used GWAS-VCF summary statistics plus EUR LDGM precision matrices to estimate LD-aware loadings, and `bcftools +score --use DS` applied those loadings to the 119 MBv samples.

This makes the bcftools approach a good primary score for this project because it:

- uses imputed dosage (`DS`) directly from VCF/BCF, without lossy conversion for scoring;
- matches variants by exact genomic position and alleles in the GWAS-VCF/PGS BCF representation;
- uses millions of LD-aware weights instead of a small clumped subset;
- avoids choosing a best PRSice p-value threshold inside a small target cohort;
- fits psychiatric traits, where signal is expected to be highly polygenic.

PRSice remains a strong and popular comparator. The point is not that PRSice is obsolete. The point is that PRSice implements clumping and thresholding (C+T), while the bcftools GraphPred workflow applies LD-aware genome-wide weights before target scoring.

## What bcftools +score actually computes

A polygenic score is usually a weighted allele count or dosage sum. The PGS Catalog describes this general model as the sum of genotype dosages multiplied by weights across variants: https://www.pgscatalog.org/

`bcftools +score` implements that scoring step. In VCF/BCF mode, it reads effect sizes from `FORMAT/ES` in the score BCF and, for each target sample, adds:

```text
score_sample = sum over matched variants (ES_variant * dosage_sample,effect_allele)
```

In this project, the command used:

```bash
bcftools +score --use DS --sample-header --counts \
  -o scores_bcftools/<DISORDER>.graphpred.tsv \
  qc/mbv.qc.chr.bcf <DISORDER_PGS_BCF>
```

With `--use DS`, the plugin uses target imputed dosage. For a biallelic variant, target `DS` is the ALT dosage, REF dosage is `2 - DS`, and the plugin adds the dosage for the effect allele named in the score record. Missing sample-site dosage values are skipped. The output score is therefore a direct dosage-weighted sum over matched variants.

Important precision: `bcftools +score` is not itself GraphPred and does not estimate LD-adjusted effects. It applies whatever weights are supplied. In this project those supplied weights are GraphPred/LDGM PGS BCF loadings.

## What GraphPred / LDGM means here

The upstream tool is `bcftools +pgs`, documented in the `freeseek/score` repository: https://github.com/freeseek/score

That documentation describes `+pgs` as inspired by the GraphPred algorithm. The workflow starts from GWAS-VCF summary statistics and LDGM-VCF precision matrices, then computes LD-aware PGS loadings. The documented model combines:

- a BLUP-style LD-adjusted component;
- a sparse-effect refinement inspired by SuSiE, implemented with Gibbs sampling;
- sparse LD graphical model (LDGM) precision matrices for efficient LD computations.

The strongest peer-reviewed base for this stack is the LDGM work. Nowbandegani et al. introduced LD graphical models for efficient inference with summary statistics and published BLUPx-ldgm in Nature Genetics: https://doi.org/10.1038/s41588-023-01487-8

The bcftools `+pgs` implementation is newer than PRSice and should be presented as a documented software implementation built on LDGM ideas, not as a method with the same citation depth as PRSice. The Broad software page also distributes related bcftools score binaries and LDGM resources: https://software.broadinstitute.org/software/score/

## How this differs from PRSice

PRSice-2 is a widely used C+T tool. It takes GWAS marginal effects, clumps variants using an LD reference to reduce correlated markers, then computes scores at one or more p-value thresholds. The PRSice-2 paper is Choi and O'Reilly 2019: https://doi.org/10.1093/gigascience/giz082

In this MBv project, PRSice was run as a fixed-threshold comparator:

- base files came from converted GWAS BCFs;
- the LD reference was the EUR PLINK panel at `/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos`;
- thresholds were `5e-8, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.5, 1`;
- no target-sample regression optimization was used.

Conceptual difference:

- PRSice C+T reduces LD by selecting approximately independent variants and thresholding by GWAS p-value.
- GraphPred/LDGM scoring uses LD-adjusted loadings across many more variants, then `+score` applies those loadings directly to target dosages.

This is similar in motivation to other LD-aware summary-statistic PRS methods such as PRS-CS, which uses a Bayesian regression framework with continuous shrinkage priors: https://doi.org/10.1038/s41467-019-09718-5

Best-practice PRS guidance still recognizes C+T as a common baseline and emphasizes careful ancestry matching, external validation, and avoiding overfitting threshold choice in the target data: https://doi.org/10.1038/s41596-020-0353-1

## Why this project uses bcftools GraphPred scores as primary

1. Direct use of MBv imputed dosages

`bcftools +score --use DS` scores directly from `qc/mbv.qc.chr.bcf`. This avoids a scoring-time conversion to PLINK/BGEN formats. In this project, preserving chr-prefixed contigs mattered because PLINK2 BCF export changed `chr1` style names to `1`, which broke exact matching to GWAS-VCF/PGS BCF files.

2. Exact allele/build matching

The bcftools path uses normalized VCF/BCF records and exact CHROM/POS/REF/ALT matching. This is especially useful when the score files are already GWAS-VCF/PGS BCFs with `ES` relative to ALT alleles.

3. More genome-wide information

The GraphPred/LDGM scores matched millions of target variants:

```text
BD: 4996043
MDD: 4832784
SCZD: 5128978
```

PRSice C+T intentionally retains a much smaller clumped subset at each threshold. That can be robust and interpretable, but it discards many correlated variants whose information an LD-aware model can use.

4. Avoiding threshold selection inside n=119

With only 119 MBv samples, selecting the "best" PRSice threshold by observed case-control separation would be exploratory and likely overfit. A precomputed LD-aware GraphPred score gives one primary score per disorder without choosing a threshold from the target cohort.

5. Psychiatric traits are highly polygenic

For BD, MDD, and SCZD, many small effects across the genome are expected. LD-aware genome-wide methods are designed for this setting. PRSice remains useful as a comparator and sanity check.

## Caveats and what not to claim

Avoid these claims:

- Treating `bcftools +score` itself as the GraphPred model. It is the scoring engine. `bcftools +pgs` generated the GraphPred/LDGM loadings.
- Treating GraphPred as guaranteed to outperform PRSice. Performance is trait, ancestry, sample, GWAS, LD reference, and validation-set dependent.
- Treating the PRSice thresholds as proof of a best model here. At n=119, internal threshold ranking is exploratory unless externally validated or corrected.
- Treating LDGM as a solution to all ancestry issues. LD reference ancestry still matters. The current loadings and PRSice LD reference are EUR-oriented.

Reasonable claim:

```text
For this MBv analysis, the bcftools GraphPred workflow is a well-justified primary PRS because it applies LD-aware GWAS-derived loadings directly to imputed MBv dosages with exact VCF allele matching, while PRSice provides a popular C+T comparator.
```

## Project-specific evidence from the no-23andMe run

These results describe the completed MBv run with the archived no-23andMe BD
and MDD loadings. They are historical evidence for the scoring implementation,
not results from the new full-European loadings.

The run reported:

- input samples: 119;
- final samples after QC: 119;
- final variants: 6639276;
- DS out-of-range sites: 0;
- KING related pairs >= 0.0884: 0;
- GraphPred overlaps: BD 4996043, MDD 4832784, SCZD 5128978.

Method correlations between GraphPred and PRSice were moderate to high, depending on trait and threshold. For example, Pearson correlations with PRSice threshold 0.05 were:

```text
BD: 0.814091266694973
MDD: 0.792858290899503
SCZD: 0.825026915058488
```

This is reassuring because both methods are measuring related GWAS signal. They are not identical because the models use different LD handling and variant inclusion rules.

## Sources

- PGS Catalog, general PGS definition and catalog context: https://www.pgscatalog.org/
- BCFtools score/Broad software page: https://software.broadinstitute.org/software/score/
- `freeseek/score` documentation and source for `bcftools +score`, `+pgs`, GWAS-VCF, and LDGM-VCF: https://github.com/freeseek/score
- LDGM and BLUPx-ldgm: Nowbandegani et al. 2023, Nature Genetics, https://doi.org/10.1038/s41588-023-01487-8
- PRSice-2: Choi and O'Reilly 2019, GigaScience, https://doi.org/10.1093/gigascience/giz082
- PRS-CS: Ge et al. 2019, Nature Communications, https://doi.org/10.1038/s41467-019-09718-5
- PRS tutorial and best-practice context: Choi et al. 2020, Nature Protocols, https://doi.org/10.1038/s41596-020-0353-1
- Recent PGS method benchmarking context: PGS-hub 2026, Nature Communications, https://doi.org/10.1038/s41467-026-68599-7
