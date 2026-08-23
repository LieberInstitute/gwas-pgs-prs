# Integrating the public and 23andMe GWAS components

## Short answer

For each disorder, the files are two independent GWAS results:

1. the public European meta-analysis excluding 23andMe; and
2. the European 23andMe-only GWAS.

They must not be concatenated. For each shared, allele-aligned variant, the two
log-odds estimates are combined by standard-error inverse-variance fixed-effect
meta-analysis. The combined beta, standard error, odds ratio, and P value are
new values. No individual-level genotypes are needed.

The implementation is
[`scripts/integrate_23andme_gwas.sh`](scripts/integrate_23andme_gwas.sh).
It calls the local `bcftools +metal` implementation of the METAL method rather
than reimplementing the arithmetic.

## What is recalculated

For public component `p` and 23andMe component `m`:

```text
w_p = 1 / SE_p^2
w_m = 1 / SE_m^2

BETA_full = (w_p * BETA_p + w_m * BETA_m) / (w_p + w_m)
SE_full   = sqrt(1 / (w_p + w_m))
Z_full    = BETA_full / SE_full
P_full    = 2 * Phi(-abs(Z_full))
OR_full   = exp(BETA_full)
```

The public BD odds ratio is first converted to `log(OR)`. MDD already supplies
beta. The 23andMe `effect` is already log odds. Effects are reoriented to the
same genomic allele before applying the formula.

The delivered 23andMe P value is not averaged with the public P value. The
final P value follows from the combined beta and SE. The 23andMe P value is a
likelihood-ratio-test result, whereas standard-error meta-analysis uses beta
and SE.

## Why this is the appropriate method

This is the standard METAL inverse-variance method. It is valid here because:

- both components estimate log odds for the same binary phenotype;
- their standard errors are available;
- the public files explicitly exclude 23andMe, so the samples are independent;
- fixed-effect meta-analysis is associative, so a previous no-23andMe
  meta-analysis can be combined with one independent 23andMe component; and
- both source papers used standard-error inverse-weighted fixed-effect
  meta-analysis.

The METAL paper defines the same weights and combined beta/SE calculation.
GWAS meta-analysis QC guidance additionally requires consistent alleles,
effect scales, sample definitions, variant QC, and effective sample size.

## Important finding about genomic inflation

The raw 23andMe files say their association results are not adjusted for
inflation. Their HTML reports show adjusted display results, with lambda 1.122
for BD and 1.222 for MDD. This does not mean those lambdas should be applied
again before meta-analysis.

We checked this directly against published full-European sentinel statistics:

| Trait | Published sentinels compared | Best match |
|---|---:|---|
| BD | 258 | delivered raw beta and SE |
| MDD | 615 | delivered raw beta and SE |

For BD, the mean absolute beta error was `3.06e-5` with raw SE versus
`3.52e-4` after multiplying SE by `sqrt(1.122)`. For MDD it was `3.96e-5`
with raw SE versus `3.02e-4` after multiplying by `sqrt(1.222)`. Published SEs
are rounded, but the result is decisive: do not apply an extra lambda
correction. The script does not do so.

## What the script does

1. Validates input files, tools, headers, and the assertion of no sample
   overlap.
2. Streams the 57.6 million annotation and association rows together, verifies
   every internal ID, and applies the documented release-specific coordinate
   offset before normalization. BD v7.0 requires `-1`; MDD v7.2 requires `0`.
3. Keeps `pass=Y`, biallelic autosomal SNVs with valid beta, SE, P, case count,
   and control count.
4. Computes per-variant 23andMe effective N as
   `4 * N_case * N_control / (N_case + N_control)`.
5. Converts both components to normalized GRCh37 GWAS-VCF BCFs with effects
   relative to exact REF/ALT alleles.
6. Runs inverse-variance fixed-effect meta-analysis with heterogeneity and
   effect-direction output.
7. Applies the paper-specific variant coverage rule: at least 75% of maximum
   effective N for BD or 80% for MDD.
8. Lifts the final result once from GRCh37 to GRCh38.
9. Writes a dense GRCh38 BCF for colocalization, a canonical tabular registry,
   a PRSice base table, indexes, checksums, and a Markdown run report.

The final BCF `ES` and `SE` are the combined log-odds beta and SE relative to
`ALT`. `NS`, `NC`, and `NE` are summed across the two independent components.
The tabular `VARIANT` and PRSice `SNP` use `CHROM:POS:REF:ALT`.

`bcftools +metal` combines association statistics but does not emit a combined
imputation score. Public HRC INFO and 23andMe `avg.rsqr` are component-specific;
there is no justified arithmetic rule that turns them into one paper-equivalent
INFO value. The final BCF therefore defines `FORMAT/SI` but leaves it missing.
This keeps field queries compatible without inventing data.

## Commands

MDD has matching v7.2 association and annotation releases:

```bash
scripts/integrate_23andme_gwas.sh \
  --trait MDD \
  --out-prefix /home/gpertea/work/ref/GWAS/MDD/full_eur_integration/pgc-mdd2025_eur_v3-49-24-11 \
  --confirm-no-sample-overlap \
  --threads 8
```

BD uses its matching v7.0 European association and annotation releases. The
v7.0 annotation positions are one base greater than GRCh37 VCF positions for
the retained autosomal SNVs, so the BD profile applies an explicit `-1`
offset:

```bash
scripts/integrate_23andme_gwas.sh \
  --trait BD \
  --out-prefix /home/gpertea/work/ref/GWAS/BD/full_eur_integration/bip2024_eur \
  --confirm-no-sample-overlap \
  --threads 8
```

`--allow-bd-v7.2-annotations` is retained only to reproduce the earlier
provisional mapping run. A complete cross-release audit found identical IDs,
rsIDs, chromosomes, alleles, ploidy, and strand across all 57,611,376 rows.

`--resume` reuses existing nonempty stage outputs. It is intentionally not a
force-overwrite option. Each run preserves public, 23andMe, unfiltered-meta,
GRCh37-final, and GRCh38-final BCFs so a result can be audited.

## Outputs

For an output prefix `PREFIX`:

| File | Purpose |
|---|---|
| `PREFIX.23andme.grch37.ssf.tsv.gz` | prepared 23andMe component |
| `PREFIX.public-no23.grch37.bcf` | normalized public component |
| `PREFIX.23andme.grch37.bcf` | normalized 23andMe component |
| `PREFIX.meta-unfiltered.grch37.bcf` | meta-analysis before effective-N filter |
| `PREFIX.full.grch37.bcf` | filtered combined GRCh37 result |
| `PREFIX.full.grch38.bcf` | dense coloc and registry input |
| `PREFIX.full.grch38.meta.tsv.gz` | canonical combined summary statistics |
| `PREFIX.full.grch38.prsice.tsv.gz` | PRSice base table |
| `PREFIX.integration-report.md` | provenance, thresholds, counts, caveats |

## Completed production runs

MDD completed on 2026-07-26. BD was rebuilt from the recovered matching v7.0
annotation on 2026-07-30. Both runs processed all 57,611,376 paired
annotation/association rows with zero ID-alignment, statistic, sample-size,
source, or corrected-position errors.

| Trait | 23andMe autosomal SNVs prepared | Maximum NE | Applied minimum NE | Final GRCh38 variants | Variants at P <= 5e-8 |
|---|---:|---:|---:|---:|---:|
| BD | 21,137,709 | 440,999 | 330,749.25 | 6,394,788 | 10,753 |
| MDD | 21,098,300 | 1,577,200 | 1,261,760 | 6,656,222 | 36,787 |

Stable installed files are:

```text
/home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.bcf
/home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.meta.tsv.gz
/home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.prsice.tsv.gz
/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf
/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.meta.tsv.gz
/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.prsice.tsv.gz
```

All stable full-data names omit `_no23andMe`. Stage files, checksums, and exact
run-specific counts are retained under each trait's `full_eur_integration`
directory. These generated data files are not part of the Git repository.

The corrected v7.0 BD prepared component, canonical meta table, and PRSice
table are byte-for-byte identical to the earlier v7.2-mapped products. The
final BCF record stream is also identical; only bcftools run timestamps in the
BCF header differ. The installed BD BCF was retained to preserve the exact hash
already recorded by downstream analyses. The old run is preserved under
`BD/full_eur_integration_v7.2_provisional_2026-07-26/`.

The completed EUR LDGM `+pgs` outputs are:

| Trait | PGS BCF | Records | Logged seed |
|---|---|---:|---:|
| BD | `BD/bip2024_eur.hg38.pgs.b5e-8.bcf` | 5,227,123 | 1,785,102,265 |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf` | 5,339,102 | 1,785,102,085 |

The seeds were time-derived in the completed runs. Reproduction commands make
the logged values explicit.

Shizhong's legacy BD script takes OR and applies `log(OR)`. The canonical table
contains both `BETA` and `OR`; new common code should use `BETA` directly for
both BD and MDD. PRSice should be run with `--beta`.

## Checks and deliberate limits

- The script requires an explicit no-overlap assertion. If participants overlap
  across components, ordinary inverse-variance weights are wrong unless the
  covariance is modeled.
- It uses exact reference alleles and normalized variants, not rsID alone.
- Symbolic 23andMe `D/I` alleles are excluded because they do not identify an
  exact sequence allele. The current public BD and MDD BCFs contain only 2 and
  1 indels, respectively, so almost all shared variants are SNVs.
- BD paper-level QC also used DENTIST after meta-analysis. DENTIST is not
  installed locally. A paper-provided GRCh37 exclusion list can be applied with
  `--exclude-list`; without it, the result is a valid fixed-effect integration
  but not an exact reconstruction of every paper QC exclusion.
- The final beta, SE, and P can be reconstructed from aggregate components, but
  original cohort-level heterogeneity cannot. The emitted heterogeneity is
  between the public aggregate and 23andMe aggregate.
- The current spatialDLPFC `loadGWAS()` significance query requires `SI >= 0.8`
  and would drop missing SI. Before activating these BCFs, change that filter to
  retain missing SI explicitly, or implement a documented component-level
  quality rule. Dense coloc extraction already retains missing SI after query.
- Significant variant counts are not counts of independent signals or loci.

Before making either result active in PRS or coloc, compare hundreds of final
sentinel beta/SE/P values with the published supplement, inspect allele
direction, check maximum effective N, count threshold exclusions, and review
the integration report.

## Best-practice and study sources

- [METAL method, Willer et al. 2010](https://doi.org/10.1093/bioinformatics/btq340)
- [GWAS meta-analysis QC protocol, Winkler et al. 2014](https://doi.org/10.1038/nprot.2014.071)
- [GWAS Catalog summary-statistics harmonization](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics)
- [O'Connell et al. BD paper](https://doi.org/10.1038/s41586-024-08468-9)
- [Adams et al. MDD paper](https://doi.org/10.1016/j.cell.2024.12.002)
- [`bcftools +metal` implementation](https://github.com/freeseek/score)

The method is established best practice. The remaining BD uncertainty is not
the meta-analysis formula or annotation mapping; it is exact reproduction of
the paper's post-meta-analysis DENTIST QC.
