# European GWAS and PGS files on GRCh38

The current BD and MDD datasets combine the public European no-23andMe
meta-analysis with the independent European 23andMe component. SCZD remains
the public European dataset. GWAS files use GWAS-VCF BCF on GRCh38; PGS files
are EUR LDGM GraphPred loadings from `bcftools +pgs`.

Large BCF, index, compressed table, LDGM, log, and generated output files stay
under `~/work/ref/GWAS` and are not committed to Git.

## Current files

| Trait | GWAS BCF | Sample | GWAS records | PGS BCF | PGS sample | PGS records |
|---|---|---|---:|---|---|---:|
| BD | `BD/bip2024_eur.hg38.bcf` | `BD_2024_FULL_EUR` | 6,394,788 | `BD/bip2024_eur.hg38.pgs.b5e-8.bcf` | `BD_2024_FULL_EUR_pgs_a0.5_b5e-08` | 5,227,123 |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf` | `MDD_2025_FULL_EUR` | 6,656,222 | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf` | `MDD_2025_FULL_EUR_pgs_a0.5_b2e-08` | 5,339,102 |
| SCZD | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf` | `SCZ_2022.EUR` | 7,658,487 | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf` | `SCZ_2022.EUR_pgs_a0.5_b2e-07` | 6,076,466 |

The full BD and MDD names deliberately omit `_no23andMe`. Each BCF has a
matching `.bcf.csi` index. The full integrations also have stable tabular
exports:

```text
BD/bip2024_eur.hg38.meta.tsv.gz
BD/bip2024_eur.hg38.prsice.tsv.gz
MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.meta.tsv.gz
MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.prsice.tsv.gz
```

The canonical `meta.tsv.gz` table contains `CHROM`, `POS`, `RSID`, canonical
`CHROM:POS:REF:ALT`, alleles, beta, SE, P, OR, and sample-size fields. The
`prsice.tsv.gz` table contains `CHR BP SNP A1 A2 BETA P` and is ready for
PRSice with `--beta`. Shizhong's OR-based method can use `OR` from the
canonical table; new shared code should use `BETA` directly.

## Source data

| Trait | Component | Local source | Release |
|---|---|---|---|
| BD | Public EUR, no 23andMe | `BD/bip2024_eur_no23andMe.gz` | O'Connell et al., doi:10.1038/s41586-024-08468-9; Figshare doi:10.6084/m9.figshare.27216117 |
| BD | 23andMe-only EUR | `23andMe_MDD_BD/Bipolar-Disorder-O_Connell-2025/OConnell_2025_bipolar_european-7.0/bipolar.dat.gz` | association release 7.0 |
| MDD | Public EUR, no 23andMe | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz` | Adams et al., doi:10.1016/j.cell.2024.12.002; Figshare doi:10.6084/m9.figshare.27061255 |
| MDD | 23andMe-only EUR | `23andMe_MDD_BD/MDD-Adams-2025/Adams_2025_mdd_european-7.2/mdd.dat.gz` | association release 7.2 |
| BD and MDD | 23andMe EUR annotation | `23andMe_MDD_BD/7.2-Annotations/v7.2_europe/all_snp_info.txt.gz` | annotation release 7.2 |
| SCZD | Public EUR | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz` | Trubetskoy et al., doi:10.1038/s41586-022-04434-5; Figshare doi:10.6084/m9.figshare.19426775 |

The BD v7.0 association file was delivered without v7.0 European annotation.
The integration used the supplied v7.2 annotation only after explicit risk
acknowledgement. Published sentinel comparisons support the mapping, but BD
remains provisional until 23andMe supplies v7.0 annotation or confirms ID
stability. The BD paper's DENTIST exclusions were not available and therefore
were not replicated.

## BD and MDD integration

The components were not concatenated. For each exact allele-aligned variant,
the script combined log-odds estimates by standard-error inverse-variance
fixed-effect meta-analysis with `bcftools +metal`, applied the paper's
effective-N coverage threshold, and lifted the combined result once to GRCh38.

| Trait | Maximum NE | Minimum retained NE | GRCh38 records | Variants at P <= 5e-8 |
|---|---:|---:|---:|---:|
| BD | 440,999 | 330,749.25 (75%) | 6,394,788 | 10,753 |
| MDD | 1,577,200 | 1,261,760 (80%) | 6,656,222 | 36,787 |

Significant-variant counts are not counts of independent signals or loci.
Details, equations, validation, commands, checksums, and known limits are in
`mbv-prs/GWAS_23andMe_integration.md` and each
`full_eur_integration/*.integration-report.md` external run report.

## GWAS BCF fields

All current GWAS BCF effects are relative to `ALT`. Match genotype data by
GRCh38 `CHROM`, `POS`, `REF`, and `ALT`, not rsID alone.

The integrated BD and MDD BCFs contain these FORMAT fields:

| FORMAT field | Meaning |
|---|---|
| `NS` | variant-specific total sample count |
| `NC` | variant-specific case count |
| `ES` | log-odds beta relative to `ALT` |
| `SE` | standard error of `ES` |
| `LP` | `-log10(P)` |
| `NE` | effective sample size |
| `I2` | Cochran heterogeneity I squared |
| `CQ` | Cochran Q `-log10(P)` |
| `ED` | effect direction across components |
| `SI` | combined imputation quality; missing for integrated BD and MDD |

Public HRC INFO and 23andMe `avg.rsqr` are component-specific and cannot be
combined by a justified arithmetic rule. `SI` is therefore missing rather
than fabricated. Code requiring `SI >= 0.8` must explicitly retain missing
`SI` for these integrated files or use a separately documented component QC
rule.

Current rsID coverage is 6,360,316 of 6,394,788 BD variants and 6,645,867 of
6,656,222 MDD variants. Canonical coordinates remain the primary identifier.

The SCZD BCF retains `NS:SI:NC:ES:SE:LP:NE`. PGS BCFs are reduced outputs
containing only the GraphPred loading in `FORMAT/ES`. Liftover can add `FLIP`
and `SWAP` INFO annotations. `IFFY` and `REF_MISMATCH` are the defined site
filters.

The shared source-header mapping is
`mbv-prs/scripts/gwas_meta_colheaders.tsv`. It maps common SNP, coordinate,
allele, beta/OR, SE, P, INFO, sample-size, heterogeneity, direction, and
frequency headers to `bcftools +munge` fields. Case-only and control-only
frequencies are intentionally not treated as one overall effect-allele
frequency.

## PGS generation

The same explicit parameters used for the old no-23andMe loadings were retained
so that the GWAS input is the intended change:

| Trait | EUR LDGM | `--beta-cov` | `--max-alpha-hat2` | Seed | Exclusion |
|---|---|---:|---:|---:|---|
| BD | `1kg_ldgm.EUR.bcf` | `5e-8` | `0.001` | 1,785,102,265 | `FILTER="IFFY"` |
| MDD | `1kg_ldgm.EUR.bcf` | `2e-8` | `0.0005` | 1,785,102,085 | `FILTER="IFFY"` |
| SCZD | `1kg_ldgm.EUR.bcf` | `2e-7` | `0.002` | see archived log | `FILTER="IFFY"` |

Exact integration, installation, and `+pgs` commands are in
`tutorial_ldgm_gwas_hg38_pgs.md`.

## Archived no-23andMe files

These remain available for reproduction and direct before/after comparison,
but they are not the current BD/MDD defaults:

| Trait | GWAS BCF | Records | PGS BCF | PGS records |
|---|---|---:|---|---:|
| BD | `BD/bip2024_eur_no23andMe.hg38.bcf` | 6,938,764 | `BD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf` | 5,660,257 |
| MDD | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf` | 7,362,678 | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf` | 5,737,561 |

The full files can have fewer rows despite larger cohorts because the full
paper workflows impose high effective-N coverage and the integration retains
exactly mapped, QC-passing, shared autosomal SNVs.

## Inspection

```bash
bcftools view -h BD/bip2024_eur.hg38.bcf
bcftools query -l MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf
bcftools index -n BD/bip2024_eur.hg38.pgs.b5e-8.bcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES]\n' \
  MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf | head
```
