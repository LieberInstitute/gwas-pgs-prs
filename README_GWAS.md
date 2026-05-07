# GWAS Summary Statistics on GRCh38

This directory contains psychiatric GWAS summary statistics converted to GWAS-VCF BCF on GRCh38 and matched PGS loading outputs from the BCFtools `+pgs` plugin.

The source summary statistics were normalized with `bcftools +munge`, lifted from GRCh37/hg19 to GRCh38 with `bcftools +liftover`, sorted, and indexed. PGS loadings were computed from the converted GWAS BCFs with the EUR LDGM file and `bcftools +pgs`.

All main outputs are BCF files with matching `.bcf.csi` indexes. They use GRCh38 contig names with `chr` prefixes.

## Main Files

| Disorder | Converted GWAS BCF | Sample | Records | PGS BCF | PGS sample | PGS records |
|---|---|---:|---:|---|---:|---:|
| BPD | `BPD/bip2024_eur_no23andMe.hg38.bcf` | `BIP_2024.EUR` | 6938764 | `BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf` | `BIP_2024.EUR_pgs_a0.5_b5e-08` | 5660257 |
| MDD | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf` | `MDD_2025` | 7362678 | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf` | `MDD_2025_pgs_a0.5_b2e-08` | 5737561 |
| SCZD | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf` | `SCZ_2022.EUR` | 7658487 | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf` | `SCZ_2022.EUR_pgs_a0.5_b2e-07` | 6076466 |

The converted GWAS BCFs preserve the association summary statistics as GWAS-VCF FORMAT fields. The PGS BCFs are reduced outputs from `+pgs` and contain only the PGS loading/effect score in FORMAT/ES.

## VCF Columns and IDs

Each BCF uses standard VCF columns:

| Column | Meaning |
|---|---|
| `CHROM` | GRCh38 chromosome or contig, usually `chr1` to `chr22` for these autosomal outputs |
| `POS` | 1-based GRCh38 position |
| `ID` | source variant identifier from the GWAS input |
| `REF` | GRCh38 reference allele after liftover |
| `ALT` | alternate/effect allele used by GWAS-VCF FORMAT fields |
| `QUAL` | unset in these files |
| `FILTER` | site filter status |
| `INFO` | liftover metadata, when applicable |
| `FORMAT` | per-sample GWAS or PGS fields |
| sample column | one synthetic sample holding summary statistics |

If a dbSNP rsID is present, it is in VCF column 3, `ID`. Not every source identifier is an rsID: MDD and SCZD include some `chr:pos_ref_alt` style IDs. Current ID coverage:

| File type | BPD rsIDs | MDD rsIDs | SCZD rsIDs |
|---|---:|---:|---:|
| Converted GWAS BCF | 6938764 / 6938764 | 7362102 / 7362678 | 7637489 / 7658487 |
| PGS BCF | 5660257 / 5660257 | 5737542 / 5737561 | 6075044 / 6076466 |

For harmonization with genotype VCF/BCF files, match on GRCh38 `CHROM`, `POS`, `REF`, and `ALT`. Treat `ES` as relative to `ALT`.

## Converted GWAS BCF FORMAT Fields

The converted `.hg38.bcf` files store GWAS summary statistics in FORMAT fields. Field values are per alternate allele (`Number=A`).

| FORMAT field | Meaning |
|---|---|
| `NS` | variant-specific number of samples or individuals with called genotypes |
| `SI` | imputation accuracy score |
| `NC` | variant-specific number of cases |
| `ES` | effect size estimate relative to `ALT`; BETA sources stay beta, OR sources are converted by `+munge` to log effect |
| `SE` | standard error of `ES` |
| `LP` | `-log10(P)` for the effect estimate |
| `NE` | variant-specific effective sample size |
| `AF` | alternate allele frequency in the trait subset; present in BPD only |
| `I2` | Cochran heterogeneity I squared; present in BPD and MDD |
| `CQ` | Cochran Q `-log10(P)`; present in BPD only |
| `ED` | effect direction across studies; present in BPD only |

Per-disorder FORMAT layouts:

| Disorder | FORMAT fields in converted GWAS BCF |
|---|---|
| BPD | `NS:SI:NC:ES:SE:LP:AF:NE:I2:CQ:ED` |
| MDD | `NS:SI:NC:ES:SE:LP:NE:I2` |
| SCZD | `NS:SI:NC:ES:SE:LP:NE` |

## PGS BCF FORMAT Fields

The `.hg38.pgs.*.bcf` files are the `bcftools +pgs` outputs. They contain:

| FORMAT field | Meaning |
|---|---|
| `ES` | GraphPred PGS loading/effect score relative to `ALT` |

The PGS outputs were created with EUR LDGM and these options:

| Disorder | `--beta-cov` | `--max-alpha-hat2` | Input exclusion |
|---|---:|---:|---|
| BPD | `5e-8` | `0.001` | `FILTER="IFFY"` |
| MDD | `2e-8` | `0.0005` | `FILTER="IFFY"` |
| SCZD | `2e-7` | `0.002` | `FILTER="IFFY"` |

## INFO and FILTER Fields

The converted GWAS and PGS BCFs share liftover INFO definitions:

| INFO field | Meaning |
|---|---|
| `FLIP` | allele was strand-flipped during liftover |
| `SWAP` | alternate allele became reference during liftover; `-1` means a new reference allele was added |

Most records have `INFO=.`. Records with liftover changes may have `FLIP`, `SWAP=1`, `SWAP=-1`, or both.

The shared FILTER definitions are:

| FILTER | Meaning |
|---|---|
| `IFFY` | reference allele could not be determined |
| `REF_MISMATCH` | reference does not match any allele |

The inspected outputs currently have all records unfiltered as `FILTER=.`.

## colheaders.tsv

`colheaders.tsv` is the shared two-column header map used by `bcftools +munge` to translate source summary-statistic headers into canonical GWAS-VCF inputs. It is not disorder-specific.

Key mappings:

| Source headers | Canonical field |
|---|---|
| `SNP`, `ID` | `SNP` |
| `CHR`, `CHROM`, `#CHROM` | `CHR` |
| `BP`, `POS` | `BP` |
| `A1`, `EA` | `A1` |
| `A2`, `NEA` | `A2` |
| `OR` | `OR` |
| `BETA` | `BETA` |
| `SE` | `SE` |
| `P`, `PVAL` | `P` |
| `INFO`, `IMPINFO` | `INFO` |
| `NCAS`, `Nca` | `N_CAS` |
| `NCON`, `Nco` | `N_CON` |
| `NEFF` | `NEFF` |
| `NEFFDIV2`, `Neff_half` | `NEFFDIV2` |
| `HETI`, `HetISqt` | `HET_I2` |
| `HETPVAL`, `HetPVa` | `HET_P` |
| `Direction` | `DIRE` |
| `HRC_FRQ_A1` | `FRQ` |

`FCAS` and `FCON` are intentionally not mapped because they are case and control allele frequencies, not one overall effect allele frequency.

## Inspection Commands

Useful checks:

```bash
bcftools view -h BPD/bip2024_eur_no23andMe.hg38.bcf
bcftools query -l MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf
bcftools index -n SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES]\n' BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf | head
```
