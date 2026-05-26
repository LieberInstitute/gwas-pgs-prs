# gwas-pgs-prs

Reusable tools for preparing GWAS summary statistics and calculating PRS/PGS from genotype data.

The preferred workflow is to standardize each GWAS into GRCh38 GWAS-VCF/BCF first, then use that same BCF for both scoring paths:

1. `bcftools +munge` converts source summary statistics to GWAS-VCF/BCF.
2. `bcftools +liftover` remaps the GWAS BCF to GRCh38.
3. `bcftools +pgs` estimates LD-aware GraphPred/LDGM PGS loadings.
4. `bcftools +score` applies those loadings to target genotype dosages.
5. PRSice-2 uses a base file extracted from the same standardized GWAS BCF for clumping and thresholding.
6. Comparison scripts summarize and correlate the resulting score sets.

## Inputs

Typical project inputs are:

- source GWAS summary statistics, usually TSV or TSV.GZ;
- a column map such as `colheaders.tsv`, or a `bcftools +munge` preset;
- source reference FASTA, target GRCh38 FASTA, and liftover chain when the GWAS is not already GRCh38;
- target sample genotypes with imputed dosages, preferably VCF/BCF with `FORMAT/DS`;
- LDGM files for `bcftools +pgs`;
- PLINK LD reference for PRSice-2;
- optional phenotype or demographic table for QC and score summaries.

The standardized GWAS BCF stores effect sizes in `FORMAT/ES` relative to `ALT` and p-values as `FORMAT/LP`, where `LP = -log10(P)`.

## Reference Files

The BPD, MDD, and SCZD GWAS BCFs documented in `README_GWAS.md` were converted from GRCh37/hg19 summary statistics with these local files:

```text
GRCh37 source FASTA: /dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta
GRCh38 liftover FASTA: /dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
UCSC liftover chain: /dbdata/cdb/ref/hg19ToHg38.over.chain.gz
```

These paths resolve the older tutorial's `~/work/ref/...` paths. The chain file is the standard UCSC `hg19ToHg38.over.chain.gz`.

Target genotype QC may use a different but compatible GRCh38 reference. The MBv workflow uses:

```text
/dbdata/cdb/genotyping/ref/hg38c.fa
```

Use the GWAS liftover FASTA for `bcftools +munge/+liftover`, and the target QC FASTA for normalizing the target genotype VCF.

## Scripts

Generic scripts are in `scripts/`.

| Script | Purpose |
|---|---|
| `prepare_gwas_vcf.sh` | Convert source GWAS summary statistics with `bcftools +munge`, optionally lift to GRCh38 with `bcftools +liftover`, then sort and index. |
| `compute_pgs_loadings.sh` | Run `bcftools +pgs` on standardized GWAS BCF files to create GraphPred/LDGM loading BCFs. |
| `prepare_target_qc.sh` | QC target genotypes, normalize variants, preserve a chr-prefixed BCF for bcftools scoring, and export BGEN for PRSice. |
| `score_bcftools.sh` | Run `bcftools +score --use DS --sample-header --counts` for sample-level dosage-weighted scores. |
| `build_prsice_base_from_gwas_vcf.sh` | Build PRSice base files from standardized GWAS BCF using `A1=ALT`, `A2=REF`, `BETA=ES`, and `P=10^-LP`. |
| `run_prsice_scores.sh` | Run PRSice-2 fixed-threshold C+T scoring with `--print-snp` enabled. |
| `export_prsice_model_tables.R` | Join PRSice retained SNPs back to the base file and export SNP-level model tables by threshold. |
| `summarize_prs_scores.R` | Merge arbitrary score files into long, raw-wide, and z-score-wide tables. |
| `compare_prs_methods.R` | Compare score methods by label, threshold, and pairwise correlation. |
| `run_prs_pipeline.sh` | Stage-based wrapper around the scripts above. |

The old MBv-specific scripts under `mbv-prs/scripts/` remain useful as a worked historical run, but new projects should use the generic top-level scripts and manifests.

## Standard Outputs

Recommended output layout:

```text
work/<project>/
  gwas/                 standardized GRCh38 GWAS BCFs
  pgs/                  bcftools +pgs loading BCFs
  target/qc/            QC target PGEN/BCF files
  target/prsice/        PRSice BGEN target files
  scores_bcftools/      sample scores from bcftools +score
  prsice/base/          PRSice base files
  prsice/out/           PRSice score outputs
  prsice_model_export/  retained SNP and threshold model tables
  summary/              merged score tables
  comparison/           method comparison tables
  tmp/                  temporary files
```

## bcftools Path

The bcftools path has two separate steps:

- `bcftools +pgs` computes LD-aware variant loadings from GWAS summary statistics and LDGM files.
- `bcftools +score` applies those loadings to target sample dosages.

The score for a sample is a weighted dosage sum over matched variants. In this repository the target dosage field is `FORMAT/DS`.

## PRSice-2 Path

PRSice-2 receives a plain base table built from the same standardized GWAS BCF:

```text
CHR BP SNP A1 A2 BETA P
```

For these BCFs:

- `A1` is `ALT`;
- `A2` is `REF`;
- `BETA` is `FORMAT/ES`;
- `P` is reconstructed as `10^-FORMAT/LP`.

`run_prsice_scores.sh` enables `--print-snp`, so the retained clumped SNPs can be exported as explicit model tables by `export_prsice_model_tables.R`.

## Method Comparison

Use `summarize_prs_scores.R` and `compare_prs_methods.R` to compare:

- GraphPred/LDGM scores from `bcftools +score`;
- PRSice C+T scores at fixed p-value thresholds;
- optional external C+T scores, including Shizhong-style outputs.

The comparison output is sample-level and method-level. It does not require publishing SNP-level weights.

## Example

The directory `examples/mbv_bpd/` contains a BPD plus MBv example config and manifests. See `README_tutorial.md` for the full workflow.

The example manifest now includes the concrete reference paths used by the earlier BPD/MDD/SCZD GWAS conversion.
