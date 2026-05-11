# Updated MBv PRS Calculation Plan

## Summary

Compute PRS for 119 MBv genotypes for SCZD, MDD, and BPD with two methods:

1. `bcftools +score` using existing LDGM/GraphPred PGS BCF loadings.
2. PRSice-2 using clumping and thresholding from converted GWAS summary statistics.

This updated plan includes corrections learned from the real run:

- Apply variant missingness before sample missingness because the input VCF is a union across imputation batches/arrays.
- Preserve `chr` contig names for `bcftools +score`; PLINK2 BCF export changes contigs to numeric names.
- Deduplicate PRSice base SNP IDs before running PRSice.
- Provide a dummy phenotype file for PRSice BGEN input while keeping `--no-regress`.
- Normalize PRSice duplicate-form sample IDs back to MBv brain IDs in merged outputs.

## Inputs And Outputs

Inputs:

- Target VCF: `genotypes/merged_R.8_MAF.01.RSann.vcf.gz`
- Demographics: `MBv_demographics_n119.tab`
- Reference FASTA: `/dbdata/cdb/genotyping/ref/hg38c.fa`
- EUR LD panel for PRSice: `/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos`
- GraphPred PGS BCFs:
  - BPD: `/home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.pgs.b5e-8.bcf`
  - MDD: `/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf`
  - SCZD: `/home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf`
- GWAS BCFs for PRSice base files:
  - BPD, MDD, SCZD `.hg38.bcf` files under `/home/gpertea/work/ref/GWAS`

Main outputs:

- `qc/target_qc_report.tsv`
- `scores_bcftools/mbv_graphpred_scores.tsv`
- `scores_bcftools/mbv_graphpred_scores_z.tsv`
- `prsice/out/mbv_prsice_scores.tsv`
- `prsice/out/mbv_prsice_scores_z.tsv`
- `mbv_prs_merged_scores.tsv`
- `README_MBV_PRS.md`
- `logs/commands_successful.tsv`

## Updated Target QC

1. Record versions and checksums.
2. Confirm VCF samples match demographics exactly.
3. Build normalized autosomal SNP BCF:
   - chr1-22 only
   - biallelic SNPs only
   - `FILTER=PASS` or `.`
   - `INFO/R2 >= 0.8`
   - `INFO/MAF >= 0.01`
   - REF check against `/dbdata/cdb/genotyping/ref/hg38c.fa` with error on mismatch
   - remove exact duplicate records
   - set IDs to `CHROM:POS:REF:ALT`
4. Audit dosage:
   - require usable `DS`
   - write `qc/ds_out_of_range_sites.tsv`
   - exclude any site with `DS < 0` or `DS > 2`
5. Convert to PLINK2 PGEN with dosage for QC.
6. Remove duplicate target variant IDs with `--rm-dup force-first list`.
7. Apply target QC in this order:
   - variant missingness: `--geno 0.02`
   - MAF: `--maf 0.01`
   - sample missingness: `--mind 0.02` after variant filtering
   - control-only HWE: `--hwe 1e-6 0 keep-fewhet`
8. Run LD pruning, PCA, and KING relatedness.
9. Export:
   - PLINK/PGEN QC set for QC bookkeeping
   - BGEN for PRSice
   - chr-prefixed BCF rebuilt from the normalized VCF for `bcftools +score`

Reason for the chr-prefixed BCF: PLINK2 export changes contig names to numeric chromosomes, which breaks exact matching to the GWAS-VCF/PGS files that use `chr` prefixes.

## Method 1: bcftools GraphPred Scores

Use `bcftools +score` with:

```bash
bcftools +score --use DS --sample-header --counts \
  -o scores_bcftools/<DISORDER>.graphpred.tsv \
  qc/mbv.qc.chr.bcf <DISORDER_PGS_BCF>
```

Requirements:

- Use `qc/mbv.qc.chr.bcf`, not PLINK2-exported `qc/mbv.qc.bcf`.
- Score only exact CHROM/POS/REF/ALT matches.
- Fail or flag if any disorder has zero matched markers.
- Merge disorder outputs and z-score columns with `scripts/summarize_mbv_prs.R`.

Expected successful overlap from the real run:

```text
BPD: 4996043
MDD: 4832784
SCZD: 5128978
```

## Method 2: PRSice-2 C+T Scores

Prepare PRSice base files from converted GWAS BCFs:

- columns: `CHR BP SNP A1 A2 BETA P`
- `SNP = CHROM:POS:REF:ALT`
- `A1 = ALT`
- `A2 = REF`
- `BETA = FORMAT/ES`
- `P = 10^(-FORMAT/LP)`
- remove ambiguous A/T and C/G SNPs by PRSice default
- deduplicate by first `SNP` occurrence before gzip

Run PRSice with:

```bash
PRSice_linux \
  --base prsice/base/<DISORDER>.base.tsv.gz \
  --target prsice/target/mbv_qc \
  --type bgen \
  --ld /dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos \
  --ld-type bed \
  --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 \
  --stat BETA --pvalue P --beta \
  --ignore-fid \
  --pheno prsice/target/dummy.pheno --pheno-col DUMMY \
  --no-regress \
  --all-score --score sum --fastscore \
  --bar-levels 5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1 \
  --clump-kb 250 --clump-r2 0.1 --clump-p 1
```

Dummy phenotype rule:

- Create `prsice/target/dummy.pheno` with columns `IID DUMMY`.
- Use IDs matching BGEN sample interpretation, e.g. `Br1092_Br1092`.
- Keep `--no-regress`; the dummy phenotype is only for PRSice input validation.

After PRSice:

- Merge `*.all_score`.
- Convert duplicate-form IDs `BrXXXX_BrXXXX` back to `BrXXXX`.
- Z-score PRSice score columns.

Expected successful clumped variant counts from the real run:

```text
BPD: 224365
MDD: 231809
SCZD: 237872
```

## Checks And Acceptance Criteria

Required checks:

- `comm -3 qc/demog.samples qc/vcf.samples` is empty.
- `qc/ds_out_of_range_sites.tsv` is empty or excluded sites are documented.
- `qc/target_qc_report.tsv` reports 119 input and 119 final samples.
- KING relatedness table has no pairs at or above 0.0884 unless explicitly handled.
- GraphPred matched marker counts are nonzero and reported.
- PRSice logs report nonzero clumped variants for each disorder.
- Final output row counts are 119 samples plus header:
  - `scores_bcftools/mbv_graphpred_scores.tsv`
  - `prsice/out/mbv_prsice_scores.tsv`
  - `mbv_prs_merged_scores.tsv`
- `bash -n scripts/run_mbv_prs.sh` passes.
- `Rscript -e 'parse("scripts/summarize_mbv_prs.R")'` passes.
- No `bcftools`, `plink2`, or `PRSice` jobs remain running at the end.

## Production Considerations

Add explicit production guards for:

- missingness filter order on union-imputed data
- duplicate variant IDs in both target and base files
- contig naming drift after PLINK export
- allele orientation and `CHROM:POS:REF:ALT` consistency
- ambiguous SNP exclusion counts
- zero-marker score output
- sample ID rewriting by BGEN/PRSice
- external LD panel ancestry, build, and ID convention
- stale files from failed runs, such as `.valid` files
- score table row-count mismatches after merging methods

Do not treat PRSice and GraphPred scores as interchangeable. GraphPred scores apply LDGM-adjusted loadings; PRSice scores are clumped-and-thresholded scores from GWAS marginal effects.
