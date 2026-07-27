# Shizhong clump.pl Method And Available Score Files

Everything needed for score comparison is available.

## Files Found

Shizhong scores in the current directory:

- `PRS_Bipolar.csv`
- `PRS_MDD.csv`
- `PRS_SCZ.csv`

Each has 119 samples plus a header and 18 p-value threshold score columns. All sample IDs match our score files exactly.

Our bcftools GraphPred scores:

- `../scores_bcftools/mbv_graphpred_scores.tsv`
- `../scores_bcftools/mbv_graphpred_scores_z.tsv`

These have 119 samples and 3 score columns: BPD, MDD, and SCZD.

Our PRSice scores:

- `../prsice/out/mbv_prsice_scores.tsv`
- `../prsice/out/mbv_prsice_scores_z.tsv`

These have 119 samples and 27 PRSice threshold columns. They are under `../prsice/out/`, not `../scores_bcftools/`.

Existing merged output:

- `../mbv_prs_merged_scores.tsv`

This contains our GraphPred plus our PRSice scores, but not Shizhong scores.

## File Caveat

Use the top-level Shizhong CSVs for comparison:

- `PRS_Bipolar.csv`
- `PRS_MDD.csv`
- `PRS_SCZ.csv`

`PRS_MDD.csv` matches `shan-prs-keri/mdd/out/PRS_scores_by_threshold_mdd.csv`.

`PRS_Bipolar.csv` and `PRS_SCZ.csv` match the current `score.*.sscore` files exactly.

Older-looking files under `shan-prs-keri/bp/out/PRS_Bipolar.csv` and `shan-prs-keri/scz/out/PRS_SCZ.csv` differ slightly from the top-level CSVs and current `.sscore` files.

## Was clump.pl Used?

Yes. `shan-prs-keri/{bp,mdd,scz}/clump.pl` is the clumping script. Each trait has a `clump_jobs.txt` file that launches:

```bash
perl clump.pl --chr 1
perl clump.pl --chr 2
...
perl clump.pl --chr 22
```

The per-chromosome logs confirm PLINK 1.9 `--clump` was run.

## clump.pl Method

For each chromosome, `clump.pl` intersects GWAS variants with:

- 1KG EUR no-indel hg38 MAF 0.01 reference BIM files;
- target `merged_maf01_snps.pvar`;
- the trait-specific GWAS weight file.

Variants are matched by chromosome, position, and unordered allele pair.

The script then runs PLINK 1.9 clumping:

```bash
plink \
  --bfile EUR_noIndels_chr<chr>_maf0.01_hg38 \
  --clump chr<chr> \
  --clump-p1 1 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --out chr<chr>
```

Then it builds:

- `out/profile` with columns `SNP A1 BETA`;
- `out/profile_p` with columns `SNP P`.

Final scoring is done with PLINK2:

```bash
plink2 \
  --pfile /dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps \
  --score ./out/profile header ignore-dup-ids \
  --q-score-range prange ./out/profile_p header \
  --threads 1 \
  --out ./out/score
```

The `prange` thresholds are:

```text
1e-08
1e-07
1e-06
1e-05
1e-04
0.001
0.01
0.05
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
```

Final CSV values are PLINK2 `SCORE1_AVG` values per threshold.

## Processed Variant Counts

PLINK2 scoring logs report:

- BPD: 249285 variants processed.
- MDD: 257109 variants processed.
- SCZ: 140529 variants processed.

## Relationship To Our Methods

Shizhong scores are LD-aware clumping-and-thresholding scores produced by custom PLINK 1.9 clumping plus PLINK2 scoring. Here, LD is used for variant selection: PLINK clumping keeps index variants and removes nearby correlated variants above the configured `r2` cutoff.

Our PRSice scores are also LD-aware clumping-and-thresholding scores, but they use PRSice-2, 250 kb clump windows, r2 0.1, 9 thresholds, and validated canonical GWAS base files.

Our bcftools GraphPred scores are different in how they use LD: LDGM/GraphPred uses LD structure upstream to estimate adjusted genome-wide loadings, then `bcftools +score` applies those loadings directly to target dosages. It produces one score per disorder, not a p-value threshold grid.
