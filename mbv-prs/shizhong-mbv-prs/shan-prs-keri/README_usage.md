# Shizhong C+T PRS Workflow Usage

This guide documents the clumping-and-thresholding PRS workflow in this directory tree, using `mdd/` as the worked example. The same structure is used for `bp/` and `scz/`, with disorder-specific GWAS weight files and column indices inside each `clump.pl`.

The workflow is documentation-only here. The scripts are not made portable by this README. Several paths are hardcoded to the original JHPCE environment.

## Directory Layout

```text
shan-prs-keri/
  data/
    merged_maf01.*
    merged_maf01_snps.*
    readme
  mdd/
    clump.pl
    00_clump_jobs.pl
    01_clump_combine.pl
    02_score_jobs.pl
    clump_jobs.pbs
    prange
    out/
```

`bp/` and `scz/` have the same workflow structure.

## What The Workflow Does

For each disorder, the workflow:

1. Reads target variants from a PLINK2 PGEN dataset.
2. Reads the 1000 Genomes EUR no-indel hg38 LD reference for each chromosome.
3. Reads a disorder-specific GWAS weight table.
4. Keeps variants present in all three sources, matching by chromosome, position, and unordered allele pair.
5. Runs PLINK 1.9 LD clumping per chromosome.
6. Builds PLINK2 score files:
   - `profile`: SNP, effect allele, beta
   - `profile_p`: SNP, p-value
7. Uses PLINK2 `--score` plus `--q-score-range` to create one PRS per p-value threshold.
8. Optionally runs R plotting scripts to make wide CSV summaries and PDFs.

For MDD, the final wide score table is:

```text
mdd/out/PRS_scores_by_threshold_mdd.csv
```

This is an LD-aware C+T workflow. "LD-aware" here means LD is used for clumping/variant selection before scoring. This differs from LDGM/GraphPred-style scoring, where LD is used upstream to estimate adjusted genome-wide weights before those weights are applied to target dosages.

## Role Of `data/`

The `data/` folder contains the target genotype data used for scoring.

`data/readme` says the data were copied from:

```text
/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/subsets/Nikos_Keri_MBv/plink2/
```

Current files:

```text
merged_maf01.pgen
merged_maf01.pvar
merged_maf01.psam
merged_maf01.log
merged_maf01_snps.pgen
merged_maf01_snps.pvar
merged_maf01_snps.psam
merged_maf01_snps.log
```

`merged_maf01.*` is the MAF-filtered target PGEN dataset:

- 119 samples
- 7707704 variants
- created from PLINK bed/bim/fam with `--maf 0.01`

`merged_maf01_snps.*` is the SNP-only target PGEN used for scoring:

- 119 samples
- 7313747 SNP variants
- no duplicate variant IDs were observed in `merged_maf01_snps.pvar`

Important caveat: `mdd/clump.pl` and `mdd/score_jobs.txt` do not use a relative `./data/merged_maf01_snps` path. They hardcode the original path:

```text
/dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps
```

In this copied project location, those original absolute paths are not available. To rerun the workflow here, update the hardcoded paths in the scripts or restore the original directory layout.

## Required Inputs

For the MDD workflow, the scripts expect:

1. Target genotype PGEN trio:

```text
/dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps.pgen
/dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps.pvar
/dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps.psam
```

2. 1000 Genomes EUR no-indel hg38 LD reference, one PLINK bed/bim/fam prefix per chromosome:

```text
/dcs04/lieber/statsgen/shizhong/database/1KG/QC/EUR_noIndels/hg38/EUR_noIndels_chr<chr>_maf0.01_hg38
```

3. MDD GWAS weight table:

```text
/dcs04/lieber/statsgen/shizhong/database/GWAS/PGC/mdd2025/gwas_hg38/pgc-mdd2025_no23andMe_eur_hg38
```

4. PLINK 1.9 binary for clumping:

```text
/dcl02/lieber/shan/shizhong/software/plink/plink
```

5. PLINK2 binary for scoring:

```text
plink2
```

6. Optional phenotype table for plotting only:

```text
/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/genotyped_brnum_pheno_3125.tab
```

No phenotype is required for PRS scoring itself.

## Target Genotype QC Expectations

Before running this PRS workflow, prepare the target genotype data so that:

- genome build matches GWAS and LD reference, here hg38;
- variants are autosomal, because `clump.pl` loops over chromosomes 1 through 22;
- variants are biallelic SNPs; indels are not expected by the no-indel LD reference;
- target MAF is at least 0.01;
- variant IDs are stable and usable by PLINK2 scoring, here `chr:pos:REF:ALT`;
- target `.pvar` has no duplicate variant IDs;
- `.psam` sample IDs are final sample IDs for PRS output;
- REF/ALT and effect allele conventions are checked against the GWAS weights;
- target, GWAS, and LD reference use compatible chromosome naming and positions;
- standard sample and variant QC has already been performed.

Recommended pre-PRS QC includes sample missingness, variant missingness, allele frequency checks, HWE checks in controls when appropriate, relatedness checks, ancestry/PC checks, and confirmation that the final sample list is expected.

## MDD `clump.pl` Details

`mdd/clump.pl` uses these GWAS column indices:

```perl
$pvalc=7; ## p-value column
$chrc=0;  ## chromosome column
$posc=1;  ## hg38 position column
$weic=5;  ## beta/effect-size column
$a1c=3;   ## effect allele
$a2c=4;   ## other allele
```

For MDD, beta is used directly:

```perl
$weight{$bar} = $tokens[$weic]; ## beta
```

This differs from the BP and SCZ scripts, where the corresponding line uses `log(...)` around the weight value. Treat the column-index and transformation choices as trait-specific.

The clumping parameters are:

```perl
$clump_r2=0.1;
$clump_dist=1000;
```

The resulting PLINK 1.9 clumping command is:

```bash
/dcl02/lieber/shan/shizhong/software/plink/plink \
  --bfile /dcs04/lieber/statsgen/shizhong/database/1KG/QC/EUR_noIndels/hg38/EUR_noIndels_chr${chr}_maf0.01_hg38 \
  --clump chr${chr} \
  --clump-p1 1 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --out chr${chr}
```

`--clump-r2 0.1` is an LD pruning threshold, not a PRS score. Variants above this LD threshold with a selected index variant within the 1000 kb window are not kept as independent clumped variants.

## MDD Step-By-Step Commands

Run these commands from the `mdd/` directory.

```bash
cd mdd
```

Generate the per-chromosome clumping command file:

```bash
perl 00_clump_jobs.pl
```

This writes `clump_jobs.txt` with one command per chromosome:

```bash
perl clump.pl --chr 1
perl clump.pl --chr 2
...
perl clump.pl --chr 22
```

Run clumping on the cluster:

```bash
sbatch clump_jobs.pbs
```

Or run locally or interactively:

```bash
for chr in {1..22}; do
  perl clump.pl --chr "$chr"
done
```

Combine per-chromosome score profiles:

```bash
perl 01_clump_combine.pl
```

This creates:

```text
mdd/out/profile
mdd/out/profile_p
```

Generate the PLINK2 scoring command:

```bash
perl 02_score_jobs.pl
```

This writes `score_jobs.txt`. The generated command is:

```bash
plink2 \
  --pfile /dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps \
  --score ./out/profile header ignore-dup-ids \
  --q-score-range prange ./out/profile_p header \
  --threads 1 \
  --out ./out/score
```

Run scoring:

```bash
bash score_jobs.txt
```

Optional: create the wide score CSV and plots:

```bash
Rscript 04_boxplot_allGroups.R
Rscript 04_boxplot_mdd_control.R
Rscript 04_boxplot_bipolar_control.R
```

These R scripts require `data.table`, `stringr`, `ggplot2`, and `ggforce`, and they hardcode the phenotype table path listed above.

## P-Value Thresholds

The threshold file is `mdd/prange`:

```text
1e-08 0 1e-08
1e-07 0 1e-07
1e-06 0 1e-06
1e-05 0 1e-05
1e-04 0 1e-04
0.001 0 0.001
0.01 0 0.01
0.05 0 0.05
0.1 0 0.1
0.2 0 0.2
0.3 0 0.3
0.4 0 0.4
0.5 0 0.5
0.6 0 0.6
0.7 0 0.7
0.8 0 0.8
0.9 0 0.9
1 0 1
```

PLINK2 writes one `score.<threshold>.sscore` file per range.

## Outputs

Per chromosome:

```text
mdd/<chr>/chr<chr>
mdd/<chr>/chr<chr>.clumped
mdd/<chr>/chr<chr>.log
mdd/<chr>/chr<chr>.nosex
mdd/<chr>/profile
mdd/<chr>/profile_p
```

Combined clumped score inputs:

```text
mdd/out/profile
mdd/out/profile_p
```

Scoring outputs:

```text
mdd/out/score.<threshold>.sscore
mdd/out/score.log
```

Summary and plotting outputs:

```text
mdd/out/PRS_scores_by_threshold_mdd.csv
mdd/out/PRS_MDD.pdf
mdd/out/mdd_control_mddPRS_boxplot.pdf
mdd/out/bipolar_control_mddPRS_boxplot.pdf
```

The `.sscore` files contain:

```text
#IID
ALLELE_CT
NAMED_ALLELE_DOSAGE_SUM
SCORE1_AVG
```

The wide CSV has one row per sample and one score column per p-value threshold.

## Observed MDD Run Facts

From the existing logs and outputs:

- target scoring set: 119 samples and 7313747 SNP variants;
- `mdd/out/profile`: header plus 257109 variants;
- `mdd/out/profile_p`: header plus 257109 variants;
- PLINK2 processed 257109 variants;
- 18 threshold-specific `.sscore` files were produced;
- each `.sscore` file has 119 samples plus header;
- `mdd/out/PRS_scores_by_threshold_mdd.csv` has 119 samples plus header.

Example from `mdd/out/score.log`:

```text
--q-score-range: 18 ranges and 257109 variants loaded.
--score: 257109 variants processed.
```

## Basic Checks After Running

Check that all chromosome profile files exist:

```bash
for chr in {1..22}; do
  test -s "$chr/profile" || echo "missing $chr/profile"
  test -s "$chr/profile_p" || echo "missing $chr/profile_p"
done
```

Check combined profile sizes:

```bash
wc -l out/profile out/profile_p
```

Check scoring output count:

```bash
find out -maxdepth 1 -name 'score.*.sscore' | wc -l
```

Check sample counts:

```bash
for f in out/score.*.sscore; do
  printf "%s\t" "$f"
  wc -l < "$f"
done
```

Expected `.sscore` row count here is 120: 119 samples plus header.

## Porting Checklist

Before rerunning this workflow in a new location:

1. Replace hardcoded `/dcs05/.../data/merged_maf01_snps` with the intended target PGEN prefix.
2. Replace hardcoded 1KG EUR LD reference paths if needed.
3. Replace the PLINK 1.9 binary path if needed.
4. Confirm `plink2` is on `PATH`.
5. Confirm GWAS weight file path and column indices in `clump.pl`.
6. Confirm the GWAS beta transformation is correct for the trait.
7. Recreate `prange` if different thresholds are needed.
8. Update phenotype paths in R plotting scripts if plots are needed.

## Limitations

- Paths are hardcoded and environment-specific.
- The scripts assume hg38 positions and chromosome labels compatible with numeric chromosomes in the target PGEN and LD reference.
- `clump.pl` uses trait-specific GWAS column indices; copying it to a new GWAS requires checking every index.
- MDD uses beta directly, while BP and SCZ use `log(...)` around the weight value.
- The plotting scripts run t-tests and display p-values; those should be treated as exploratory in this small cohort.
- Choosing the best p-value threshold from the same target samples is exploratory unless validated externally or corrected with an appropriate resampling/permutation strategy.

## Relationship To The Top-Level Comparison

The top-level project compared Shizhong C+T scores with our PRSice C+T scores and our bcftools GraphPred scores. This README documents how the Shizhong MDD score source was produced.

Use `../PRS_MDD.csv` or `mdd/out/PRS_scores_by_threshold_mdd.csv` as the MDD Shizhong score table for downstream comparison.
