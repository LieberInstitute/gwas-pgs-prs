# GWAS gene and variant list provenance

This file documents the gene and variant list TSVs added beside the BPD, MDD, and SCZD GWAS summary statistics. These are paper/supplement-derived lists, not raw scans of every public summary-statistic row with P <= 5e-8.

The derived tables named below are kept in the external generated-output
archive documented in `README_generated_outputs.md`. Their repository-relative
names are retained below for provenance.

## Shared extraction rules

- Gene files use the common schema in `*_gene_lists.tsv`.
- Variant files use the common schema in `*_variant_lists.tsv` and related split files.
- Source coordinates were preserved as published. For BPD and MDD, source variant coordinates are hg19 where present. For SCZD, prioritized FINEMAP variant coordinates are GRCh37/hg19.
- `hg38_*` fields were filled by matching source rsIDs to the local GRCh38 BCFs with `bcftools query`; allele matching was used when source alleles were available.
- `primary_list=yes` marks the broadest or most central paper-level list for that disorder, not necessarily the smallest list.

## BPD 2024

Source paper: OConnell et al., Nature, doi:10.1038/s41586-024-08468-9.

Source workbook: `41586_2024_8468_MOESM4_ESM.xlsx`
URL: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-08468-9/MediaObjects/41586_2024_8468_MOESM4_ESM.xlsx`

Files:
- `BPD/bpd2024_gene_lists.tsv`: 116 rows from Table S31, "Prioritization of credible genes for multi-ancestry meta-analysis including self-report data".
- `BPD/bpd2024_prioritized_credible_genes.tsv`: the same 116 Table S31 rows as a separate focused file.
- `BPD/bpd2024_variant_lists.tsv`: 934 published locus/signal rows from Tables S5, S6, S8, S11, and S12.
- `BPD/bpd2024_finemapped_credible_variants.tsv`: 295 fine-mapped credible variant rows from Tables S28 and S29.

23andMe/self-report caveat:
- The local public BPD summary statistics BCF is `bip2024_eur_no23andMe.hg38.bcf`, so it excludes 23andMe/self-report data.
- The paper supplement is not limited to that public no23andMe file. It provides both including-self-report/23andMe and excluding-self-report/23andMe lists.
- The most extensive BPD paper lists are the including-self-report/23andMe lists. These are marked primary where applicable, and the no23andMe/EUR-specific lists are retained separately.
- Because some source fine-mapping rows have rsIDs but no source base-pair coordinate or allele columns, their `hg38_match_status` is usually `id_found_allele_not_checked`.

## MDD 2025

Source paper: Adams et al., Cell, doi:10.1016/j.cell.2024.12.002.

Source GWAS summary Figshare: `https://doi.org/10.6084/m9.figshare.27061255`

Source downstream results Figshare: `https://doi.org/10.6084/m9.figshare.27089614`

Files:
- `MDD/mdd2025_gene_lists.tsv`: 9896 rows combining fastBAT, Hi-C, and DrugTargetor/MAGMA gene lists.
- `MDD/mdd2025_fastBAT_bonferroni_significant_genes.tsv`: 1568 genes from `Online Results (fastBAT).xlsx`, sheet `fastBAT Results`, where `Bonf signif = YES`.
- `MDD/mdd2025_hic_significant_all_tissues_genes.tsv`: 1034 genes from `Online Results (hiC).xlsx`, sheet `HiC Gene Associations`, where `Significantinalltissues = YES`.
- `MDD/mdd2025_drugtargetor_magma_qBH_le_0_05_genes.tsv`: 7294 genes from `Online Results (DrugTargetor).xlsx`, sheet `G GENE_results`, where `q_valueBH <= 0.05`. This is marked primary because it is the largest significant gene set.
- `MDD/mdd2025_variant_lists.tsv`: 31588 rows combining COJO independent signals and fine-mapped credible causal variants.
- `MDD/mdd2025_cojo_independent_signals.tsv`: 1319 rows from `Online Results (COJO).xlsx`, sheets `COJO Multi-ancestry` and `COJO European ancestry`.
- `MDD/mdd2025_finemap_credible_causal_variants.tsv`: 30269 rows from the Figshare fine-mapping credible-causal files for diverse and EUR analyses.

23andMe caveat:
- The local public MDD BCF is based on `pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz`, so it excludes 23andMe.
- The gene and variant lists here were extracted from the paper's downstream Figshare results, not derived from the public no23andMe BCF.
- These MDD downstream files are paper-level results and are not sheet-labeled as no23andMe public-summary-statistic subsets. Treat them as less crippled by the public-summary restriction than the local no23andMe BCF, but preserve their exact source labels (`multi-ancestry`, `European ancestry`, `full_div`, `full_eur`) rather than relabeling them as with23andMe.

## SCZD 2022

Source paper: Trubetskoy et al., Nature, doi:10.1038/s41586-022-04434-5.

Source supplement zip: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-04434-5/MediaObjects/41586_2022_4434_MOESM11_ESM.zip`

Files:
- `SCZD/sczd2022_gene_lists.tsv`: 1072 rows combining Table S12 all-criteria genes, Table S12 focused prioritized genes, and SMR/Hi-C support tables.
- `SCZD/sczd2022_prioritized_genes.tsv`: 120 rows from `Supplementary Table 12.xlsx`, sheet `Prioritised`. This is the smaller focused prioritized gene list.
- `SCZD/sczd2022_prioritized_variants.tsv`: 1678 rows from `Supplementary Table 13.xlsx`, sheet `ST13 Prioritised NS+UTR pp>=.1`, and `Supplementary Table 14.xlsx`, sheet `ST14 Prioritised Single Gene`.

SCZD prioritized-list clarification:
- Yes, SCZD has a smaller focused prioritized gene list: Table S12, sheet `Prioritised`, 120 rows.
- SCZD also has smaller prioritized FINEMAP variant tables: Table S13 for nonsynonymous/UTR variants with posterior probability >= 0.1, and Table S14 for credible sets spanning a single gene.

## Verification summary

Current extracted row counts:

| File | Rows |
|---|---:|
| `BPD/bpd2024_gene_lists.tsv` | 116 |
| `BPD/bpd2024_prioritized_credible_genes.tsv` | 116 |
| `BPD/bpd2024_variant_lists.tsv` | 934 |
| `BPD/bpd2024_finemapped_credible_variants.tsv` | 295 |
| `MDD/mdd2025_gene_lists.tsv` | 9896 |
| `MDD/mdd2025_fastBAT_bonferroni_significant_genes.tsv` | 1568 |
| `MDD/mdd2025_hic_significant_all_tissues_genes.tsv` | 1034 |
| `MDD/mdd2025_drugtargetor_magma_qBH_le_0_05_genes.tsv` | 7294 |
| `MDD/mdd2025_variant_lists.tsv` | 31588 |
| `MDD/mdd2025_cojo_independent_signals.tsv` | 1319 |
| `MDD/mdd2025_finemap_credible_causal_variants.tsv` | 30269 |
| `SCZD/sczd2022_gene_lists.tsv` | 1072 |
| `SCZD/sczd2022_prioritized_genes.tsv` | 120 |
| `SCZD/sczd2022_prioritized_variants.tsv` | 1678 |

GRCh38 matching status for split variant files:

| File | Match summary |
|---|---|
| `BPD/bpd2024_finemapped_credible_variants.tsv` | 289 rsIDs found, alleles not checked; 6 rsIDs not found in local BCF |
| `MDD/mdd2025_cojo_independent_signals.tsv` | 1311 rsID+allele matches; 8 rsIDs not found in local BCF |
| `MDD/mdd2025_finemap_credible_causal_variants.tsv` | 30128 rsID+allele matches; 141 rsIDs not found in local BCF |
| `SCZD/sczd2022_prioritized_variants.tsv` | 1623 rsID+allele matches; 54 rsIDs not found in local BCF; 1 row without rsID |
