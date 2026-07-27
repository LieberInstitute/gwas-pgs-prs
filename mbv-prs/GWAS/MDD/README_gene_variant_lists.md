# MDD 2025 gene and variant lists

Source paper: Adams et al. Cell 2025 doi:10.1016/j.cell.2024.12.002.
Source GWAS summary Figshare: `https://doi.org/10.6084/m9.figshare.27061255`.
Source downstream results Figshare: `https://doi.org/10.6084/m9.figshare.27089614`.

The derived tables named below are preserved in the external archive
documented in `../../README_generated_outputs.md`.

Files:
- `mdd2025_gene_lists.tsv`: fastBAT Bonferroni-significant genes, Hi-C genes significant in all listed tissues/cell types, and DrugTargetor/MAGMA genes with qBH <= 0.05. The DrugTargetor/MAGMA qBH <= 0.05 list is marked as primary because it is the largest significant gene set.
- `mdd2025_variant_lists.tsv`: COJO selected independent SNPs and PolyFun credible causal variants. Multi-ancestry COJO and diverse credible-causal variants are marked primary; EUR-specific lists are retained.

Coordinates: source variant coordinates are hg19. `hg38_*` columns were filled by rsID and allele matching against `pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf` when possible.
