# BPD 2024 gene and variant lists

Source paper: OConnell et al. Nature 2025 doi:10.1038/s41586-024-08468-9.
Source supplement: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-08468-9/MediaObjects/41586_2024_8468_MOESM4_ESM.xlsx`.

The derived tables named below are preserved in the external archive
documented in `../../README_generated_outputs.md`.

Files:
- `bpd2024_gene_lists.tsv`: Table S31 multi-ancestry prioritized credible genes. This is the primary BPD gene list.
- `bpd2024_variant_lists.tsv`: published locus/signal rows from Tables S5, S6, S8, S11, and S12. `multiancestry_with23andMe` is marked as the primary variant list; EUR and no-23andMe lists are retained.

Coordinates: source variant coordinates are hg19. `hg38_*` columns were filled by rsID and allele matching against `bip2024_eur_no23andMe.hg38.bcf` when possible.
