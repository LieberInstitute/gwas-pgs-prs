# BPD 2024 public significant locus/signal extracts

Source: Nature supplementary workbook for OConnell et al. 2025, doi:10.1038/s41586-024-08468-9. Coordinates are hg19 as reported by the supplement.

The extracts named below are preserved in the external archive documented in
`../README_generated_outputs.md`.

| Extract | Rows/signals | Unique loci | File |
|---|---:|---:|---|
| EUR including self-report/23andMe, Table S5 | 261 | 229 | `bpd2024_eur_with23andMe_gws_loci_signals.csv` |
| EUR excluding self-report/23andMe, Table S6 | 94 | 88 | `bpd2024_eur_no23andMe_gws_loci_signals.csv` |
| multi-ancestry including self-report/23andMe, Table S8 | 337 | 298 | `bpd2024_multiancestry_with23andMe_gws_loci_signals.csv` |
| multi-ancestry excluding self-report/23andMe, Table S11 | 116 | 105 | `bpd2024_multiancestry_no23andMe_gws_loci_signals.csv` |
| self-report only, Table S12 | 126 | 126 | `bpd2024_selfreport_only_gws_loci_signals.csv` |

`bpd2024_public_gws_signals_minimal.tsv` combines the main signal columns across these extracts.
`bpd2024_multiancestry_prioritized_credible_genes.csv` extracts Table S31 credible gene prioritization.

These are lead/significant loci and independent signal tables, not full GWAS summary statistics and not every raw P <= 5e-8 variant.
