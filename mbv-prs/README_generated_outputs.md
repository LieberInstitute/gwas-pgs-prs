# Generated outputs kept outside Git

Generated datasets, figures, reports, logs, and similar analysis work products
are excluded from the repository. New workflows should write repo-adjacent
work products under the ignored directory:

```text
mbv-prs/generated/
```

Existing reproducibility scripts retain their documented output paths, such as
`qc/`, `scores_bcftools/`, `prsice/out/`, and `prs_diagnosis_assoc/`. These are
also designated work-product directories and are ignored by Git.

Large or shared preserved MBv files are stored under:

```text
/home/gpertea/work/ref/GWAS/mbv-prs-generated-outputs/
```

That directory preserves repository-relative paths and contains `SHA256SUMS`
for integrity checks. It includes legacy outputs removed from Git when this
policy was adopted. Git retains reusable code, documentation, configuration,
and small test fixtures.
