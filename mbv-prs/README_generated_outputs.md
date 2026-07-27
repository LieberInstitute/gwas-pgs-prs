# Generated outputs kept outside Git

Generated result tables, figures, reports, and similar analysis work products
are excluded from the repository. New repo-adjacent work products should be
written under the ignored directory:

```text
mbv-prs/generated/
```

Large or shared preserved MBv files are stored under:

```text
/home/gpertea/work/ref/GWAS/mbv-prs-generated-outputs/
```

That directory preserves repository-relative paths and contains `SHA256SUMS`
for integrity checks. Git retains reusable code, documentation, configuration,
and small test fixtures.
