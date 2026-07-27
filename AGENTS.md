# Repository instructions

Treat this repository primarily as a reusable GWAS/PGS protocol, guidance, and
template project. Proactively tell the user when a proposed change mixes
generic workflow material with cohort-specific analysis or work products, and
recommend a clearer boundary.

Keep reusable code, documentation, configuration examples, and small test
fixtures in Git. Put generated datasets, figures, reports, logs, and other work
products in a designated ignored directory. Use `mbv-prs/generated/` for local
repo-adjacent products and documented external reference storage for large or
shared artifacts.

Never commit large binary or data files unless the user specifically asks for
those files to be committed. Keep GWAS, BCF, index, compressed data, LDGM, and
generated analysis outputs outside Git.
