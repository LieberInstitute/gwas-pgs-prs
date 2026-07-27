#!/usr/bin/env bash
set -euo pipefail

## verify both the 23andMe parser and inverse-variance meta-analysis arithmetic
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=/tmp/test-gwas-meta-integration.$$
mkdir -p "$work"
ref=/dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta
columns="$root/scripts/gwas_meta_colheaders.tsv"
export BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-/opt/sw/bcf-plugins}

cat > "$work/annotation.tsv" <<'EOF'
all.data.id	gt.data.id	im.data.id	assay.name	scaffold	position	alleles	ploidy	cytoband	gene.context	is.v1	is.v2	is.v3	is.v4	is.v5	h550	omni	strand
1	NA	1	rs3094315	chr1	752566	A/G	A	1p36	[]	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	+
2	NA	2	rs3131972	chr1	752721	A/G	A	1p36	[]	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	+
3	NA	3	test3	chr3	1000000	A/T	A	3p	[]	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	+
EOF

cat > "$work/association.tsv" <<'EOF'
all.data.id	src	pvalue	effect	stderr	pass	im.num.0	dose.b.0	im.num.1	dose.b.1	AA.0	AB.0	BB.0	AA.1	AB.1	BB.1
1	I	0.0455	-0.2	0.1	Y	200	0.2	100	0.3	NA	NA	NA	NA	NA	NA
2	I	0.0455	-0.2	0.1	Y	200	0.2	100	0.3	NA	NA	NA	NA	NA	NA
3	I	0.0455	-0.2	0.1	Y	200	0.2	100	0.3	NA	NA	NA	NA	NA	NA
EOF

paste "$work/annotation.tsv" "$work/association.tsv" | \
    awk -f "$root/scripts/prepare_23andme_gwas.awk" > "$work/23me.tsv"

cat > "$work/public.tsv" <<'EOF'
CHR	BP	SNP	A1	A2	BETA	SE	P	NCAS	NCON	NEFF
1	752566	rs3094315	A	G	0.1	0.2	0.6171	50	100	133.3333333
1	752721	rs3131972	A	G	0.1	0.2	0.6171	50	100	133.3333333
3	1000000	test3	A	T	0.1	0.2	0.6171	50	100	133.3333333
EOF

for component in public 23me; do
    bcftools +munge --no-version -Ob -o "$work/$component.bcf" \
        --write-index -C "$columns" -f "$ref" -s trait "$work/$component.tsv"
done

bcftools +metal --no-version -Ob -o "$work/meta.bcf" --write-index \
    "$work/public.bcf" "$work/23me.bcf"
bcftools query -f '%POS[\t%ES\t%SE\t%NE]\n' "$work/meta.bcf" > "$work/observed.tsv"

awk 'BEGIN{ok=1} {
    rows++
    expected_beta=($1==752721 ? -0.18 : 0.18)
    if (($2-expected_beta)^2 > 1e-10) ok=0
    if (($3-0.0894427191)^2 > 1e-10) ok=0
    if (($4-400)^2 > 1e-6) ok=0
} END{exit !(ok && rows==3)}' "$work/observed.tsv"

printf 'PASS: parser, allele orientation, beta, SE, and effective N\n'
