#!/usr/bin/env bash
set -euo pipefail

BCF=${BCF:-/opt/sw/bin/bcftools}
PSQL=${PSQL:-psql}
DB=${DB:-genotypes}
DBUSER=${DBUSER:-gpertea}
OUTDIR=${OUTDIR:-qc/gwas_variant_alignment}
TARGET_BCF=${TARGET_BCF:-qc/mbv.qc.chr.bcf}

mkdir -p "$OUTDIR"
TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR"

if [ "$#" -lt 2 ] || [ $(( $# % 2 )) -ne 0 ]; then
  printf 'usage: %s LABEL GWAS_BCF [LABEL GWAS_BCF ...]\n' "$0" >&2
  exit 2
fi

target_tsv="$TMPDIR/target.tsv"

## write final target variants with canonical IDs and numeric RS annotations
"$BCF" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/RS\n' "$TARGET_BCF" | \
  awk 'BEGIN{OFS="\t"} {rs=($6=="." ? "" : "rs"$6); print $1,$2,$3,$4,$5,rs}' > "$target_tsv"

summary_all="$OUTDIR/alignment_summary.tsv"
printf 'disorder\tmetric\tvalue\n' > "$summary_all"

while [ "$#" -gt 0 ]; do
  label=$1
  gwas=$2
  shift 2

  gwas_tsv="$TMPDIR/${label}.gwas.tsv"
  canonical_tsv="$OUTDIR/${label}.canonical.tsv"
  summary_tsv="$TMPDIR/${label}.summary.tsv"

  ## extract biallelic SNP summary statistics with ES relative to ALT
  "$BCF" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%LP]\n' \
    -i 'N_ALT=1 && TYPE="snp" && FILTER!="IFFY" && FILTER!="REF_MISMATCH"' "$gwas" | \
    awk 'BEGIN{OFS="\t"} $6!="." && $7!="." {print}' > "$gwas_tsv"

  "$PSQL" -U "$DBUSER" -d "$DB" -v ON_ERROR_STOP=1 <<SQL
create temp table target_raw(
  chr text,
  pos integer,
  varid text,
  ref text,
  alt text,
  target_rsid text
);
\copy target_raw from '${target_tsv}' with (format csv, delimiter E'\t', null '')

create temp table target_map as
select chr,
       public.chr_to_code(chr) as chr_code,
       pos,
       varid,
       ref,
       alt,
       nullif(target_rsid, '') as target_rsid,
       case
         when target_rsid ~* '^rs[0-9]+$' then public.rsid_to_num(target_rsid)
       end as target_rsid_num
from target_raw;
create index target_map_varid_idx on target_map(varid);
create index target_map_rsid_pos_idx on target_map(target_rsid_num, chr_code, pos);

create temp table gwas_raw(
  chr text,
  pos integer,
  source_id text,
  ref text,
  alt text,
  beta double precision,
  lp double precision
);
\copy gwas_raw from '${gwas_tsv}' with (format csv, delimiter E'\t', null '')

create temp table gwas_map as
select row_number() over () as row_id,
       '${label}'::text as disorder,
       chr,
       public.chr_to_code(chr) as chr_code,
       pos,
       chr || ':' || pos::text || ':' || ref || ':' || alt as varid,
       case
         when source_id ~* '^rs[0-9]+$' then source_id
       end as rsid,
       source_id,
       ref,
       alt,
       alt as a1,
       ref as a2,
       beta,
       power(10::double precision, -lp) as p,
       lp,
       case
         when source_id ~* '^rs[0-9]+$' then public.rsid_to_num(source_id)
       end as rsid_num
from gwas_raw;
create index gwas_map_row_idx on gwas_map(row_id);
create index gwas_map_rsid_pos_idx on gwas_map(rsid_num, chr_code, pos);
create index gwas_map_varid_idx on gwas_map(varid);

create temp table db_match as
select g.row_id,
       bool_or(r.ref = g.ref and r.alt = g.alt) as db_exact_match,
       bool_or(r.rsid_num is not null) as db_same_position,
       bool_or(r.ref = g.alt and r.alt = g.ref) as db_ref_alt_swap
from gwas_map g
left join public.rsid_hg38 r
  on r.rsid_num = g.rsid_num
 and r.chr_code = g.chr_code
 and r.pos = g.pos
group by g.row_id;
create index db_match_row_idx on db_match(row_id);

create temp table aligned as
select g.disorder,
       g.chr as "CHR",
       g.pos as "BP",
       g.varid as "SNP",
       coalesce(g.rsid, '') as "RSID",
       g.source_id as "SOURCE_ID",
       g.ref as "REF",
       g.alt as "ALT",
       g.a1 as "A1",
       g.a2 as "A2",
       g.beta as "BETA",
       g.p as "P",
       g.lp as "LP",
       coalesce(d.db_exact_match, false) as db_exact_match,
       coalesce(d.db_same_position, false) as db_same_position,
       coalesce(d.db_ref_alt_swap, false) as db_ref_alt_swap,
       (t.varid is not null) as target_exact_match,
       coalesce(t.target_rsid, '') as target_rsid,
       (g.rsid is not null and t.target_rsid = g.rsid) as target_rsid_same_as_gwas,
       (tr.rsid_num is not null) as target_db_exact_match
from gwas_map g
join db_match d on d.row_id = g.row_id
left join target_map t on t.varid = g.varid
left join public.rsid_hg38 tr
  on tr.rsid_num = t.target_rsid_num
 and tr.chr_code = t.chr_code
 and tr.pos = t.pos
 and tr.ref = t.ref
 and tr.alt = t.alt;

\copy aligned to '${canonical_tsv}' with (format csv, delimiter E'\t', header true, null '')

create temp table alignment_summary as
select '${label}'::text as disorder, metric, value::bigint
from (
  select 'total_rows' as metric, count(*) as value from aligned
  union all select 'rsid_source_id', count(*) from aligned where "RSID" <> ''
  union all select 'non_rsid_source_id', count(*) from aligned where "RSID" = ''
  union all select 'db_exact_match', count(*) from aligned where db_exact_match
  union all select 'db_same_position', count(*) from aligned where db_same_position
  union all select 'db_ref_alt_swap', count(*) from aligned where db_ref_alt_swap
  union all select 'target_exact_match', count(*) from aligned where target_exact_match
  union all select 'target_rsid_present', count(*) from aligned where target_rsid <> ''
  union all select 'target_rsid_same_as_gwas', count(*) from aligned where target_rsid_same_as_gwas
  union all select 'target_db_exact_match', count(*) from aligned where target_db_exact_match
) s
order by metric;
\copy alignment_summary to '${summary_tsv}' with (format csv, delimiter E'\t', header false)
SQL

  gzip -f "$canonical_tsv"
  cat "$summary_tsv" >> "$summary_all"
done

rm -rf "$TMPDIR"
