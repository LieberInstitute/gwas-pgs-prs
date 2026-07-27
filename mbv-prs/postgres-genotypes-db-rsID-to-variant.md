# Postgres genotypes database rsID and variant mapping

The local Postgres database `genotypes` stores hg38 dbSNP mappings between dbSNP rsIDs and canonical variant identifiers of the form `chr:pos:REF:ALT`.

Connect with:

```bash
psql -U gpertea -d genotypes
```

## Objects

`public.chr_code_map`

Maps compact chromosome codes to canonical chromosome names.

| column | type | meaning |
|---|---|---|
| `chr_code` | `smallint` | internal code |
| `chr_name` | `text` | canonical name, for example `chr1` |
| `is_canonical` | `boolean` | true for canonical names |

`public.rsid_hg38`

Stores hg38 rsID to allele mappings. `rsid_num` is the numeric part of an rsID, so `rs12238997` is stored as `12238997`.

| column | type | meaning |
|---|---|---|
| `rsid_num` | `bigint` | numeric rsID |
| `chr_code` | `smallint` | chromosome code |
| `pos` | `integer` | 1-based hg38 position |
| `ref` | `text` | hg38 reference allele |
| `alt` | `text` | alternate allele |

Important indexes:

```text
rsid_hg38_idx_rsid    on rsid_hg38(rsid_num)
rsid_hg38_idx_chr_pos on rsid_hg38(chr_code, pos)
rsid_hg38_uq_row      on rsid_hg38(rsid_num, chr_code, pos, ref, alt)
```

Large indexes should not be built in the default tablespace because the root partition is space-constrained. Use a tablespace rooted under `/ssdata/postgres/`, for example `ts_ssdata_indexes`, for new or moved large indexes.

Recommended setup:

```bash
sudo mkdir -p /ssdata/postgres/pg_ts_indexes
sudo chown postgres:postgres /ssdata/postgres /ssdata/postgres/pg_ts_indexes
sudo chmod 700 /ssdata/postgres/pg_ts_indexes
```

Then create and use the tablespace:

```sql
create tablespace ts_ssdata_indexes location '/ssdata/postgres/pg_ts_indexes';
alter index public.rsid_hg38_idx_chr_pos set tablespace ts_ssdata_indexes;
```

`public.rsid_hg38_decoded`

View with decoded names:

| column | meaning |
|---|---|
| `rsid` | `rs` plus `rsid_num` |
| `varid` | `chr:pos:REF:ALT` |
| `chr`, `pos`, `ref`, `alt` | decoded coordinates and alleles |
| `rsid_num`, `chr_code` | stored numeric keys |

## Functions

`normalize_chr_name(p_chr text) returns text`

Normalizes chromosome input. Examples: `1` and `chr1` become `chr1`; `x` becomes `chrX`; `mt` becomes `chrM`. Empty input raises an error.

`chr_to_code(p_chr text) returns smallint`

Returns the internal chromosome code after chromosome normalization.

`code_to_chr(p_code smallint) returns text`

Returns the canonical chromosome name for an internal code.

`rsid_to_num(p_rsid text) returns bigint`

Accepts `rs123`, `RS123`, or `123`, and returns `123`. Invalid input raises `Invalid rsID`.

`num_to_rsid(p_num bigint) returns text`

Returns `rs` plus the numeric rsID.

`make_varid(p_chr_code smallint, p_pos integer, p_ref text, p_alt text) returns text`

Builds `chr:pos:REF:ALT` from stored fields.

`lookup_rsid(p_rsid text)`

Returns all hg38 canonical variants for one rsID:

```sql
select *
from public.lookup_rsid('rs12238997');
```

Returns columns: `rsid`, `varid`, `chr`, `pos`, `ref`, `alt`.

`lookup_interval(p_chr text, p_start integer, p_end integer)`

Returns all rsID mappings in a genomic interval:

```sql
select *
from public.lookup_interval('chr1', 758000, 759000);
```

## Single conversions

rsID to canonical variant:

```sql
select rsid, varid, chr, pos, ref, alt
from public.lookup_rsid('rs12238997');
```

Canonical variant to rsID:

```sql
with v(chr, pos, ref, alt) as (
  values ('chr1', 758351, 'A', 'G')
)
select public.num_to_rsid(r.rsid_num) as rsid,
       public.make_varid(r.chr_code, r.pos, r.ref, r.alt) as varid
from v
join public.rsid_hg38 r
  on r.chr_code = public.chr_to_code(v.chr)
 and r.pos = v.pos
 and r.ref = v.ref
 and r.alt = v.alt;
```

REF/ALT swap check:

```sql
with v(chr, pos, ref, alt) as (
  values ('chr1', 758351, 'A', 'G')
)
select public.num_to_rsid(r.rsid_num) as rsid,
       public.make_varid(r.chr_code, r.pos, r.ref, r.alt) as db_varid,
       case
         when r.ref = v.ref and r.alt = v.alt then 'exact'
         when r.ref = v.alt and r.alt = v.ref then 'ref_alt_swap'
       end as match_type
from v
join public.rsid_hg38 r
  on r.chr_code = public.chr_to_code(v.chr)
 and r.pos = v.pos
 and ((r.ref = v.ref and r.alt = v.alt) or
      (r.ref = v.alt and r.alt = v.ref));
```

## Bulk rsID to canonical variants

Prepare a one-column text file:

```text
rsid
rs12238997
rs4951859
```

Load and convert:

```sql
create temp table input_rsid(rsid text);
\copy input_rsid from 'rsids.tsv' with (format csv, delimiter E'\t', header true)

select i.rsid as input_rsid,
       public.num_to_rsid(r.rsid_num) as rsid,
       public.make_varid(r.chr_code, r.pos, r.ref, r.alt) as varid,
       public.code_to_chr(r.chr_code) as chr,
       r.pos,
       r.ref,
       r.alt
from input_rsid i
left join public.rsid_hg38 r
  on r.rsid_num = public.rsid_to_num(i.rsid)
order by i.rsid, r.chr_code, r.pos, r.ref, r.alt;
```

One rsID can map to multiple variants. Keep all rows unless a downstream step requires biallelic SNP-only records.

## Bulk canonical variants to rsIDs

Prepare a four-column file:

```text
chr	pos	ref	alt
chr1	758351	A	G
chr1	794299	C	G
```

Load and convert:

```sql
create temp table input_variant(
  chr text,
  pos integer,
  ref text,
  alt text
);
\copy input_variant from 'variants.tsv' with (format csv, delimiter E'\t', header true)

select i.chr || ':' || i.pos || ':' || i.ref || ':' || i.alt as input_varid,
       public.num_to_rsid(r.rsid_num) as rsid,
       public.make_varid(r.chr_code, r.pos, r.ref, r.alt) as db_varid
from input_variant i
left join public.rsid_hg38 r
  on r.chr_code = public.chr_to_code(i.chr)
 and r.pos = i.pos
 and r.ref = i.ref
 and r.alt = i.alt
order by input_varid, rsid;
```

To flag REF/ALT swaps in bulk:

```sql
select i.chr || ':' || i.pos || ':' || i.ref || ':' || i.alt as input_varid,
       public.num_to_rsid(r.rsid_num) as rsid,
       public.make_varid(r.chr_code, r.pos, r.ref, r.alt) as db_varid,
       case
         when r.ref = i.ref and r.alt = i.alt then 'exact'
         when r.ref = i.alt and r.alt = i.ref then 'ref_alt_swap'
       end as match_type
from input_variant i
left join public.rsid_hg38 r
  on r.chr_code = public.chr_to_code(i.chr)
 and r.pos = i.pos
 and ((r.ref = i.ref and r.alt = i.alt) or
      (r.ref = i.alt and r.alt = i.ref));
```

## Pipeline use in this repository

The MBv target genotype VCF uses canonical IDs in the VCF `ID` column and stores numeric dbSNP IDs in `INFO/RS`. The PRS pipeline now writes alignment outputs under:

```text
qc/gwas_variant_alignment/
```

For each disorder, `${label}.canonical.tsv.gz` contains canonical GWAS rows with:

```text
disorder CHR BP SNP RSID SOURCE_ID REF ALT A1 A2 BETA P LP db_exact_match db_same_position db_ref_alt_swap target_exact_match target_rsid target_rsid_same_as_gwas target_db_exact_match
```

PRSice base files are rebuilt from rows where the GWAS rsID maps exactly to the same hg38 `chr:pos:REF:ALT` in Postgres and no REF/ALT swap candidate is present. `A1` remains `ALT`, `A2` remains `REF`, and `BETA` remains the GWAS-VCF `ES` value relative to `ALT`.

## 2026-05-12 implementation log

Today the broken `rsid_hg38_idx_chr_pos` index was repaired and moved to a dedicated non-root tablespace under `/ssdata/postgres/`. The PRS pipeline was also updated to validate GWAS rsIDs against the Postgres hg38 mapping table before creating PRSice base files.

Files added or changed:

```text
postgres-genotypes-db-rsID-to-variant.md
scripts/align_gwas_variants.sh
scripts/run_mbv_prs.sh
plan_26-05-12_09-19_master.md
qc/gwas_variant_alignment/
prsice/base/*.base.tsv.gz
prsice/out/*
mbv_prs_merged_scores.tsv
```

Postgres database and index checks:

```bash
psql -U gpertea -d genotypes -Atc "select n.nspname||'.'||p.proname||'('||pg_get_function_arguments(p.oid)||') returns '||pg_get_function_result(p.oid) from pg_proc p join pg_namespace n on n.oid=p.pronamespace where n.nspname not in ('pg_catalog','information_schema') and (p.proname ilike '%rs%' or p.proname ilike '%snp%' or p.proname ilike '%variant%' or p.proname ilike '%dbsnp%' or p.proname ilike '%chr%') order by 1;"
psql -U gpertea -d genotypes -Atc "select table_schema||'.'||table_name||E'\t'||column_name||E'\t'||data_type||coalesce('('||character_maximum_length||')','')||E'\t'||is_nullable from information_schema.columns where table_schema='public' order by table_name, ordinal_position;"
psql -U gpertea -d genotypes -Atc "select schemaname||'.'||tablename||E'\t'||indexname||E'\t'||indexdef from pg_indexes where schemaname='public' order by tablename,indexname;"
psql -U gpertea -d genotypes -Atc "select relname, oid, relfilenode, pg_relation_filepath(oid), pg_size_pretty(pg_total_relation_size(oid)) from pg_class where relnamespace='public'::regnamespace order by relname;"
psql -U gpertea -d genotypes -Atc "select * from public.lookup_rsid('rs12238997') limit 5;"
```

Initial repair attempt and issue:

```bash
psql -U gpertea -d genotypes -c "REINDEX INDEX public.rsid_hg38_idx_chr_pos;"
```

That failed because the old tablespace path for `rsid_hg38_idx_chr_pos` was missing:

```text
could not create directory "pg_tblspc/61718/PG_17_202406281/61672": No such file or directory
```

The index was then dropped and recreated, then moved off the root-backed default tablespace:

```bash
psql -U gpertea -d genotypes -c "DROP INDEX IF EXISTS public.rsid_hg38_idx_chr_pos; CREATE INDEX rsid_hg38_idx_chr_pos ON public.rsid_hg38 USING btree (chr_code, pos);"
psql -U gpertea -d genotypes -c "ALTER INDEX public.rsid_hg38_idx_chr_pos SET TABLESPACE ts_sdb1;"
```

After `/ssdata/postgres/pg_ts_indexes` was created as root with owner `postgres:postgres`, the final tablespace and move were run:

```bash
psql -U gpertea -d genotypes -v ON_ERROR_STOP=1 -c "CREATE TABLESPACE ts_ssdata_indexes LOCATION '/ssdata/postgres/pg_ts_indexes';" -c "ALTER INDEX public.rsid_hg38_idx_chr_pos SET TABLESPACE ts_ssdata_indexes;"
```

Final tablespace verification:

```bash
psql -U gpertea -d genotypes -Atc "select oid, spcname, pg_tablespace_location(oid) from pg_tablespace order by oid; select c.relname, t.spcname, pg_relation_filepath(c.oid), pg_size_pretty(pg_relation_size(c.oid)) from pg_class c left join pg_tablespace t on t.oid=c.reltablespace where c.relname='rsid_hg38_idx_chr_pos'; select * from public.lookup_rsid('rs12238997') limit 5;"
```

Final verified state:

```text
ts_ssdata_indexes -> /ssdata/postgres/pg_ts_indexes
rsid_hg38_idx_chr_pos -> ts_ssdata_indexes
index size -> 25 GB
rs12238997 -> chr1:758351:A:G
```

GWAS and genotype validation commands:

```bash
for f in /home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf /home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf /home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf; do printf '%s\n' "$f"; /opt/sw/bin/bcftools norm -f /dbdata/cdb/genotyping/ref/hg38c.fa -c e -Ou "$f" >/dev/null; done
/opt/sw/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' qc/mbv.qc.chr.bcf | awk 'BEGIN{bad=0;n=0} {n++; if($3 != $1":"$2":"$4":"$5) bad++} END{print "target_rows",n,"bad_canonical_ids",bad}'
```

Validation results:

```text
BPD GWAS REF check: 0 mismatches
MDD GWAS REF check: 0 mismatches
SCZD GWAS REF check: 0 mismatches
target_rows 6639276 bad_canonical_ids 0
```

Full GWAS-to-Postgres alignment command:

```bash
rm -rf qc/gwas_variant_alignment
BCF=/opt/sw/bin/bcftools TARGET_BCF=qc/mbv.qc.chr.bcf scripts/align_gwas_variants.sh BPD /home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf MDD /home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf SCZD /home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf
cat qc/gwas_variant_alignment/alignment_summary.tsv
```

Alignment summary:

```text
disorder	metric	value
BPD	db_exact_match	6937793
BPD	db_ref_alt_swap	0
BPD	db_same_position	6937796
BPD	non_rsid_source_id	0
BPD	rsid_source_id	6938758
BPD	target_db_exact_match	6065249
BPD	target_exact_match	6065997
BPD	target_rsid_present	6065249
BPD	target_rsid_same_as_gwas	6062290
BPD	total_rows	6938758
MDD	db_exact_match	7348190
MDD	db_ref_alt_swap	0
MDD	db_same_position	7348190
MDD	non_rsid_source_id	576
MDD	rsid_source_id	7362098
MDD	target_db_exact_match	5859109
MDD	target_exact_match	5859844
MDD	target_rsid_present	5859109
MDD	target_rsid_same_as_gwas	5845877
MDD	total_rows	7362674
SCZD	db_exact_match	7623093
SCZD	db_ref_alt_swap	0
SCZD	db_same_position	7623108
SCZD	non_rsid_source_id	20934
SCZD	rsid_source_id	7637433
SCZD	target_db_exact_match	6258203
SCZD	target_exact_match	6259027
SCZD	target_rsid_present	6258203
SCZD	target_rsid_same_as_gwas	6243657
SCZD	total_rows	7658367
```

PRSice base regeneration command used after alignment:

```bash
for label in BPD MDD SCZD; do
  out=prsice/base/${label}.base.tsv
  printf 'CHR\tBP\tSNP\tA1\tA2\tBETA\tP\n' > "$out"
  gzip -cd qc/gwas_variant_alignment/${label}.canonical.tsv.gz | awk -F '\t' 'BEGIN{OFS="\t"} NR>1 && $14=="t" && $16!="t" {chr=$2; sub(/^chr/,"",chr); if(!seen[$4]++) print chr,$3,$4,$9,$10,$11,$12}' >> "$out"
  gzip -f "$out"
done
```

Validated PRSice base row counts, excluding the header:

```text
prsice/base/BPD.base.tsv.gz	6937793
prsice/base/MDD.base.tsv.gz	7348190
prsice/base/SCZD.base.tsv.gz	7623093
```

PRSice and score summarization command:

```bash
THREADS=${THREADS:-8}
PRSICE=/opt/sw/bin/PRSice_linux
LDREF=/dbdata/cdb/ref/1000g_hg38/1KG_EUR_chrpos
THRESHOLDS=5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1
for label in BPD MDD SCZD; do
  "$PRSICE" --base "prsice/base/${label}.base.tsv.gz" --target prsice/target/mbv_qc --type bgen --ld "$LDREF" --ld-type bed --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 --stat BETA --pvalue P --beta --ignore-fid --pheno prsice/target/dummy.pheno --pheno-col DUMMY --no-regress --all-score --score sum --fastscore --bar-levels "$THRESHOLDS" --clump-kb 250 --clump-r2 0.1 --clump-p 1 --thread "$THREADS" --out "prsice/out/${label}"
done
Rscript scripts/summarize_mbv_prs.R
```

PRSice run results:

```text
BPD observed in base: 6937793
BPD ambiguous excluded: 1046953
BPD variants after clumping: 224355
MDD observed in base: 7348190
MDD ambiguous excluded: 950512
MDD variants after clumping: 231761
SCZD observed in base: 7623093
SCZD ambiguous excluded: 1151139
SCZD variants after clumping: 237767
```

Output table checks:

```bash
for f in prsice/out/BPD.all_score prsice/out/MDD.all_score prsice/out/SCZD.all_score prsice/out/mbv_prsice_scores.tsv prsice/out/mbv_prsice_scores_z.tsv mbv_prs_merged_scores.tsv; do printf '%s\t' "$f"; awk 'END{print NR}' "$f"; done
```

All checked PRSice and merged score tables had 120 lines: one header plus 119 MBv samples.
