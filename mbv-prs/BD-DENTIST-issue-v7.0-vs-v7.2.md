# BD Delivery: DENTIST and v7.0/v7.2 Audit

## Bottom Line

The delivered O'Connell 2025 European BD association file is the 23andMe-only
component, not the final full-European meta-analysis described in the paper's
data-availability statement. The matching European v7.0 annotation archive was
available in a separate, undocumented `7.0 Annotations` folder in the shared
Drive. It was recovered on 2026-07-30.

Our reconstructed fixed-effect integration can be used for a PRS sensitivity
analysis, but it is not yet an exact copy of the paper's final statistics. The
remaining material reason is that the paper's final post-meta-analysis
HRC-based DENTIST filtering has not been applied to the reconstruction.

The v7.0/v7.2 mapping concern has been resolved empirically across all
57,611,376 rows. IDs, rsIDs, chromosomes, alleles, ploidy, and strand are
identical. v7.0 positions require the coordinate conversion documented below.

## Archive Audit

The delivered European BD TAR is:

```text
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/
  Bipolar-Disorder-O_Connell-2025/
  OConnell_2025_bipolar_european-7.0.tar
```

It contains only:

```text
bipolar.dat.gz
bipolar.html
bipolar_hits_1.html
bipolar_hits_2.html
bipolar_hits_3.html
```

The BD study TAR does not contain annotations. The shared Drive separately
contains `7.0 Annotations/v7.0_europe.tar.gz`, but provides no manifest linking
the study folder to that annotation folder. The recovered archive is:

```text
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/
  7.0-Annotations/v7.0_europe.tar.gz
SHA-256: fcd43848d62d97c89f0eac2a0cce148675f7ea07f8e9a83a1eb37b1e7d1dba7d
```

It passed a complete `7z t` archive test and was extracted without overwriting
v7.2. The large text members were individually compressed in place:

```text
7.0-Annotations/v7.0_europe/all_snp_info.txt.gz
7.0-Annotations/v7.0_europe/gt_snp_stat.txt.gz
7.0-Annotations/v7.0_europe/im_snp_stat.txt.gz
7.0-Annotations/v7.0_europe/23andMe_GWAS_Results_v7.0.docx
7.0-Annotations/v7.0_europe/23andMe_Platform_Annotations_v7.0.docx
```

At the start of this audit, no byte-identical v7.0 document or data payload
existed under `23andMe_MDD_BD`. Existing generic annotation files decompress to
different v7.2 content. The current filesystem cannot establish whether v7.0
files had existed earlier and were overwritten. No supplied file is a DENTIST
exclusion list.

Pre-extraction hashes, raw/compressed payload hashes, and the complete
cross-release comparison are retained under
`7.0-Annotations/audit_2026-07-30/`.

## DOCX Audit

The matching documentation is:

```text
7.0-Annotations/v7.0_europe/23andMe_GWAS_Results_v7.0.docx
7.0-Annotations/v7.0_europe/23andMe_Platform_Annotations_v7.0.docx
```

The platform document states that the merged annotation has the same variant
count and ordering as the corresponding association result, and defines
`all.data.id` as the unique integer key. It does not specify whether its
`position` field is zero- or one-based and does not mention DENTIST.

A full v7.0-versus-v7.2 row comparison found:

| Check | Rows |
|---|---:|
| Annotation rows | 57,611,376 |
| Internal-ID mismatches | 0 |
| rsID mismatches | 0 |
| Chromosome mismatches | 0 |
| Allele mismatches | 0 |
| Ploidy or strand mismatches | 0 |
| Position differences | 54,940,299 |
| BD `pass=Y` autosomal SNVs retained by integration | 21,137,709 |
| Retained SNVs with position differences | 21,137,709 |

For every retained SNV, v7.0 is one base greater than the correct GRCh37 VCF
position. For example, v7.0 places `rs3131972` at 1:752722, while the GRCh37
reference/public BD BCF and v7.2 place it at 1:752721. The integration script
therefore applies an explicit v7.0 annotation-position offset of `-1` before
normalization. The v7.2 override remains only to reproduce the earlier
provisional run.

## DENTIST

The public no-23andMe BD summary file says DENTIST-detected variants were
filtered before release. O'Connell et al. describe a second relevant step for
the full analysis: after inverse-variance meta-analysis and effective-sample
size filtering, DENTIST was applied with ancestry-matched HRC v1.0 LD to remove
variants whose observed z scores were inconsistent with neighboring-LD
predictions. See the [paper methods](https://pmc.ncbi.nlm.nih.gov/articles/12163093/)
and [DENTIST documentation](https://github.com/Yves-CHEN/DENTIST).

Our reconstructed full BD meta-analysis combines the public no-23andMe
aggregate with the delivered 23andMe aggregate and applies the paper's 75%
effective-N rule. It has not received a new DENTIST pass after that
meta-analysis. The public component's prior filtering does not substitute for
post-meta-analysis QC of the newly combined statistics.

The paper and supplied files do not provide enough detail to reproduce the
final DENTIST pass exactly: the exclusion list, exact software version,
parameters, HRC files/build, and variant-matching details are absent.

## Delivery Discrepancy

The delivered European 23andMe BD component has 72,682 cases and 1,541,394
controls. The paper's full European analysis has 131,969 cases and 2,322,416
controls. The delivered file is therefore not the full meta-analysis.

The paper states that approved researchers can access full summary statistics
including all analyzed SNPs and samples in the GWAS meta-analyses. See the
[Nature data-availability statement](https://www.nature.com/articles/s41586-024-08468-9).
The current delivery does not match that description.

Until corrected, use these labels:

- `BD no23andMe`: public release, already DENTIST filtered;
- `BD full pre-DENTIST`: reconstructed public plus 23andMe fixed-effect result;
- do not label the reconstruction `final paper GWAS`.

## Draft Email to 23andMe

**Subject:** O'Connell 2025 BD request: final full European meta-analysis

Hello,

We received the approved O'Connell 2025 bipolar-disorder files. The European
file `OConnell_2025_bipolar_european-7.0/bipolar.dat.gz` appears to contain the
23andMe-only component: 72,682 cases and 1,541,394 controls. The paper reports
131,969 cases and 2,322,416 controls for the full European meta-analysis, and
its data-availability statement says approved researchers can access full
summary statistics including all analyzed SNPs and samples.

Could you please provide the final DENTIST-filtered full-European meta-analysis
summary statistics used in O'Connell et al. 2025, including genomic build,
effect/non-effect alleles, beta or odds ratio, standard error, p value, and
variant-specific effective sample size?

Please also confirm whether any DENTIST exclusion list or other final
post-meta-analysis QC file should have been included in the delivery.

Thank you.
## Draft Email to PGC/Paper Team

**Subject:** Request for final European BD summary statistics or DENTIST QC details

Hello,

We are using the European O'Connell et al. 2025 bipolar-disorder GWAS for PRS
and eQTL/GWAS colocalization. Our approved 23andMe delivery contains the
23andMe-only European component rather than the final full-European
meta-analysis. We can reconstruct the inverse-variance meta-analysis with the
public no-23andMe statistics, but cannot reproduce the paper's final DENTIST
QC exactly.

Could you please provide either:

1. the final DENTIST-filtered full-European summary-statistics file used for
   the paper; or
2. the DENTIST exclusion list plus the exact software version, command-line
   parameters, HRC v1.0 reference files and build, variant-matching rules, and
   any additional post-meta-analysis filters?

Thank you.
