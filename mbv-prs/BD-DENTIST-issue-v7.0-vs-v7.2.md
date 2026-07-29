# BD Delivery: DENTIST and v7.0/v7.2 Issues

## Bottom Line

The delivered O'Connell 2025 European BD association file is the 23andMe-only
component, not the final full-European meta-analysis described in the paper's
data-availability statement. It is labeled release v7.0, but the delivery does
not include the corresponding European v7.0 annotation tables. The only
European annotation tables supplied are v7.2.

Our reconstructed fixed-effect integration can be used for a PRS sensitivity
analysis, but it is not yet an exact copy of the paper's final statistics for
two independent reasons:

1. cross-release `all.data.id` and allele stability is not formally documented;
2. the paper's final post-meta-analysis HRC-based DENTIST filtering has not been
   applied to the reconstructed meta-analysis.

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

The original delivery ZIPs contain the European v7.2 annotation package and
BD association TARs. They do not contain European v7.0 `all_snp_info`,
`gt_snp_stat`, or `im_snp_stat` annotation files. A downloaded v7.0 Africa
annotation archive is irrelevant to the European association file.

The available European annotation files are:

```text
7.2-Annotations/v7.2_europe/all_snp_info.txt.gz
7.2-Annotations/v7.2_europe/gt_snp_stat.txt.gz
7.2-Annotations/v7.2_europe/im_snp_stat.txt.gz
```

No delivered file is a DENTIST exclusion list.

## DOCX Audit

The only relevant annotation documentation is also v7.2:

```text
7.2-Annotations/v7.2_europe/23andMe_GWAS_Results_v7.2.docx
7.2-Annotations/v7.2_europe/23andMe_Platform_Annotations_v7.2.docx
```

These documents state that association rows align with their corresponding
merged annotation release and that `all.data.id` is the unique integer key.
They do not state that `all.data.id`, coordinates, or alleles are stable across
v7.0 and v7.2. They do not mention DENTIST.

The reconstructed BD integration therefore requires an explicit
`--allow-bd-v7.2-annotations` acknowledgement. Matching checked paper sentinels
supports the mapping empirically, but does not replace release provenance.

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

**Subject:** O'Connell 2025 BD request: full European meta-analysis and v7.0 annotations

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

If that final file is not available through 23andMe, please provide the
European v7.0 annotation files corresponding to `bipolar.dat.gz`. The delivery
contains European v7.2 annotations only. Alternatively, please provide written
confirmation that `all.data.id`, coordinates, allele labels, and association
row ordering are unchanged between European v7.0 and v7.2.

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

The 23andMe association file is labeled v7.0 while the supplied European
annotations are v7.2. Any confirmation of the annotation release used in the
paper or cross-release ID/allele stability would also be helpful.

Thank you.
