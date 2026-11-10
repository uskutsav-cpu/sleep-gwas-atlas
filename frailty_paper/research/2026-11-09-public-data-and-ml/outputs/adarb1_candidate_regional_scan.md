# Regional association context for the ADARB1-labeled candidate

I scanned published summary-statistic rows within 250 kb of rs2838788 in each of four available GWAS. **The selected candidate is not the lowest-p row in this window in any of these files.** None of the smallest displayed p-values reaches the conventional genome-wide threshold of 5 × 10⁻⁸. This makes the candidate useful for follow-up, but these regional marginal associations do not establish a genome-wide significant independent locus, a shared causal variant, or ADARB1 as its effector gene.

| Source | Valid p-value rows in window | rs2838788 p | Smallest p in that source window | Marker for smallest p | Distance from candidate |
|---|---:|---:|---:|---|---:|
| Mak FI | 2,061 | 2.62 × 10⁻⁵ | 2.07 × 10⁻⁵ | rs8129315 | −43,511 bp |
| Mak physical phenotype | 2,061 | 0.00813 | 2.38 × 10⁻⁵ | rs35424628 | +159,064 bp |
| MVP insomnia | 4,880 | 0.001884 | 1.80 × 10⁻⁴ | rs188632112 | −19,299 bp |
| FinnGen R13 insomnia | 4,772 | 0.000345 | 1.16 × 10⁻⁶ | rs1320256433 | −99,660 bp |

These are the numerically smallest **marginal source p-values**, not independently confirmed or QC-certified lead variants. Variant coverage, imputation and quality filters can differ; this scan required only a finite p in [0,1] and did not compare matching SNP sets, inspect per-variant quality, estimate LD with the candidate, or run conditional analysis. In particular, a smaller p-value in a different file is not evidence for a different causal signal. None of the comparisons was prespecified independently of the candidate selection.

### Reference LD for the two Mak minimum-p markers

A subsequent allele-specific check used 503 European-labeled 1000 Genomes phase 3 individuals. For the **specific tested allele pairs**, rs2838788 T and the Mak-FI minimum-p marker rs8129315 G had phased-haplotype **r² = 0.951** (r = +0.975); rs2838788 T and the Mak-physical-phenotype minimum-p marker rs35424628 G had **r² = 0.042** (r = −0.204). All three position/ID/biallelic allele-pair records matched in this reference and all 503 individuals had phased calls. The two minimum-p alleles are the Mak A1 alleles. This is consistent with substantial tagging between the FI minimum-p marker and rs2838788 in this panel, and less pairwise tagging between the physical-phenotype minimum-p marker and rs2838788.

Neither comparison proves that the signals share or differ in their causal variants. In particular, low pairwise LD does not establish conditional independence, and the 1000 Genomes panel is not the Mak discovery sample. Different variants and other intermediates can influence a regional signal. The corresponding insomnia minimum-p markers were **not** tested in this LD check. [Reference LD receipt](adarb1_mak_regional_ld.json); calculation: `work/adarb1_frailty_pair_ld.py`. Source rows for rs8129315 and rs35424628 were checked against Mak's released summary files. This check used reference genotypes only and saved no individual-level data.

### Both frailty outcomes at the same three markers

Direct lookup in both Mak files shows nominal p < 0.05 for each of the three specific allele contrasts in **both** frailty outcomes. Thus the fact that different markers minimize p in each file does not mean that either marker is specific to only one phenotype. The table uses the same source A1 allele across the two outcomes; FI and physical-phenotype coefficients use different scales and must not be compared by magnitude.

| Marker, effect / other allele | Mak FI beta (SE), p | Mak physical beta (SE), p |
|---|---:|---:|
| rs8129315, G / A | +0.07135 (0.01676), 2.07 × 10⁻⁵ | +0.00528 (0.00188), 0.00494 |
| rs2838788, T / C | +0.07063 (0.01680), 2.62 × 10⁻⁵ | +0.00498 (0.00188), 0.00813 |
| rs35424628, G / A | −0.09745 (0.03200), 0.00233 | −0.01515 (0.00358), 2.38 × 10⁻⁵ |

These tests and marker choices are post hoc and uncorrected. The identical direction across these two frailty outcomes for a given allele does not establish association with insomnia or a shared causal signal. At rs35424628 the G allele is the more frequent Mak effect allele (AF 0.930); its negative coefficient is a contrast with A, not proof that manipulating G protects anyone. Allele frequency, different LD and phenotype scales complicate interpretation of raw effect magnitudes. [Released Mak files and checksums](mak2023_same_cohort_fi_fp_ldsc_review.md).

The regional insomnia scans also contain rs8129315 at the corresponding reported chromosome and position: MVP p = 0.00216 and FinnGen p = 0.000462. For rs35424628 the corresponding p-values are 0.116 and 0.177. In particular, the rs8129315 and rs2838788 observations are highly correlated in the reference panel and are not independent confirmations.

I subsequently checked the raw insomnia rows and harmonized both markers to **G versus A**: Mak and MVP report G as effect allele, and FinnGen reports A as alternate/effect allele, so its beta is reversed for the comparison. The assumed Mak GRCh37 coordinates map uniquely to the reported GRCh38 coordinates on the positive strand, and the allele pairs match; these checks do not independently authenticate all source coding or reference representations.

| Marker | MVP insomnia OR on G (p) | FinnGen insomnia OR on G, approximate 95% CI (p) | Direction versus both Mak frailty outcomes |
|---|---:|---:|---|
| rs8129315 | 1.018 (0.00216) | 1.024 (1.011–1.038; 0.000462) | Positive point estimates in all four results. |
| rs35424628 | 0.9783 (0.116) | 0.982 (0.956–1.008; 0.177) | Negative point estimates in all four results, but the insomnia associations are imprecise. |

FinnGen odds ratios and intervals are exponentiated from the G-aligned log-odds beta and its standard error; the MVP OR is as rounded in the source. The MVP file reports imputation r² of 0.9576 for rs8129315 and 0.6812 for rs35424628; assess the original study's QC and intended thresholds before relying on the latter. These nominal, post-selected comparisons do not establish trait specificity, a significant difference in effects, an absence of insomnia association, or a shared causal variant. The [aligned values](adarb1_nearby_markers_aligned_four_sources.csv) retain estimates and standard errors; source input checksums and release locations are documented in the linked reports.

### Rare insomnia minimum-p markers: limited reference resolution

The MVP minimum-p row for rs188632112 reports effect-allele frequency **0.9988** and imputation r² **0.652**. The FinnGen minimum-p row for rs1320256433 reports alternate-allele frequency **0.00211**. These are rare contrasting alleles, so both minimum-p rows need particular source-specific imputation, allele-count and association-quality review before being prioritized; no universal quality cutoff is applied here.

After coordinate conversion and exact reference-allele checks, the MVP marker's rare G allele appeared on only **2 of 1,006** EUR-labeled reference haplotypes. Its descriptive r² with rs2838788 was 0.0013, but two rare copies are far too few to support a stable LD interpretation. The FinnGen marker had no matching biallelic rsID/position/allele record in this reference at the mapped coordinate, so its pairwise LD was **not calculated**, not assumed to be zero. Neither observation rules out association, proves an independent signal, or measures LD in the source cohorts. The [receipt](adarb1_insomnia_regional_ld.json) preserves source rows, mapping, counts and results; script: `work/adarb1_insomnia_pair_ld.py`. Only the position was lifted over, with a unique positive-strand mapping; source-to-reference allele equivalence still requires study-specific verification.

### Exploratory frequency and quality sensitivity

To examine sensitivity to the rare minimum-p rows, I repeated the p-value ranking after requiring reported minor-allele frequency ≥1%. For MVP only, I also required the source's reported imputation r² ≥0.8. These cutoffs were selected **after** inspecting the region; they are a sensitivity exercise, not prespecified QC or proof that the remaining results are reliable. FinnGen's scanned result fields did not provide an equivalent imputation-quality statistic, so the filtering is not equivalent across studies.

| Source and filter | Passing / window rows | Smallest-p passing row | p | rs2838788 passes? |
|---|---:|---|---:|---|
| MVP: MAF ≥1%, r² ≥0.8 | 1,649 / 4,880 | rs8134789 | 0.0002604 | Yes: p = 0.001884, reported r² = 0.9004 |
| FinnGen: MAF ≥1% only | 2,099 / 4,772 | rs11701589 | 1.53 × 10⁻⁶ | Yes: p = 0.000345 |

The minimum-p marker changes in both datasets, and neither filtered minimum reaches 5 × 10⁻⁸. This supports caution about ranking purely by marginal p-value; it does not establish a shared or independent causal signal. The candidate's survival of these filters is not independent replication. This comparison does not harmonize variants across files, address LD or apply source-specific full QC. The [MVP receipt](adarb1_mvp_frequency_quality_sensitivity.json) and [FinnGen receipt](adarb1_finngen_frequency_quality_sensitivity.json) preserve the top ten passing rows, fields and counts. Script: `work/adarb1_insomnia_quality_sensitivity.py`.

In a subsequent allele-specific reference check, the filtered MVP minimum-p marker **rs8134789 G** and rs2838788 T had phased-haplotype **r² = 0.695** (r = −0.834) in 503 EUR-labeled 1000 Genomes samples. The filtered FinnGen minimum-p marker **rs11701589 T** and rs2838788 T had **r² = 0.025** (r = −0.157) in the same panel. Each tested SNP had a unique positive-strand GRCh38-to-GRCh37 mapping, a matching biallelic reference record and 503 phased calls. The signs describe co-occurrence of the specified allele indicators, not the direction of a disease effect.

These results show substantial reference LD with the filtered MVP marker and low pairwise LD with the filtered FinnGen marker, but they do not prove that either association is the same or independent after conditioning. The pooled general-European panel is not a Finnish-specific reference, and source-specific imputation and allele identity still need review. [Reference-LD receipt](adarb1_filtered_peak_reference_ld.json); script: `work/adarb1_filtered_peaks_ld.py`.

### Do the filtered insomnia markers also associate with frailty?

I checked the four raw source rows at each marker's mapped coordinate and aligned the same allele pair. Both Mak rows and the MVP row use the allele shown below as their effect allele; FinnGen's reported alternate-allele beta was reversed. These are post-selected marginal lookups, and the estimated effects for FI, physical frailty and insomnia are on different scales.

| Marker / aligned allele | Mak FI p | Mak physical p | MVP insomnia OR (p) | FinnGen insomnia OR, approx. 95% CI (p) |
|---|---:|---:|---:|---:|
| rs8134789 G / A | 0.0165 | 0.00228 | 0.9796 (0.000260) | 0.980 (0.966–0.993; 0.00362) |
| rs11701589 C / T | 0.644 | 0.726 | 1.020 (0.0289) | 1.052 (1.030–1.073; 1.53 × 10⁻⁶) |

For rs8134789 the estimated direction on G is negative in all four files, and each unadjusted p is below 0.05. It is correlated with rs2838788 in the reference panel and is not independent confirmation. At rs11701589 both Mak frailty estimates are positive but imprecise; the nominal p-values do not establish that a frailty effect is absent or differs significantly from the insomnia effect. Stronger insomnia evidence at this marker does not establish a shared FI/physical-frailty signal. Each result is affected by selection, LD, different trait scales and power, phenotype overlap, potential sample overlap and allele/source QC. Neither row establishes causality.

The [aligned two-marker table](adarb1_filtered_peaks_aligned_four_sources.csv) records betas, SEs and source p-values. FinnGen OR intervals are Wald intervals from the aligned beta; MVP ORs are as rounded in the source. For these two rows the MVP reported imputation r² values are 0.9973 and 0.9273, respectively. Position, allele and reference checks do not settle every source-specific coding or QC issue.

Mak coordinates are treated as GRCh37 under the existing source-matching assumption; the two clinical insomnia files use GRCh38. The ±250 kb intervals were centered separately on the matched candidate coordinates in each build and are not an exact lift-over of one common interval. Distances are within the respective source build. The four source populations and outcome definitions differ, the two Mak outcomes share UK Biobank, and participant overlap with other analyses may exist. The scanned window is a convenient local view, not a biologically defined locus boundary.

The original candidate was shortlisted from prior insomnia–FI analyses and post-selected follow-up; both its nominal evidence and the region's small p-values should be interpreted accordingly. Neither this scan nor the [mouse ADARB1 circadian background](shared_candidate_gene_background.md) demonstrates a mechanism in human insomnia or frailty. The next locus analysis needs source allele/build verification, comparable QC, dense ancestry-appropriate LD and sufficiently informative fine-mapping before any colocalization claim.

## Data and checks

The [regional rows](adar_region_marginal.csv) retain source, build, coordinate, variant ID, p-value and distance; the [receipt](adar_region_marginal.json) records scan counts and source-specific five smallest p-values. We checked that rs2838788 appeared at the expected position in each source, its p-values match the existing candidate matrix, the CSV row totals match the receipt, and its SHA-256 matches the receipt. Source files were previously checksum-verified as documented in [the Mak report](mak2023_same_cohort_fi_fp_ldsc_review.md) and [the candidate follow-up](song124_mvp_hfrs_followup.md); the FinnGen insomnia source was checked in [the six-source follow-up](22_candidate_six_trait_followup.md). The script is `work/adarfocus_region_scan.py`. No individual data were analyzed.
