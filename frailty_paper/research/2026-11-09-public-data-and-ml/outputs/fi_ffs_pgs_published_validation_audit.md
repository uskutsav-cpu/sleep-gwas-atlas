# FI- and FFS-derived score validation: primary supplement versus catalog

**General-factor score follow-up:** The same source discrepancy also affects PGS005221 (GF-PRS). For FI, the PGS Catalog reports ELSA/PISA/LBC R² of 0.009/0.006/0.0002, while the paper's Supplementary Table 44 reports 0.01/0.02/0.01 (the LBC row is Wave 1). ELSA is consistent with rounding, but the PISA and LBC pairs are not. The reason is unresolved; see the [general-factor comparator review](stage6_frailty_general_pgs_comparator.md). These are all results for published scores rather than this project's prototypes.

## Physical frailty outcome

[Foote et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12626885/) evaluated both aggregate FI-derived and Fried-phenotype-derived polygenic scores against observed frailty phenotypes. In their [Supplementary Table 44](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM4_ESM.xlsx), the LBC1936 Wave 1 **frail versus non-frail** multinomial model reports the following odds ratios per standard deviation of score:

| Predictor | OR | 95% CI | Reported FDR q |
|---|---:|---:|---:|
| FI-derived PRS | 1.47 | 1.17–1.86 | 0.02565 |
| Fried-phenotype-derived PRS | 1.46 | 1.16–1.85 | 0.02640 |

These are associations with a categorical physical phenotype in an external cohort, not absolute risks, a validated joint insomnia–frailty score, or proof the two scores perform equivalently. The similar numbers have not been subjected here to a paired comparison. The paper's supplementary note describes its categories and multinomial regression; participants with a pre-frail classification are a separate category. [Supplementary note](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM1_ESM.pdf). Sample composition and exact analyzable counts should be checked before re-use; the note's phenotype category counts sum to more than the 1,005 stated in the main paper for the LBC1936 cohort, so they should not be treated as confirmed regression sample counts.

## An unresolved discrepancy in FI-outcome statistics

My earlier briefing quoted values from the [PGS Catalog publication entry](https://www.pgscatalog.org/publication/PGP000735/). The published Table 44 now provides a second, primary-source set with the same predictor labels and cohort names, but several values **do not match even allowing for reported rounding**:

| Cohort, observed FI | Source | FI-derived score R² | FFS-derived score R² |
|---|---|---:|---:|
| ELSA | PGS Catalog | 0.020 | 0.0138 |
| ELSA | Paper Table 44 | 0.02 | 0.01 |
| PISA | PGS Catalog | 0.021 | 0.030 |
| PISA | Paper Table 44 | 0.02 | 0.02 |
| LBC1936 | PGS Catalog | 0.005 | 0.019 |
| LBC1936 Wave 1 | Paper Table 44 | 0.02 | 0.03 |

Published point estimates and q-values also differ for these same labels. The reason is unknown: the records may encode different analysis specifications or contain an error, and the sources inspected do not resolve it. Do not present the catalog numbers as verified results from Table 44 or use them to claim a ranking of FI versus FFS performance. The R² entries are reported as incremental variance explained above covariates in the paper's methods; no comparison of absolute event risk follows. The table values are rounded as published, and equal rounded values do not demonstrate identical performance.

The [PGS Catalog API record PPM022628](https://www.pgscatalog.org/rest/performance/PPM022628) independently confirms that the catalog associates its ELSA metrics with PGS005228 and identifies the publication and covariates, but leaves phenotype details and performance comments blank. Its publication record lists catalog release on 24 July 2025 and journal publication on 4 August 2025. These dates do not establish that a version change caused the disagreement; no definitive source explanation was found.

## Additional sample-metadata discrepancy

The [FFS score’s catalog page](https://www.pgscatalog.org/score/PGS005229/) lists the LBC sample (PSS012045, n=1,005) with mean age 69.6 years and SD **6.98**. The paper’s [supplementary note](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM1_ESM.pdf), “Cohort Samples Used for the Main Polygenic Risk Score Analyses,” states LBC1936 Wave 1 mean 69.60 and SD **0.83**, also n=1,005. The same note reports PISA mean 60.34 and SD **6.98**, matching the PISA age SD shown in the catalog. These distinct published descriptions point to a sample-metadata inconsistency; the matching PISA SD raises a possible copy error but does not prove its cause. It does not itself explain the conflicting association or R² values above. Use the primary paper’s LBC age summary with the source explicitly identified, and seek clarification before relying on the catalog values.

A targeted search of Supplementary Tables 39–54 for the catalog’s distinctive rounded values (0.0138, 2.48e−6, 9.79e−7 and 3.93e−6; numeric match tolerance 0.01%) did not find those values. This limited search cannot rule out different reporting or an earlier analysis outside the workbook.

## Evaluation phenotype overlap

I checked all rows of Tables 52–54 for sleep-, insomnia-, tiredness-, fatigue-, energy- and exhaustion-related terms. PISA lists **seven** matching items: sleep problems, lethargy, sleep apnea, medicine for sleep or mood, longer sleep, poor sleep, and overall sleep quality. ELSA lists **one**: restless sleep during the past week. LBC1936 has **no explicit match** under those search terms; that does not establish that none of its broader health or psychological measures relates to sleep. These are documented item labels, not six insomnia diagnoses in PISA, and the mixed sleep-or-mood medicine item does not isolate insomnia. The PISA and ELSA FI evaluations therefore include sleep-related content. All matching source rows and their original descriptions are in `foote_fi_evaluation_sleep_items.csv`; outcome scales and coding should be checked for each evaluation. The observed outcomes here are distinct from the FI and FFS discovery GWAS phenotypes.

## Traceability

`foote_supplement_44_fi_ffs_score_results.csv` transcribes all 12 FI-PRS and FP-PRS rows from Table 44: FI outcomes in three cohorts and three LBC1936 physical-frailty contrasts. It retains the source row number, q-value and the outcome-specific metric. The publisher XLSX SHA-256 is `d170d2596ed43cb168d4463a722894969198aa6d8e493f73cdf7cd6d5637b5b4`; the publisher supplementary note PDF SHA-256 is `fc37c0f05285511520a19806049ee81658b5e962f39e07c498b75f901a32826f`. Neither source resolves the performance of our selected 15-marker prototype. The study uses European-ancestry scores and evaluations, so broader population transferability remains unproven.
