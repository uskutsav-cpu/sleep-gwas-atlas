# Public muscle expression and sarcopenia resource: GSE167186

Checked 29 October 2026. [Perez et al. (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9792217/) studied lower-limb skeletal muscle from young, older non-sarcopenic and older sarcopenic men. The authors use “sarcopenia” and “frail” interchangeably for their functionally impaired group; this is **not necessarily the project's FI or Fried frailty score**. They report bulk RNA sequencing on 72 men and a separate 17-biopsy single-nucleus analysis of six young and 11 older samples. Sample selection and participant characteristics make this a limited muscle-context source, not a general-population replication.

The [GEO series GSE167186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167186) publishes a small [bulk counts file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/GSE167186_counts.csv.gz) and [series metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/matrix/GSE167186-GPL20301_series_matrix.txt.gz). The released bulk matrix has 76 named sample columns: 24 labeled sarcopenia, 29 old healthy, 19 young healthy, and four unclassified. The 72 classified labels match the paper's reported group totals, but that observation alone does not prove the exact paper analysis sample or why four columns are unclassified. We do not treat unclassified samples as controls. Public per-sample metadata include physical-performance measures; we did not use them in an association analysis.

## Candidate expression feasibility

In a descriptive, postselected screen of four previously discussed gene symbols, all four had a unique matching row. Across the 53 samples labeled old healthy or sarcopenia, FOXP2, SMAD5 and RBM6 each had at least one raw read count in every sample; MEIS1 had at least one count in 52/53. With the arbitrary stricter threshold of ≥10 counts, MEIS1 was detected in 10/53 older samples. [Per-group detection table](stage6_gse167186_bulk_candidate_detection.csv); [input hashes and sample counts](stage6_gse167186_bulk_candidate_detection.json); reproduce with `python3 work/gse167186_candidate_detection.py` after downloading the two GEO files above. The table also reports young and unclassified groups separately.

These are **unnormalized bulk counts**: they reflect cell mixture, sequencing depth and assay properties. The thresholds are only a screening convention and cannot establish true expression absence, differences between groups, specific expressing cell types, eQTL effects or the biological role of these genes. The detection screen itself performed no differential-expression test; the later exploratory analyses are documented below. The seven cell clusters reported in the separate single-nucleus experiment should not be attributed to these bulk detections. Gene-symbol matches do not make nearby genes causal.

This is a useful, openly accessible molecular and phenotype resource for a possible prespecified muscle follow-up, but the files inspected do not establish linked participant genotypes at the project variants or insomnia outcomes. They cannot validate a genetic score, establish that genetic risk acts through muscle, or show that sleep disturbance causes sarcopenia. A proper expression study would verify sample inclusion and processing, use appropriate normalization and covariates, account for limited sample size and gene selection, and distinguish cell-type composition from within-cell expression. Any future linkage of data must be authorized and documented.

## Exploratory age-adjusted comparison of four nominated genes

We compared the 24 samples labeled `Sarcopenia` with the 29 labeled `Old Healthy`; young and unclassified samples were excluded. For each of the four genes identified in the preceding detection screen, an ordinary linear model estimated the group difference in `log2(counts per million + 0.5)` after adjustment for continuous age. CPM used each sample's sum across all 36,694 matrix rows as the library size. These are **simple library-size-normalized** counts, not TMM, DESeq2 normalization or a full RNA-seq differential-expression workflow. Median ages were 77 years (sarcopenia; range 67–92) and 71 years (old healthy; range 64–87); age adjustment is linear and cannot rule out residual confounding. The paper reports an all-male cohort; no sex-stratified analysis is possible from these groups.

| Gene | Sarcopenia − old healthy coefficient on log2(CPM + 0.5), age-adjusted (95% CI) | Nominal p | BH-adjusted p across these four genes |
|---|---:|---:|---:|
| FOXP2 | +0.387 (0.083 to 0.690) | 0.0135 | 0.0539 |
| SMAD5 | −0.013 (−0.141 to 0.114) | 0.836 | 0.836 |
| RBM6 | +0.194 (0.003 to 0.384) | 0.0461 | 0.0923 |
| MEIS1 | +0.316 (−0.106 to 0.737) | 0.139 | 0.185 |

None passes BH 0.05 across these four tests. FOXP2 and RBM6 were nominally positive, but the genes were chosen from earlier project analyses, and the four-test adjustment does not correct the wider research search. Intervals assume the usual linear-model error assumptions; sequencing composition, library preparation or batch, low counts (especially MEIS1), unequal variance, sample choice, and cell composition can affect results. Because of the pseudocount, coefficients should not be interpreted as exact gene-expression fold changes. A leave-one-sample-out descriptive check retained positive coefficient signs for FOXP2, RBM6 and MEIS1; this is not a substitute for an independent replication or uncertainty analysis. No conclusions about no effect should be drawn for nonsignificant genes.

This comparison is a project-specific exploratory analysis, not the original paper's test, a test of insomnia, a cell-type comparison, a variant-to-gene association, or validation of either score. Its selected sarcopenia groups are not the project FI or FFS, and direct linkage of the 76 metadata columns to the precise 72-sample published analysis remains unverified. Future work should establish technical covariates and sample identity, use a standard count-model workflow with appropriate filtering and normalization, and seek independent validation before biological claims.

[Results](stage6_gse167186_candidate_age_adjusted.csv); code `Rscript work/gse167186_candidate_comparison.R`; R 4.6.1. Input SHA-256: counts `50e42a02981dc991845e3ee124eb653a64557811bc369b1e691799cf030f5a67`, metadata `d692c5e3b4210200493745d026efa93413289bbdfe9a46be0072923dae9ef969`. The code checks sample ordering, group sizes, unique target-gene rows and nonnegative counts.

## Count-model follow-up

We subsequently fit an [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) quasi-likelihood negative-binomial model on the same 53 labeled older men, with sarcopenia group and linear age as covariates. We applied `filterByExpr`, TMM library-size normalization, robust dispersion estimation and a robust quasi-likelihood fit. Of 36,694 input count rows, 8,304 passed filtering; all four previously nominated target genes were retained. A BH correction was applied **only across these four selected genes**.

| Gene | Log2 fold change, sarcopenia versus old healthy | Quasi-likelihood p | BH p across four genes |
|---|---:|---:|---:|
| FOXP2 | +0.402 | 0.00730 | 0.0292 |
| SMAD5 | −0.0346 | 0.495 | 0.495 |
| RBM6 | +0.179 | 0.0467 | 0.0933 |
| MEIS1 | +0.264 | 0.195 | 0.260 |

Under this count-based model, FOXP2 meets the 5% BH threshold for this four-gene family, whereas it did not in the previous simple log-CPM model (BH p=0.0539). The four genes were selected after prior project results and the correction does not account for the wider search; the choice of follow-up and model was exploratory, and results depend on normalization and error assumptions. Age differences, small sample size, unmeasured batch and other covariates, possible differences in muscle-cell composition, and uncertain mapping to the paper's exact analysis sample remain. The edgeR result is a **bulk expression association**, not evidence that FOXP2 is a causal gene or that our variants regulate it in muscle. No insomnia, genotype–expression link, individual genetic score or external replication was tested; neither FI nor FFS was measured by this definition. Nonsignificant results do not establish zero effects.

[Machine-readable results](stage6_gse167186_edger_candidate_results.csv); run `Rscript work/gse167186_edger_comparison.R` with `work/rna_Rlib` populated. R 4.6.1, edgeR 4.10.7, limma 3.68.5, locfit 1.5.9.12 and statmod 1.5.2; input hashes are listed in the preceding section. Script SHA-256 `8e3909b21223c158988a6359a658c93da1cd7c015c9db6fd27854e1772501137`; output SHA-256 `599afde0762958623d50f456ac1e3dc384317826183a7e2ffbd119976f46fc3c`. Public sample identity or covariate verification and a genuinely independent replication remain necessary before a mechanistic claim.
