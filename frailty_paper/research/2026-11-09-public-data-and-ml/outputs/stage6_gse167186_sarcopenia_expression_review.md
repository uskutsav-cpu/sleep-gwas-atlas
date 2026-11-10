# Open skeletal-muscle sarcopenia expression follow-up

[GSE167186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167186) releases bulk muscle RNA-seq and a separate single-nucleus dataset. We analyzed the **bulk** counts for 24 participants labeled sarcopenia and 29 labeled old healthy; 19 young healthy and four unclassified samples were excluded. GEO metadata list per-sample age and physical measures, but no sex field was found. This endpoint is sarcopenia, not necessarily FI or the Fried frailty phenotype, and the dataset does not establish insomnia or linked genotypes for score testing. The 17 single-nucleus samples were not used and cannot be counted as independent bulk samples.

We selected the same four existing candidate genes as the prior two muscle checks and compared log₂ counts per million using OLS with age adjustment and HC3 standard errors. Sex, cell composition and technical batch were not adjusted; the submitted metadata did not supply sex. BH values correct only the four postselected genes, not the broader project search. A separate negative-binomial RNA-seq analysis and detailed covariate audit would be needed for stronger inference.

| Gene | Age-adjusted sarcopenia minus healthy difference in log₂ CPM | HC3 SE | Nominal p | BH q within four genes |
|---|---:|---:|---:|---:|
| FOXP2 | +0.403 | 0.160 | 0.015 | 0.061 |
| SMAD5 | −0.013 | 0.069 | 0.848 | 0.848 |
| RBM6 | +0.197 | 0.097 | 0.048 | 0.097 |
| MEIS1 | +0.366 | 0.208 | 0.085 | 0.113 |

None passes BH 0.05 in this small post hoc panel. FOXP2 and RBM6 nominal p-values are not sufficient for a claim of replication, and the direction of the RBM6 estimate differs from the previous confounded microarray comparison. MEIS1 has a low median raw count (6); FOXP2's is 14, so low-count sensitivity is also relevant. Differences in study populations, definitions, normalization, platforms, cell composition and unmeasured sex limit comparison with other sources; participant overlap has not been verified. Null results do not exclude effects, and bulk gene expression does not establish a genetic variant–gene mechanism.

The series has 76 bulk metadata columns: the four unclassified samples explain the distinction from the 72 classified young/old/sarcopenia samples often reported. Reproduce with `work/geo_muscle/analyze_gse167186_candidates.py`; results are `outputs/stage6_gse167186_candidate_expression.json`. Files: `work/geo_muscle/GSE167186-GPL20301.txt.gz` and `work/geo_muscle/GSE167186_counts.csv.gz`. The single-nucleus GEO series metadata were inspected only to confirm it is a separate 17-sample modality.
