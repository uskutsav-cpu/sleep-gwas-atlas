# FI versus FFS weights: conditional Taiwanese insomnia screen

This postselected sensitivity combines the **12 available markers** of the 15-marker prototype with each of its three existing weight sets, using rounded TPMI browser estimates and 1000 Genomes East Asian reference LD. It does not test the complete 15-marker scores or compare prediction in people.

| Weight set | Insomnia 327.4 weighted Z | Organic/persistent insomnia 327.41 weighted Z | Full vs 12-marker score correlation in EAS reference |
|---|---:|---:|---:|
| FI | 2.79 | 1.04 | 0.970 |
| FFS | 2.18 | 1.31 | 0.936 |
| Insomnia | 2.92 | 1.29 | 0.980 |

These statistics are conditional on browser beta being an ALT-allele effect. The available [source comparison](stage6_tpmi_published_allele_comparison.md) supports that convention in examples but does not verify these particular markers. The FFS weight directions were inferred and remain unresolved. All tests are post hoc; target overlap, complete target alleles, and score selection independence have not been established. Their nominal p-values are in the [receipt](stage6_tpmi_12_weighted_z.json), but should not be interpreted as adjusted replication evidence. The two outcomes and the weight sets are correlated. Differences in their Z values do not demonstrate that one score predicts better, and the less precise 327.41 results do not establish zero effects.

For each included marker, the displayed beta/SE is oriented to the prototype's counted allele under the ALT assumption. With equal target samples and compatible models assumed, `a = weight × EAS-reference dosage SD`; we calculate `aᵀz / sqrt(aᵀRa)` using reference dosage correlations within chromosome and zero cross-chromosome correlation. Rounded values, reference LD, heterogeneous imputation, and any target-sample or model differences can distort the statistic. The 504 EAS-reference genotype checks matched all 15 GRCh37 biallelic allele pairs, but this is not target validation; the VCF was accessed by indexed regions without full-file checksum verification. Three target markers were not found in the browser, and multiallelic contrasts remain a concern.

The final column compares weighted genotype sums only in those 504 reference individuals. It does not establish similar associations, calibration or predictive performance in TPMI, and the three subsets correspond to different scores. FI, FFS and hospital-coded insomnia are different phenotypes; neither TPMI outcome is a direct FI or FFS measure. A valid FI-versus-FFS comparison still requires a sufficiently powered, independent cohort with linked genotypes and matched outcomes, predefined weights, resolved direction, and paired evaluation against the same baseline.

Run `python3 outputs/stage6_tpmi_12_weighted_z.py` to rebuild the receipt; it records input hashes and the source assumptions. TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).
