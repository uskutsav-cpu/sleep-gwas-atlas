# Reference-panel comparison of the 22-marker count

The 22-marker equal-weight research count was calculated using the exact scored alleles in two 1000 Genomes phase 3 reference groups: 504 EAS-labeled and 503 EUR-labeled individuals. All 22 variants matched one biallelic hg19 record by position, rsID and allele pair, and all individuals had called diploid genotypes at each site. These are pooled reference groups, not representative samples of all people in either ancestry group and not the TPMI target cohort.

| Reference group | Full 22-marker count mean | Full count SD | Markers with sample MAF ≤ 1% | Correlation of full count with the 17-marker subset |
|---|---:|---:|---:|---:|
| EAS | 27.56 | 2.36 | 5 | 0.9991 |
| EUR | 23.84 | 3.04 | 0 | 0.9443 |

The 17-marker subset is the set for which TPMI displayed rounded insomnia associations in the [browser screen](stage6_tpmi_shared22_insomnia_review.md); this subset is not a separately developed or validated score. In this small EAS reference sample the five omitted markers contributed nearly constant allele counts, so the full and subset sums were highly correlated. That observation concerns only genotype counts in this reference sample. It does not prove that the subsets have equivalent associations, prediction, calibration or clinical usefulness in TPMI. Target genotype and LD data, effect-allele confirmation, and full-precision outcome associations remain necessary. The different group means are allele-count summaries and must not be interpreted as between-population differences in disease risk.

The [per-marker table](stage6_shared22_cross_population_frequencies.csv) and [summary receipt](stage6_shared22_cross_population_reference.json) are produced by `python3 outputs/stage6_shared22_cross_population_reference.py`. The EUR SD agrees with the previously reported count variance (3.0369² ≈ 9.223), an internal consistency check. The referenced VCF was accessed through indexed remote regions without full-file checksum verification; the analysis includes no phenotype information. Differences in allele frequencies, marker coding, LD, ascertainment and discovery ancestry all require attention in any transfer study.
