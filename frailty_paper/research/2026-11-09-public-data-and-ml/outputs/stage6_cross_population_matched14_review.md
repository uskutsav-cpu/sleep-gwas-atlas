# Matched 14-marker insomnia and knee-strength comparison

Fourteen of the 22 count-prototype variants have one position-and-allele-pair match in both the Taiwanese insomnia browser and the Japanese knee-strength GWAS. Holding that **same 14-marker set and counted alleles** fixed, the exploratory weighted-Z statistics are:

| Outcome | Weighted Z | Nominal two-sided p |
|---|---:|---:|
| Taiwan clinical insomnia | 2.27 | 0.023 |
| Japan lower knee strength | −0.38 | 0.70 |

The source knee-strength effects were sign-reversed so positive values correspond to lower strength. The insomnia association is nominal under the conditional analysis, whereas the lower-strength estimate is statistically uncertain. The two p-values are **not** a formal test of a difference between traits, and non-significance does not establish a zero effect or inconsistent biology. The 14 markers were selected based on availability after seeing the data; these and other screens are postselected, and p-values are not multiplicity-adjusted. The results do not establish a joint insomnia–frailty association or validate the 22-marker score.

Taiwan browser effects are rounded and interpreted under the assumption that the alternate allele is the effect allele. The Japanese file is unharmonized, although [source and reference checks](stage6_japanese_knee_allele_frequency_audit.md) support the tested allele pairs. Reference LD and genotype scaling differ across the analyses: the insomnia calculation uses 504 EAS reference individuals for both, whereas the knee-strength calculation uses 104 JPT individuals for LD and target study allele frequencies for scaling. Both set cross-chromosome correlations to zero. Target LD and target-sample/relative overlap with the discovery and selection studies remain unverified. Ancestry, clinical ascertainment, sample sizes and outcomes differ. Japanese knee strength is not FI, FFS or a clinical frailty diagnosis.

The [machine-readable receipt](stage6_cross_population_matched14.json) records the exact marker list and `python3 outputs/stage6_cross_population_matched14.py` reproduces the analysis using the project score and both source inventories. The prior full-available-marker screens are [Taiwan 17](stage6_tpmi_17_marker_exploratory_z_review.md) and [Japan 14](stage6_japanese_knee_weighted_z_review.md). One earlier execution failed due to a remote reference seek error; the successful run independently retrieved the 14 required regions. The reference VCF was not full-file checksum verified. Neither result gives individual discrimination, calibration or incremental utility.

TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x). The Japanese source is [Ito et al.](https://doi.org/10.1038/s42003-024-06108-6).
