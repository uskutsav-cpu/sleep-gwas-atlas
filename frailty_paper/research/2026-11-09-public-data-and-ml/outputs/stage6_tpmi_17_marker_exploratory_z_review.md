# Conditional 17-marker insomnia summary-statistic screen

For the 17 markers of the postselected 22-marker equal-weight count with displayed TPMI insomnia results, a conditional, approximate reference-LD-adjusted weighted-Z screen gave the following results:

| TPMI phenotype | Weighted Z | Nominal two-sided p | Leave-one-marker-out Z range |
|---|---:|---:|---:|
| Insomnia, 327.4 | 2.44 | 0.0146 | 2.04 to 2.68 |
| Organic or persistent insomnia, 327.41 | 1.28 | 0.201 | 1.02 to 1.55 |

**These are exploratory sensitivity statistics, not replication or a measured score effect.** The 17-marker subset was defined by available browser results after the candidate set and outcomes were inspected. The two phenotypes overlap. The nominal p-values are unadjusted for multiple analyses and cannot be assumed well calibrated given the approximations. The broader insomnia result is nominally positive under the stated assumptions; the more specific result is statistically uncertain, which is not evidence of no effect.

The calculation assumes that browser beta is an ALT-allele effect and reverses signs where the score counts REF. The [published-table comparison](stage6_tpmi_published_allele_comparison.md) supports that assumption for selected records, including an insomnia record, but does not independently verify the 17 project variants or all multiallelic contrasts. The analysis uses manually transcribed, rounded effects and SEs; `z=beta/SE` may differ from a more precise source calculation. Each count weight is 1 and `a` is the corresponding counted-allele dosage SD in 504 1000 Genomes phase 3 EAS reference individuals. Using the reference dosage-correlation matrix `R`, with cross-chromosome correlations set to zero, the statistic is `aᵀz / sqrt(aᵀRa)`. This relies on compatible target samples, outcome models and allele orientations, and on reference LD adequately approximating the target Z-score covariance. TPMI-specific LD and full-precision summary statistics were not used; reference source VCF files were not independently full-file checksum verified.

The leave-one-out range reflects only changes in the selected marker set with the same input and reference assumptions. It is neither a confidence interval nor evidence that removing any marker improves prediction. The [reference comparison](stage6_shared22_cross_population_reference_review.md) showed that the available and full-count genotypes were highly correlated in this EAS reference sample, but this does not establish equivalence of association or prediction in TPMI. Five score markers remain unverified, as does exact participant/relatedness overlap between discovery, selection and target datasets. Neither phenotype measures FI or FFS, and the test cannot establish calibration, discrimination, incremental benefit, or a joint sleep–frailty outcome.

A subsequent [reference-subset sensitivity](stage6_tpmi_reference_population_sensitivity.md) returned Z values of 2.36–2.51 for 327.4 across selected EAS subsets, with TPMI associations held fixed. Those subsets overlap the original reference and do not constitute independent replication or capture target-study uncertainty.

The [receipt](stage6_tpmi_17_marker_exploratory_z.json), [component table](stage6_tpmi_17_marker_exploratory_components.csv), source inventories and `python3 outputs/stage6_tpmi_17_marker_exploratory_z.py` record the computation and inputs. TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).
