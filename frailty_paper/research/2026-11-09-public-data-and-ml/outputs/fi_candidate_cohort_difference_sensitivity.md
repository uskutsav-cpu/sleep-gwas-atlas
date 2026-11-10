# FI candidate effects: exploratory differences between cohorts

This is a **post-selection diagnostic**, not a replicated heterogeneity finding or a comparison of insomnia effects. Using the 15 existing candidates, I contrasted TwinGene FI and UK Biobank FI effects on the same listed allele. Under the simplifying assumption of independent cohort estimates, the estimated difference is `beta_TwinGene − beta_UKB`, with approximate standard error `sqrt(SE_TwinGene² + SE_UKB²)`; the 95% interval and two-sided p-value use a normal approximation. The full results are in `fi_candidate_cohort_difference_sensitivity.csv`.

| Variant; counted allele | UKB FI beta (SE) | TwinGene FI beta (SE) | Difference, TwinGene − UKB (95% CI) | Unadjusted p |
|---|---:|---:|---:|---:|
| rs12536335; G | −0.01889 (0.00345) | +0.03922 (0.01454) | +0.05811 (0.02882 to 0.08740) | 0.000101 |
| rs323509; A | +0.01931 (source SE in CSV) | −0.02735 (source SE in CSV) | −0.04665 (−0.07783 to −0.01548) | 0.00336 |

These are the two smallest p-values among 15 candidates; both happen to be below 0.05 without adjustment. **Neither is a valid post-selection discovery claim.** The candidate list was selected from earlier analyses that include the FI meta-analysis, which itself includes both cohorts. We have not verified cross-cohort relatedness, and selection, ancestry, allele representation, differing phenotype composition and sample size can affect this comparison. The UKB and TwinGene FI outcomes were separately quantile-normalized, so an equal standardized association is an imperfect target when the underlying indices differ. These contrasts cannot establish why estimates differ or isolate the influence of the insomnia item. For rs323509 the source rsID is multiallelic, adding unresolved allele-representation uncertainty.

The TwinGene beta is aligned using the project's independently audited sign adapter; the archive's original beta and the full allelic and source context remain in `candidate_twingene_fi_cohort_lookup.csv` and `candidate_twingene_sensitivity.md`. The source issue has not been confirmed by the original authors. The calculation assumes zero covariance; participant overlap or relatedness would alter its standard error. Values are not adjusted for the selection process or multiple comparisons; a nonsignificant difference is not evidence of equivalence.

A direct difference test is more informative than comparing one cohort's p-value to another's, but these assumptions and the small TwinGene sample prevent causal interpretation. This single-variant diagnostic does not relax the previously failed TwinGene LDSC heritability-precision gate, validate local genetic sharing, or turn the candidate list into a risk score.
