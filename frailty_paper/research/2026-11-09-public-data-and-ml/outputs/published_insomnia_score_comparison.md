# Benchmarking the 15 candidates against a published insomnia score

## An existing sparse score

[PGS003859](https://www.pgscatalog.org/score/PGS003859/) is a published 129-variant insomnia genetic score. Gibson and colleagues selected genome-wide significant variants from a 23andMe insomnia-symptom GWAS, used European 1000 Genomes LD clumping, and evaluated the main score in UK Biobank. The discovery phenotype and evaluation phenotype were defined differently. The score was constructed for insomnia, not for frailty or joint insomnia–frailty risk. Its original catalog file reports effect weights as log odds ratios and does not report a genome build or variant coordinates. [Gibson et al.](https://link.springer.com/article/10.1186/s12916-023-02832-8)

In 336,975 unrelated white-British UK Biobank participants, the authors reported OR **1.08 (95% CI 1.07–1.09)** per one standard deviation of this score for usually versus never/rarely/sometimes having trouble falling or staying asleep. This is a score–insomnia association in that cohort, not an absolute risk, a frailty estimate, or proof of usefulness in other ancestries. The score was originally used as an instrument in a phenome-wide MR investigation; an MR instrument is not automatically an optimized predictor. The authors discuss overfitting concerns for separate sensitivity scores selected using UK Biobank; the main 129-variant score used 23andMe discovery, so one should not conflate those constructions.

## Comparison with the project's candidates

I compared exact rsID strings in the complete 129-row PGS scoring file with the 15-row project prototype and the 124 published Song candidate rows:

| Comparison | Exact rsID overlap |
|---|---:|
| PGS003859 versus the 15-marker prototype | 1 of 15 |
| PGS003859 versus the full 124 Song candidates | 1 of 124 |

The shared ID is **rs113851554**, with T as effect/counted allele and G as the other allele in the score file and prototype. PGS003859 gives T a log-odds weight of approximately **0.09011**; the project insomnia input's T-aligned estimate is approximately **0.10975**. Their different discovery phenotypes, datasets, samples and selection procedures make the weights non-interchangeable. In the project inputs this candidate's reported FI p=0.1987 and FFS p=0.19, so the strong sleep evidence does not by itself establish a frailty effect; these p-values also do not prove no frailty effect.

This is a literal ID-and-allele-string lookup, not a complete locus or LD comparison. The catalog file has no position or reported build; this rsID has multiple alleles in the reference, and the precise allelic representations require verification before harmonization or scoring. Different correlated variants can tag the same locus. Thus 1/15 does **not** mean the scores share only one biological signal, nor does absence mean a candidate has zero insomnia effect. Neither candidate list is an independent test set.

Reproducibility: source [original PGS003859 scoring file](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS003859/ScoringFiles/PGS003859.txt.gz), downloaded 27 October 2026; SHA-256 `c22ffe5c06eb8a6957d7551718cfed445ab1a3c3c14b2acd477e0332e5814e30`. The file contains 129 data rows with 129 unique rsIDs. The lookup uses exact strings only and retains nonmatches; no scoring was performed. `insomnia_score_candidate_id_overlap.csv` records the complete lookup against the prototype.

## A testable shared-risk model

Define the outcome before choosing predictors. A practical first target is **incident frailty by a fixed follow-up horizon** among participants non-frail at baseline; evaluate FI and Fried definitions as separately specified outcomes in the same cohort. If the actual research question is the joint onset of insomnia and frailty, that needs a distinct endpoint with sufficient cases and a clear temporal definition.

1. Assemble independent genotyped development and held-out validation cohorts with both sleep and frailty measurements. Audit overlap and relatedness with *all* discovery GWAS, not only overlap between the held-out and development cohorts. If UK Biobank contributes to a discovery GWAS, a random split of UK Biobank evaluation participants is insufficient unless discovery excludes the test individuals and close relatives. Match and verify allele/build/imputation QC before scoring.
2. Freeze an age/sex/ancestry-PC baseline, then compare (a) baseline + FI score, (b) baseline + FFS score, (c) baseline + insomnia score, (d) baseline + all three with coefficients estimated in training, and (e) any separately developed sparse 15–30-marker candidate. Evaluate adding insomnia score to each frailty-score baseline and vice versa. Treat the current prototype as exploratory; selecting markers or tuning weights on validation data invalidates the test. Study clinical sleep measures in a separate prespecified model, because predicting using measured insomnia answers a different question from a genetic-only model.
3. Use training-only standardization and model fitting. For fixed-horizon incident frailty, account for loss to follow-up and death before frailty using an appropriate competing-risk or censoring-aware method; prespecify the estimator and how frailty is ascertained. Report event counts, baseline risk, calibration-in-the-large, calibration slope, observed versus predicted risk by group, a time-appropriate Brier score, discrimination with uncertainty, and paired incremental performance compared with the baseline. Include subgroup evaluation where sample size supports it and external validation across relevant populations.
4. Reconstruct FI without insomnia, and then without insomnia and tiredness, keeping the denominator and missing-item rules explicit. Compare these with the original FI and FFS within the same people; a changed FI has a changed interpretation. If the desired quantity is genetic correlation rather than score prediction, each alternate FI requires its own valid GWAS summary statistics followed by carefully matched LDSC.

No cohort genotypes or outcomes have been used for this comparison, and no score has been independently fitted or validated in this project. A fixed number of markers cannot substitute for an endpoint, representative data or demonstrated performance.

## A second warning for multivariate frailty scores

[Foote et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12626885/) used insomnia as one of the 30 deficits in their genetic factor model; it loaded on the general frailty factor. Their multi-factor score is an important published frailty comparator, and their paper reports evaluation for FI across three cohorts and the physical phenotype in LBC1936. However, an association of an insomnia-containing factor with insomnia is not an independent demonstration of a newly shared sleep–frailty mechanism. A leave-insomnia-out factor would be a sensitivity analysis requiring the model and GWAS to be refit, not a simple subtraction from the published score. The published factor scores also differ from the separate FI and FFS PGS files audited in this project.
