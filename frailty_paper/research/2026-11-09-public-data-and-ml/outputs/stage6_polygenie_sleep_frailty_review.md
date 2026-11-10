# Published frailty and insomnia scores in the GCAT cohort

[PolyGenie](https://polygenie.igtp.cat/) is a public results interface for a [2026 analysis](https://doi.org/10.1093/nargab/lqag056) of 135 genome-wide polygenic scores and 1,483 clinical, questionnaire and other traits in the Spanish GCAT cohort. Its paper describes 4,988 genotyped participants of predominantly Iberian, White Western European ancestry, aged 40–65 at cohort recruitment, with individual-level data requiring controlled access. Exact target-specific numbers and overlap or relatedness with every project cohort are not established by the public table. The cohort is a distinct named target, not proof of full individual independence.

The tool identifies its Fried-frailty score as derived from the Ye et al. UK Biobank FFS GWAS (386,565 participants), and its insomnia score as derived from Jansen et al. UK Biobank insomnia summary statistics (386,533 participants). These broad scores are different from this project's 15-marker weighted prototypes and 22-marker count; results for them cannot validate our scores. Because the sources overlap the GWAS used in our project, these are external target-cohort PheWAS observations, not independent discovery of the genetic effects.

We used the interface's comparison of the top score decile against the bottom decile, separately for the two scores and ICD diagnosis outcomes. The paper describes logistic regression adjusted for age, sex and ten genetic principal components. The following are displayed odds ratios and nominal p-values, not per-standard-deviation effects:

| Existing broad score | GCAT ICD outcome | Top vs bottom decile OR | Nominal p |
|---|---|---:|---:|
| Fried frailty | G47: sleep disorders | 1.528 | 0.0360 |
| Fried frailty | F51: sleep disorders not due to a substance or known physiological condition | 1.141 | 0.490 |
| Fried frailty | R53: malaise and fatigue | 1.739 | 0.00151 |
| Insomnia | G47: sleep disorders | 0.986 | 0.945 |
| Insomnia | F51: sleep disorders not due to a substance or known physiological condition | 1.598 | 0.0258 |
| Insomnia | R53: malaise and fatigue | 1.004 | 0.981 |

These checks were chosen after seeing the available results. The nominal p-values are uncorrected for the many outcomes, scores and analyses examined; they do not establish a robust shared effect. G47 and F51 are broad clinical codes, not identical definitions of insomnia; R53 is not FI or FFS and overlaps conceptually with exhaustion. Different nominal significance between the two scores or phenotypes is not a test that their effects differ. These cross-sectional score–diagnosis associations do not establish causation, discrimination, calibration or incremental predictive utility.

The public browser also has Phecode and questionnaire results, exported alongside the ICD rows. Interpretation of questionnaires requires units, case definitions, denominators and missingness, which were not verified. The underlying individual genotypes and outcomes cannot be recovered from this aggregate table. The interface requests citation of PolyGenie and the original GWAS for publication use.

A [sleep-duration follow-up](stage6_polygenie_sleep_duration_followup.md) summarizes six displayed questionnaire estimates for the two published scores. None has nominal p below 0.05, but this is not evidence of no association; exact units, sample sizes and confidence intervals remain unverified, and these results do not test our prototypes.

Reproducibility: the interface's documented read-only Dash result callback was queried on 29 October 2026 UTC for `Fried Frailty` and `Insomnia`, 10 groups, bottom reference, and the ICD/Phecode/questionnaire tabs. The [10 extracted ICD and Phecode results](stage6_polygenie_sleep_frailty_results.csv), [receipt and hashes](stage6_polygenie_sleep_frailty_receipt.json), `work/stage6_polygenie_sleep_frailty_extract.py`, source metadata `work/polygenie_gwas.json` and target metadata `work/polygenie_target_inventory.json` preserve the inspection. This is exploratory and not prespecified replication.

**Multiple-testing display:** an [audit](stage6_polygenie_multiple_testing_audit.md) found that the live “domain-wise” threshold corresponds to dividing by the number of domains, rather than the number of outcomes within each domain. Under within-domain Bonferroni corrections based on the returned rows, neither G47 nor R53 passes its domain threshold. The live global FDR setting does flag R53, but that is not prespecified validation of our scores.

A [separate PheWeb outcome inventory](stage6_gcat_pheweb_endpoint_audit.md) returned ten public SNP-GWAS phenotypes and no insomnia or frailty label. The sleep outcomes shown by PolyGenie cannot therefore be assumed to have full GWAS results in that PheWeb.
