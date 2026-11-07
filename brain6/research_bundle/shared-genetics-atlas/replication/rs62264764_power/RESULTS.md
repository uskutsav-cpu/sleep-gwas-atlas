# rs62264764: illustrative ADHD replication power

A small registry cohort could miss the published chromosome 3 association even if its true effect is nonzero. Under one simplified scenario—**5,559 hypothetical ADHD cases, 100,000 controls**, the discovery effect and the 1000 Genomes EUR reference frequency—estimated power for the ADHD association is **69% at two-sided α=0.05** and **43% at α=0.05/6**. The six-test correction is illustrative, not a frozen replication family. These are planning scenarios, not estimates of the actual FinnGen GWAS power.

| Assumed true ADHD effect | Power at α=0.05 | Power at α=0.05/6 |
|---|---:|---:|
| Half the discovery estimate | 23% | 8% |
| Equal to the discovery estimate | 69% | 43% |
| 1.5 times the discovery estimate | 96% | 85% |

The input discovery estimate is |log odds ratio|=0.06470 for the A allele of rs62264764; the 503-person EUR reference A frequency is 0.1640. The published signal was selected from discovery data, so winner's curse and phenotype differences can make the actual replication effect smaller or different. The [FinnGen R13 F5_ADHD endpoint](https://r13.risteys.finngen.fi/endpoints/F5_ADHD) displays 5,559 endpoint individuals, but this has **not been verified as the GWAS case count**; 100,000 controls is an illustrative assumption, not a verified FinnGen control count. Frequency, case count, covariate information, imputation and effective sample size must be checked in an authorized GWAS release before making a study-specific power estimate. Registry ADHD and the discovery clinical phenotype may differ.

With the same effect, frequency and simplified model, reaching 80% power requires an effective sample size of about 27,348 at α=0.05, or 42,194 at α=0.05/6. With infinitely many controls these correspond to theoretical lower bounds of about 6,838 and 10,549 cases, respectively; finite controls require more. These numbers are not cohort recruitment targets without a validated design and parameter estimates. For fixed 5,559 cases, even the infinite-controls scenario yields 71% and 46% power at the two respective thresholds in this model.

The calculations use the normal approximation to an additive logistic-regression Wald test, genotype variance 2f(1−f), and N_eff=4 Ncase Ncontrol/(Ncase+Ncontrol). They assume unrelated samples, complete genotypes, homogeneous effects and no information loss from covariates. They do not address overlap with discovery cohorts, LD mismatch or proxy selection, phenotype definition, multiple independent signals, or uncertainty in the discovery effect. Actual performance may differ in either direction. Because a two-trait replication requires the ADHD test to pass as well as the sleep test, its probability cannot exceed the marginal ADHD power under the same rule; no independence between the two tests is assumed. Neither marginal nor joint replication was performed here.

[Protocol](PROTOCOL.md), [scenarios](power_scenarios.tsv), [sample-size calculations](target_samples.tsv), [source receipt](receipt.json) and [code](calculate.py) record the assumptions. The zero-effect scenario returns the chosen α and each 80% sample-size solution is obtained numerically from the same two-sided power expression. This is a planning calculation, not validation of the chromosome 3 model, a repair of LAVA QC, or evidence of a causal gene.

A subsequent [public-browser count update](finngen_browser_counts/RESULTS.md) verifies 5,559 cases and 489,493 controls in the R13 phenotype table. The earlier assumed 100,000-control scenario above remains historical; the updated calculation is still illustrative and has no variant-level result.
