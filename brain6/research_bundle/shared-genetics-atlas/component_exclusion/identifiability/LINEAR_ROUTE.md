# A conditional summary-statistics route around individual-level access

Checked 6 November 2026. There is a mathematically valid narrow route, but the current releases have not met its prerequisites. This note corrects an overbroad sentence in the construction protocol: correlated items and shared covariates do not themselves prevent subtraction in a common linear model.

## When reconstruction is exact

Let F be the untransformed five-item count, E the exhaustion indicator, and P = F − E. For the same individuals and the same design matrix X (including genotype and all covariates), ordinary least squares gives

`b_P = (X'X)^−1 X'(F−E) = b_F − b_E`.

This identity allows correlated items. It requires identical per-SNP samples, coding, allele orientation, covariates and estimation weights. Separately fitted logistic models, trait-specific mixed-model estimators, different complete-case samples, rank normalization and study meta-analyses do not inherit this guarantee.

The SNP-effect uncertainty must include sampling covariance:

`Var(b_P) = Var(b_F) + Var(b_E) − 2 Cov(b_F,b_E)`.

Two marginal SE columns do not supply the missing covariance. A genome-wide LDSC intercept is not an exact per-SNP covariance replacement. An approximate summary-based reconstruction would need a declared covariance model and sensitivity analysis. For standardized linear phenotypes, first restore raw units using the actual phenotype SDs, subtract, then divide by the modified-score SD. Standardized betas cannot simply be subtracted. Rank-based inverse-normal transformation has no such SD correction.

For a raw 49-item FI average A and raw item S, the corresponding identity is `b_A_without_S = (49 b_A − b_S)/48`; removing S and T gives `(49 b_A − b_S − b_T)/47`. These formulas require the actual matching raw score and item GWAS. They cannot reconstruct the separately transformed Atkins phenotype.

At the additive-genetic covariance level, with all quantities on consistent raw scales, `Cov_g(sleep,P) = Cov_g(sleep,F) − Cov_g(sleep,E)` and `Var_g(P) = Var_g(F)+Var_g(E)−2Cov_g(F,E)`. Divide the covariance by `sqrt(Var_g(sleep) Var_g(P))` to obtain the modified-score genetic correlation. This requires the full relevant covariance information, compatible populations and joint uncertainty; subtracting genetic correlations alone is invalid. It defines an additive genetic linear combination and does not reconstruct binary cases.

## Why the binary outcome cannot generally be recovered

Here is an exact counterexample, not simulated biology. Let K be the number of four physical deficits and E exhaustion. Conditional on K, distribute the K positive items uniformly across the four positions. Consider these two genotype groups:

| Group | Joint (K,E) distribution |
|---|---|
| G=0 | probability 1/2 at (0,0), 1/4 at (2,1), 1/4 at (4,0) |
| G=1 | probability 1/2 at (0,0), 1/4 at (3,1), 1/4 at (3,0) |

Each physical item has probability 0.375 in both groups; exhaustion has probability 0.25 in both. The full five-item count has the identical distribution in both groups: 0 with probability 1/2, 3 with probability 1/4, and 4 with probability 1/4. Thus neither the separate item GWAS nor the full count or full binary GWAS identifies a genotype association in this example. Nevertheless, the exhaustion-deleted binary case probability `P(K≥3)` is 0.25 versus 0.50, an odds ratio of 3. A second world using the G=0 distribution in both groups has the same marginal and full-score distributions but an exhaustion-deleted odds ratio of 1. The unavailable joint item distribution distinguishes these worlds. More computation or embeddings cannot identify it from those marginal summaries alone.

## Computational checks

Run `python3 check_linear_route.py` from this directory, or use its full path. The [script](check_linear_route.py) uses only NumPy and the Python standard library; its [saved output](linear_route_check.json) contains synthetic results only. In a 4,000-row fixed-seed example with correlated binary components and a shared covariate, it verifies:

- Every fitted coefficient obeys raw-score subtraction to 3.2×10⁻¹⁵.
- The reconstructed physical-score SNP SE equals the direct regression SE, 0.0289831. Incorrect independence gives 0.0353631.
- Separately inverse-normalized phenotypes break the naive identity: direct physical-score SNP beta 0.152154 versus subtraction 0.103318.
- The exact binary counterexample preserves all listed marginals while changing the deleted-score case probability.

These are algebra checks, not GWAS results or validation of real data. They establish what additional information would be sufficient and what is not identified.

## What to use now

The existing [explicit physical-function genetic index](../physical_factor_workaround/INDEX_RESULTS.md) remains the usable public-data route: a defined linear combination of walking, grip and activity genetic components with aligned block-jackknife uncertainty. It is a different construct from exhaustion-deleted Fried frailty. The raw-score reconstruction above becomes actionable only if the source provides compatible untransformed full/item results, scale information and sampling covariance, or runs the modified GWAS directly.

A further [Ma et al. source](https://pmc.ncbi.nlm.nih.gov/articles/PMC11979021/) was checked on 6 November. It uses the five-component Fried phenotype, including exhaustion, and states that its UKB data are restricted. The main article does not identify an account-free full modified-frailty GWAS. The PMC workbook link encountered a browser challenge, but the publisher download subsequently succeeded. A [complete 26-sheet audit](../ma_2025/REVIEW.md) found phenotype definitions, exposure analyses and attributable-fraction summaries rather than a full SNP-level GWAS table. This is not a newly verified replacement dataset.
