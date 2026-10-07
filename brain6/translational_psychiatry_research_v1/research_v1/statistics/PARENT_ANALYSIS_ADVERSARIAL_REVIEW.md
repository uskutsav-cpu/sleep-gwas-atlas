# Separate adversarial review of additional parent analyses

This is a separate automated review using independently written numerical implementations. It is not independent human validation, participant-level verification or independent cohort replication. No parent frozen code, protocol or outputs were modified.

## Global psychiatric contrasts

The full72-family and its input/code hashes were verified. All **48 unique psychiatric marginal rg rows** used by the **72 contrasts** satisfy the frozen positive h²/SE and h²-Z≥4 rule. Separate stdlib `erfc` tails and `NormalDist` quantiles reproduce all72 differences, maximum SEs, P values, Bonferroni72 values and simultaneous confidence bounds, with maximum absolute numerical discrepancy **3.34e-16**. Exactly **16** contrasts cross the frozen bound. The complete sensitivity grid never produces a larger P than the stated maximum-covariance P.

The mathematical claim is correct under calibrated joint asymptotic normality: `Var(rg1-rg2)=s1²+s2²-2c≤(s1+s2)²`. Even if pairwise error correlation -1 is incompatible with every pair simultaneously in a global covariance matrix, the larger pairwise admissible set gives a conservative bound for each true covariance. Bonferroni72 then controls a normal-model family irrespective of dependence. A bound based on estimates and plug-in SE is asymptotic, not an exact finite-sample theorem.

The reported sufficient 80% normal-power difference `(z_family+z_0.8)*(s1+s2)` is a sufficient normal-model signal magnitude: at worst-case SE, the correct-direction rejection tail alone has ≥80% probability, with the opposite tail adding probability. It is not actual cohort prediction power or a power estimate after accounting for ascertainment and LD/source uncertainty.

The16 contrasts compare source-defined rg parameters, principally ADHD/MDD versus schizophrenia/bipolar across correlated sleep phenotypes. They do not constitute16 independent biological discoveries. Marginal h² QC does not prove rg/SE calibration, remove source phenotype/cohort heterogeneity or validate disorder-specific causal architecture. No-significant-bound contrasts are not evidence that the disorders are identical. The historical outcome selection and same-source nature are correctly disclosed. Novelty remains a separate literature question.

## Original molecular posterior sensitivity

The four input hashes, retained coloc5.2.3 reference-source hash, frozen original-code hash and corrected-code amendment were verified. Independent least squares through the origin reproduces sdY. A bounded explicit **all i≠j configuration sum** reproduces H3 without relying on the parent's prefix/suffix arithmetic. All **60** posterior conditions agree within **3.45e-15** maximum absolute error. Default conditions reproduce all four historical posteriors within the frozen1e-10 tolerance. Maximum H4 remains **0.6175746925**, with **zero conditions ≥0.8**.

The unchanged retained R reference functions `sdY.est`, `approx.bf.estimates` and `combine.abf` were also executed directly in base R on the same60 conditions. Maximum cross-language posterior error is **8.39e-15**. This is native reference-function numerical execution on actual inherited molecular slices; the complete coloc package's orchestration/checking code was not installed or rerun. The R warnings about estimating sdY are retained in the execution log.

The serialization amendment changes explicit float/int conversion and hash-amendment validation. The three mathematical helper functions are AST-identical. It does not change priors, inputs, the posterior model or the60-condition family.

The arithmetic follows the pinned coloc equations, conditional on log-OR effect scale for the ADHD GWAS, a quantitative molecular trait, the assumed sdY relationship and one causal variant per trait. The p12 grid reaches an extreme p12=p1=p2 condition; it is a prior diagnostic, not an empirical probability calibration or a biological positive. Posterior H3 dominance does not prove distinct causal variants if multiple signals, missing variants or other model assumptions fail. All tests concern only the ADHD component of one previously defined region; there is no insomnia molecular posterior and no independently established sleep–ADHD effector.

ABFs depend on beta². Therefore allele-sign reversal leaves this posterior unchanged: a successful posterior replay cannot certify effect-allele/strand harmonization or genetic effect direction. The retained slices omit the upstream records needed to verify exact build/liftover and source-strand identity independently.

## Palindromic-key stress test

The historical protocol says to exclude "ambiguous or mismatched alleles" while matching unique lifted GRCh38 positions and unordered alleles. Retained slices include palindromic A/T and C/G keys, including the highest conditional-H4 variant `chr5_88551455_T_A` in some defaults. This warrants disclosure and a robustness test. **Palindromic is not automatically scientifically ambiguous when explicit reference/alternate identity is established.** The wording does not prove that the historical protocol required removal of every palindromic site, and the retained slices do not prove corrupt variant identity.

The separate `HARMONIZATION_STRESS_V1` protocol and code/dependency/input hashes were frozen before filtered posterior outcomes. It excludes every palindromic key regardless of MAF, re-estimates sdY on the retained slice and retains all **4×15=60** conditions. Source N/P input QC was checked before execution. Full original inputs and all60 original posterior conditions are preserved.

Each cortex context removed **251/1624** records, retaining **1373**; each frontal context removed **236/1505**, retaining **1269**. All four contexts still meet ≥500 variants, minGWAS P≤5e-8 and minQTL P≤5e-8. No malformed keys were present. These are context-record counts, not independent SNP counts or evidence of251 wrongly mapped sites.

| Context | Original default H4 | Nonpalindromic default H4 |
|---|---:|---:|
| Cortex eQTL ENSG00000271904 | 2.17855e-6 | 2.67765e-6 |
| Cortex sQTL ENSG00000247828 | 3.15548e-8 | 4.84756e-8 |
| Frontal eQTL ENSG00000250377 | 9.38877e-6 | 1.10194e-5 |
| Frontal eQTL ENSG00000271904 | 0.0823852 | 0.0155614 |

Across the new60 conditions, maximum H4 is **0.4552597501**, and no condition reaches0.8. Thus the absence of strong H4 under this limited single-signal model persists under a stricter identity-ambiguity stress rule. It does not establish absence of molecular involvement or resolve original allelic identity. Removed-variant identities, exact SHA256s, all posterior conditions and input-QC outcomes are saved in `molecular_palindromic_diagnostic/`. The four filtered slices retain ADHD source beta/SE and are **local only** under ignored `work/ld_genotypes_research_v1/molecular_palindromic_LOCAL_ONLY/`; none belongs in Git or a public handoff. A storage-only amendment retains the original frozen code, records unchanged scientific settings/results, and redirects future raw-slice output to that local directory. The local integration check requires the actual bytes and fails if missing; it is separately selected from the source-free aggregate checks.

## Parent independent covariance-variance check

Both100000-replicate parent cases are explicitly post-result counterexample reviews. The trace/eigen identity is correct because cross-covariance C=cR is symmetric; in a general asymmetric cross-covariance case the relevant term is trace(C²), not trace(CCᵀ). The supplied model satisfies the symmetric case.

A separate stacked-Gaussian quadratic-form calculation constructs Ω for [Z1,Z2] and `H=[[0,I/2],[I/2,0]]`. It reproduces `Var(ZᵀHZ)=2trace[(HΩ)²]` and matches the parent's analytic covariance variances within **6.54e-16 relative error**. The parent's Monte Carlo sample variances are within **0.540 estimated MCSE** of these independently derived values. Its estimated MCSE from squared, mean-zero covariance draws is an appropriate large-replicate diagnostic; it is not an exact confidence interval for an empirical GWAS variance.

This corroborates the mathematical moment-variance counterexample. The simulation generation remains an ideal Gaussian eigenmode model, and finite normal-Wald tails can remain heavy even with correct variance. It neither certifies empirical Brain6 covariance P values nor estimates a historical false-positive rate.

## Decisions

All three reviewed computations pass separate numerical checks. The scientific inferences remain bounded: **16 retrospective source-parameter contrasts**, **no strong molecular H4 in either60-condition grid**, and **a corroborated implementation-specific calibration concern**. No causal psychiatric mechanism, validated clinical relevance, independently replicated two-trait locus, or robust calibrated empirical local-covariance finding follows from these checks. Independent human statistical-genetics review and native empirical source contracts remain necessary before promotion.

Reproducibility:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/review_parent_analyses.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/homebrew/bin/Rscript brain6/translational_psychiatry_research_v1/research_v1/statistics/review_coloc_reference.R /absolute/path/to/sleep-gwas-atlas
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/molecular_palindromic_diagnostic.py run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python -m pytest -q -m 'not integration' brain6/translational_psychiatry_research_v1/research_v1/statistics/test_parent_review.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python -m pytest -q -m integration brain6/translational_psychiatry_research_v1/research_v1/statistics/test_parent_review.py
```

The palindromic protocol already exists and must not be re-frozen. Commands write only their new versioned diagnostics. Frozen parent files and historical inputs remain untouched.
