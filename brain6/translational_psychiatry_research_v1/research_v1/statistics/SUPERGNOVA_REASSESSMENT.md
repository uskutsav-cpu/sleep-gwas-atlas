# SUPERGNOVA source-specific inference reassessment

The three historical covariance results remain numerically reproducible from their archived tables, but their covariance P values are **not certified as calibrated**. New source-matched experiments identify an estimator/variance mismatch in the exact implementation used historically. This changes the permissible inference: keep the three results as reported secondary candidates requiring native and statistical review; do not promote them as robust validated sharing. No historical output or threshold was changed.

## Exact source identity

The public [qlu-lab/SUPERGNOVA source](https://github.com/qlu-lab/SUPERGNOVA/tree/319e84e114a4f954005a4592756c56cfee083667) was acquired at the historically frozen commit `319e84e114a4f954005a4592756c56cfee083667`. The eight acquired files, license, pre-download manifest, byte lengths and SHA256 receipts are retained here. Applying the original Brain6 runtime patch reproduces **all six historical Python-file hashes exactly** (`historical_source_identity.tsv`). In particular, patched `calculate.py` matches `fb31e6c448e921547c571465c0413b45bd84dc8501da4503f1e845e673b1afe4`. The runtime patch preserves the relevant mathematical expressions. A fresh public `git ls-remote` observation found the current upstream HEAD still equals that commit; this is a dated observation, not a promise that upstream will remain unchanged.

The [original method paper](https://link.springer.com/article/10.1186/s13059-021-02478-w) describes eigenspace weighted regression and warns that noisy local heritability destabilizes derived correlations. The source implements the following specific path:

* `Localrho` is returned from the **full unweighted moment** `(sum(Z1*Z2)-c*m)/(meanLD*sqrt(N1*N2))`.
* Its reported variance is selected using **adaptive weighted eigenspace regression** of `y=tZ1*tZ2-c*d`, combining theoretical and empirical weighted-regression variance.
* `cur_v2/cur_v1` is computed internally, but the corresponding weighted point estimate is never assigned to returned `Localrho`.
* After variance selection, the slice `[:len(cur_d)+min_idx-1]` uses one fewer eigenmode than the cumulative selected variance. This slice affects the intercept-uncertainty `v4` term. Our primary simulations set intercept variance to zero, so they cannot evaluate the practical consequences of that additional inconsistency.
* Chromosome-split calls use chromosome-context estimates for the heritability weights; the synthetic experiment supplies four independent equal-sized blocks as context. It does not claim to recover real chromosome heritability.

## Frozen truth-known experiments

`calibration_protocol.json`, the code and source files were hashed in `calibration_freeze.json` **before** any simulation outcomes. The complete locked experiment used seed 2026100711, 160 variants per local block, four context blocks, N1=386533, N2=225534, AR(1) LD parameters 0.2/0.6/0.9/0.97, both local h² values 0.0001/0.001/0.01 and exact overlap intercept 0/0.2. True local genetic covariance was zero. All **24 scenarios × 10,000 replicates = 240,000** simulated tests were executed. No scenario was dropped; no invalid numerical estimate occurred.

These are Gaussian **summary-statistic simulations with exact deterministic LD**, not simulated participant genotypes, empirical GWAS replay, independently sampled psychiatric cohorts, or original reference validation. They exclude real LD-reference sampling noise, binary ascertainment, variant-varying meta-analysis participation, empirical genome-wide intercept estimation and its uncertainty.

For each eigenmode j, the generating variances are `A_j=d_j+N1*h1*d_j²/m` and `B_j=d_j+N2*h2*d_j²/m`, with cross-covariance `C_j=c*d_j` under the null. By Gaussian fourth moments, the exact variance of the returned full-moment estimator is:

`Var(rho_moment) = sum_j(A_j*B_j + C_j²) / (N1*N2*meanLD²)`.

A separate vectorized implementation matched the unchanged upstream numerical routine on **12 concordance fixtures** with known dense LD (`upstream_concordance.tsv`). The upstream routine was called with a synthetic LD API, bypassing real genotype reading; this is method-routine concordance, not a native full-pipeline integration test. The independent analytic variance matched simulated moment variance to within 5.05% across all scenarios; heavy tails explain substantial variance-estimation Monte Carlo uncertainty in the most concentrated scenarios.

| Illustrative scenario, no overlap | Exact moment variance / mean reported variance | Reported P rejection at α=0.05 | Exact-variance moment normal-Wald rejection | Oracle GLS normal-Wald rejection |
|---|---:|---:|---:|---:|
| AR=0.2, h²=0.0001 | 0.946 | 4.08% | 4.68% | 4.78% |
| AR=0.6, h²=0.001 | 1.842 | 14.72% | 5.61% | 5.10% |
| AR=0.9, h²=0.001 | 3.726 | 27.12% | 5.41% | 5.10% |
| AR=0.97, h²=0.01 | 19.047 | 52.16% | 6.12% | 4.99% |

The complete numerical table includes Wilson 95% intervals, every scenario and rejection counts at all three locked thresholds (`calibration_simulations.tsv`). Across the 24 scenarios, nominal 5% rejection ranged **4.08–52.16%**. These counterexamples demonstrate that the returned source-specific point estimate and uncertainty cannot be assumed calibrated merely because the software emits a finite P value. They do not establish how much bias or inflation occurred in any actual Brain6 block.

The frozen `0.05/8465` threshold was also reported without modification. In the most severe illustrative scenario, 2202/10000 null tests crossed that nominal per-test threshold. This is a synthetic per-test result, **not an estimate of historical Brain6 FWER** or the number of false historical discoveries. Ten thousand replicates also provide weak precision at a genuinely well-calibrated ~5.9e-6 tail; zero events do not certify tail calibration.

## Remaining normal-tail uncertainty

The oracle benchmarks use correct variance under the exact model, but normal Wald inference is not automatically exact. The post-result deterministic appendix `gaussian_product_tail_diagnostics.tsv` computes variance concentration and exact excess kurtosis. For the returned moment estimator, effective variance-mode count ranges **1.35–112.10**, and exact excess kurtosis ranges **0.054–4.456**. In the concentrated AR=0.97, h²=0.01 null, only about 1.35 effective modes dominate moment variance. That explains why an exact-variance normal-Wald benchmark can itself reject at 6.12%, above nominal 5%. The oracle GLS distributes variance over more modes in most scenarios but retains finite-model tail limitations. Correcting only the variance is insufficient to certify the stringent empirical family tail.

This deterministic appendix also gives exact local-h² estimator SE under the illustrative model. It is not an empirical Brain6 h² confidence interval; real LD and joint uncertainty are absent.

## Empirical block diagnostics and fragility

The archived block-table SHA256 remains `b56c4fd2d2a93c34da6b057cf46c516d51f587a60bef9f0d7a236bc1492b0249`. The original family denominator remains **8465**, and the critical absolute normal Wald Z is about **4.5297**.

| Pair | Estimated blocks | Both local h² positive | Missing derived r | Derived r outside [-1,1] | Historical family-significant blocks |
|---|---:|---:|---:|---:|---:|
| Insomnia–ADHD | 1646 | 1063 | 583 | 358 | 1 |
| Insomnia–MDD | 1658 | 1136 | 522 | 365 | 2 |
| Total | 3304 | 2199 | 1105 | 723 | 3 |

The three historical significant results would lose their original family threshold if their SEs increased by factors **1.2371** (chr11 ADHD), **1.0538** (chr6 MDD) or **1.4011** (chr11 MDD), holding the reported point estimate fixed. This is a retrospective fragility calculation, not a corrected SE or recalibrated P value. Synthetic variance ratios must **not** be substituted for the unknown empirical correction.

The insomnia–ADHD chr11 point h²/covariance matrix violates the Cauchy–Schwarz constraint and yields derived r=1.1037. Its covariance exceeds the point-h² bound by only about **0.527 reported covariance SE**. Because the joint h²/covariance uncertainty is unavailable, this quantity is a diagnostic, not a determinant significance test. No clipping, bounded-correlation claim, or assumption of a biological impossibility is warranted. Noisy unconstrained moment estimators can produce out-of-range derived r even with valid true parameters.

## N and overlap sensitivity: what is and is not executable

For fixed signed Z and fixed reference, changing **scalar** N merely rescales local/global h² and rho. The `N*h` products used in variance and intercept weights stay invariant; rho rescales by `1/sqrt(N1*N2)`, variance by `1/(N1*N2)`, while covariance Wald P and derived r remain invariant. Four numerical counterexamples verify this identity (`scalar_n_invariance.tsv`). Thus unchanged P values under total-versus-effective scalar-N substitutions **cannot validate source sample-size semantics**. This algebra does not cover SNP-varying participation, changing Z calibration, ascertainment or cohort composition.

Holding Z and LD fixed, the moment estimator responds to intercept perturbation as `delta_rho = -m*delta_c/(meanLD*sqrt(N1*N2))`. The original pair receipts, empirical intercept estimate/variance, block LD eigenvalues and meanLD are absent. Therefore an actual scientifically admissible empirical overlap-sensitivity or calibrated covariance reanalysis could not be executed. The simulations supplied exact overlap truth and cannot validate the real estimated intercept.

## Research decision and restart gate

**Three historical covariance positives: REPORTED_SECONDARY_CANDIDATES_UNCERTIFIED_CALIBRATION.** They are neither independently validated positives nor credible biological nulls. The exact implementation provides a reproducible calibration counterexample; actual Brain6 effect direction and magnitude remain unresolved beyond the archived source-defined estimates.

The scientific restart requires the original harmonized GWAS rows and per-variant source/sample semantics, exact genotype/reference and SNP ordering, original genome-wide intercept/variance and pair receipts, native reproducibility, and an independently reviewed estimator/variance pairing. Any correction belongs in a new versioned protocol with the complete frozen family, original outputs retained, properly calibrated tails and fresh independent validation. Do not select or repair only the three observed significant blocks.

Independent human statistical-genetics review remains required. Separate automated derivations can corroborate a numerical counterexample; they do not constitute independent human validation or psychiatric cohort replication.
