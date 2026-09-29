# Promotion-rule reassessment: continuous sleep duration

Date: 2026-09-26. Scope: Brain6 power-optimized sensitivity only. This reassessment is outcome-blinded to sleep–frailty and to Brain6 bivariate LAVA results; it uses the user-specified screening alternatives, frozen canonical univariate receipts, and the completed continuous-duration trait-only screen. It does not edit or promote canonical v3, v2, or roundoff outputs.

## Decision rule

The project request specified that a candidate can qualify for a full sensitivity when it materially improves local-h² estimability and gave alternative quantitative examples: at least a 10-percentage-point gain in eligible loci, at least a 25% relative reduction in low-local-h² `NOT_RUN`, or enough improvement to approach the frozen family coverage threshold. The rule is satisfied by one or more of those criteria. The later screen summary narrowed that into a single criterion of 250 additional strict univariate p-gate passes. That narrower rule was not required by the user request and conflates processability with the separate significance gate used to allow bivariate tests.

For this candidate, the qualifying criterion is the explicitly listed reduction in low-local-h² `NOT_RUN`; it passes without relying on the ambiguous phrase “eligible loci.” Long sleep has 1,271 low-local-h² `NOT_RUN` loci; continuous duration has 876, a reduction of 395 (31.08%). The candidate also increases finite local-h² output from 1,204/2,495 to 1,619/2,495 (+415; +16.63 percentage points). These are estimability changes, not evidence of genetic sharing. Its strict univariate p-gate count rises only from 1 to 7 (+6; +0.24 points), and the candidate alone still has 876/2,495 untested loci (35.11%), exceeding the frozen 5% family ceiling of 124 untested loci by 752. The canonical `FAILED_QC_NOT_PROMOTED` decision therefore remains unchanged.

## Pairwise sensitivity readiness

The continuous-duration screen is sufficient to justify a separately named full LAVA sensitivity attempt. This does not authorize canonical promotion or relax the frozen family gate. Cross-trait LAVA requires method-specific overlap matrices. The existing `longsleep__scz`, `longsleep__bipolar`, and `longsleep__parkinson` matrices are specific to the binary long-sleep GWAS and cannot be reused for continuous duration.

Using the candidate strict p-gate loci and the frozen canonical partner univariate rows, two loci (`266`, `1719`) pass both unique-trait gates for duration × schizophrenia; no loci pass both gates for duration × bipolar or duration × Parkinson disease. These are eligibility counts only; no bivariate estimates were inspected or generated. Before any pairwise execution, estimate and hash-bound the continuous-duration cross-trait LDSC intercepts and construct the three new overlap matrices. If valid matrices cannot be obtained, preserve the blocker and do not assume zero overlap or borrow the long-sleep estimates.

The full candidate trait-only result is retained under `lava_trait_screen_v1/`, with four shards covering all 2,495 frozen loci. The next valid step is an isolated, receipt-backed sensitivity family using the same loci, reference, LAVA parameters, four-worker limit, and univariate gates; pairwise testing remains subject to the same strict per-locus gate and family QC. The sensitivity must remain explicitly noncanonical and must report a failed family QC if the 5% missingness limit is not met.

## Consequences

- Recommendation: **SUPPLEMENT WITH A FULL, ISOLATED POWER-OPTIMIZED SENSITIVITY ATTEMPT**; do not replace or reclassify canonical long sleep.
- Association-based selection: none.
- Canonical results, receipts, hashes, checkpoints, v2, and roundoff runs: unchanged.
- Downstream locus follow-up: not permitted unless the sensitivity produces valid QC-passing results under its own frozen family rules.
- Immediate prerequisite: source-bound LDSC overlap estimates and a sensitivity runner/configuration that cannot read or write canonical result directories.
