# Frozen LAVA v3 failure decomposition

The source is the immutable v3 family decision and its checksum-bound 17,465-row aggregate (run `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`). The accompanying `lava_v3_failure_matrix.tsv` enumerates **17,465 observed trait-locus cells plus 12,475 derived pair-locus gate rows**. `evidence_origin` distinguishes observed univariate data from inferred pair eligibility. Pairwise LAVA was never executed in v3; no pair row is an observed local genetic correlation.

## Decisive family gate

Of 17,465 planned univariate cells, 13,745 were `TESTED`, 3,720 were `NOT_RUN`, and zero were `FAILED`. The frozen maximum is 873 untested cells (5%). Observed untested fraction is 21.30%, exceeding the count cap by **2,847**. The terminal canonical decision is `FAILED_QC_NOT_PROMOTED`; it forbids pairwise execution and promotion. All 2,495 loci and all seven traits have an aggregate status row. There is no observed missing-locus or coordinator-truncation category in the final aggregate.

| Ranked cause | Untested cells | Share of untested | Interpretation |
|---|---:|---:|---|
| Low local heritability support | 3,564 | 95.81% | LAVA did not return a supported local-h² test under the locked input, reference, and threshold. This is an estimability outcome, not proof of no heritability. |
| Fewer than minimum shared reference variants | 154 | 4.14% | Sparse reference/input intersection, concentrated in 24 loci; 20 adjacent locus IDs 950–969 are common across all traits. The source-side gap check below localizes these 20 to absent normalized GWAS variants, not an empty LD reference. |
| Fewer than minimum LAVA components | 2 | 0.05% | Bipolar and schizophrenia at locus 1416. |

Even a perfect technical repair of every 156 minimum-K cell would leave 3,564 untested, exceeding the cap by 2,691. The five disorder traits alone contribute 1,758 untested cells; hypothetically making both sleep traits perfectly testable still exceeds the full-family cap by 885. Therefore a sleep-only substitution, reference-only cleanup, or retry of the same bytes cannot qualify the fixed seven-trait family.

## Trait and pair detail

| Trait | Planned | Tested | Not run | Low h² | Shared-ref K | Component K | Strict univariate gate passes |
|---|---:|---:|---:|---:|---:|---:|---:|
| ADHD | 2,495 | 2,272 | 223 | 200 | 23 | 0 | 88 |
| Bipolar disorder | 2,495 | 2,327 | 168 | 145 | 22 | 1 | 64 |
| Insomnia | 2,495 | 1,824 | 671 | 651 | 20 | 0 | 0 |
| Long sleep | 2,495 | 1,204 | 1,291 | 1,271 | 20 | 0 | 1 |
| MDD | 2,495 | 1,906 | 589 | 565 | 24 | 0 | 14 |
| Parkinson disease | 2,495 | 1,859 | 636 | 613 | 23 | 0 | 12 |
| Schizophrenia | 2,495 | 2,353 | 142 | 119 | 22 | 1 | 167 |

The pair rows are *derived* using the fixed trait-only p gate `p < 0.05/17,465`. Each pair has 2,495 planned slots, zero executed v3 bivariate slots, and no observed pairwise convergence, rg, SE, or numerical result.

| Pair | At least one univariate NOT_RUN | Both tested but gate missed | Both strict gates passed | Bivariate executed |
|---|---:|---:|---:|---:|
| Insomnia–ADHD | 808 | 1,687 | 0 | 0 |
| Insomnia–MDD | 1,072 | 1,423 | 0 | 0 |
| Long sleep–SCZ | 1,342 | 1,152 | 1 | 0 |
| Long sleep–bipolar | 1,360 | 1,135 | 0 | 0 |
| Long sleep–Parkinson | 1,582 | 912 | 1 | 0 |

## Technical, statistical, and input classification

- **Process and numerical failure:** zero observed `FAILED` trait cells in v3. The final aggregate and decision report no univariate exception, nonconvergence or numerical failure. Singular covariance, unstable local rg/SE, and bivariate convergence cannot be assessed for v3 because the bivariate stage was not authorized or run. The bounded pilot checked source fields, positive finite N, finite Z, alleles, duplicate SNP IDs and hashes at four loci; it is not a genome-wide proof of those properties.
- **Reference processability:** 154 shared-reference minimum-K outcomes and two component-K outcomes. The diagnostic pilot at locus 950 found zero materialized sumstat variants for each trait; at locus 1416 several traits had 0–3. A follow-up read-only source/reference check (`pilot_locus_950_reference_source_gap.json`) found 1,907 sealed-reference variants in locus 950 and 43,819 across loci 950–969, while all seven normalized GWAS SQLite files had zero variants in those bounds. Immediately adjacent loci 949 and 970 have hundreds or thousands of variants in those same normalized inputs. Each normalization receipt records `exclude_regions=[]` and `input_rows=retained_rows`, so the shared gap predates this LAVA materialization and appears in the harmonized source inputs. Whether that upstream MHC removal was intentional is unresolved; filling it would require a new source/normalization protocol and cannot be assumed to repair the main low-h² failure.
- **Low local-h² support:** 3,564 cells, distributed over the traits above. Locus 1 demonstrates mixed testability with input fields/hash-valid materialization; all seven trait outputs at locus 2207 are testable controls. The pilot does not change LAVA estimates.
- **Scientific input caveats:** the seven sources are hash-bound in the rescue manifest, but MDD's original build remains recorded as `UNRESOLVED` in the existing metadata and the pairwise sample-overlap matrices are LDSC-intercept-derived rather than participant-verified. Neither is established as the cause of the v3 univariate missingness; both require separate scientific review before strong local-rg claims. No invalid effect allele, effective-N, cohort, or overlap correction was demonstrated by this audit. Unknown global validity must not be described as confirmed valid.

## Prior analyses compared diagnostically

The earlier v2 baseline finished locus receipts but failed aggregation because pair-context univariate p-values disagreed for the same trait/locus; it also had 2,833 unobserved unique trait-locus cells (16.22%), beyond the same 5% cap. A separate roundoff run handled numerical asymmetry but retained the duplicate-context aggregation issue; neither is a completed promotable family. Canonical v3 resolved the pair-context identity problem by testing each trait independently, then failed for sparse local-h² support. No completed reduced-family result meeting the frozen five-pair denominator was found; changing the family would be a new sensitivity, not a repair.

The previously completed Dashti continuous-duration *trait-only* screen increased processable long-sleep-related loci from 1,204 to 1,619 and reduced low-h² misses by 395. The hypothetical seven-trait substitution family would still have 3,305/17,465 untested (18.92%) versus the 873 maximum. Its candidate-specific overlap matrices passed integrity checks, but no bivariate family was run. These results do not license trait replacement in the frozen five pairs.

The separate `lava_rescue_v1_protocol.md` was fixed before the four-locus receipt/input diagnostic. Its pilot receipts all pass their specified checks and reveal no correctable large-scale technical cause. The conditional full rescue gate is not met, so no new LAVA statistics were generated. See `lava_rescue_v1_qc.json` for the terminal rescue classification.
