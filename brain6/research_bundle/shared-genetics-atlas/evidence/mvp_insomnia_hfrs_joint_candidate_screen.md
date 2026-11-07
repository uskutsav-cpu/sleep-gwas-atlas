# Joint candidate screen: MVP insomnia and Finnish hospital frailty

Two of the existing 15 candidate variants have strong **two-sided marginal association evidence in both** the MVP insomnia and Finnish hospital frailty risk score (HFRS) results, using a conservative descriptive joint screen. Their estimated signs agree across those two outcomes. This is a follow-up prioritization, not proof of a shared causal variant or a validated risk score.

| Variant; counted allele | Positional gene label | MVP insomnia p | FinnGen HFRS p | Larger p | Larger p × 15 |
|---|---|---:|---:|---:|---:|
| rs12536335; G | FOXP2 | 4.57 × 10⁻⁸ | 7.85 × 10⁻⁹ | 4.57 × 10⁻⁸ | 6.86 × 10⁻⁷ |
| rs1909650; G | PCDH8P1 | 1.58 × 10⁻⁵ | 0.000386 | 0.000386 | 0.00578 |
| rs4886140; A | RPP40P2 | 0.0154 | 0.00226 | 0.0154 | 0.231 |
| rs7734045; C | AC009014.3 | 0.0184 | 0.000325 | 0.0184 | 0.276 |
| rs323509; A | RP11-6N13.1 | 0.0509 | 0.00425 | 0.0509 | 0.763 |

The complete 15-row screen, alleles, signed coefficients, p-values, and flags are in `mvp_insomnia_hfrs_joint_candidate_screen.csv`. These gene names describe the closest annotated gene body in the earlier positional lookup. They do not identify an effector gene. In particular, rs1909650 has only a nominal FFS p of 0.097 in the original reviewed source, so evidence involving HFRS should not be presented as an FFS association. The estimate for rs12536335 also showed opposing FI cohort directions under an inferred TwinGene allele adapter, which needs resolution.

## What the calculation means

For each candidate, `max(p_insomnia, p_HFRS)` is an intersection–union statistic: under a null that at least one trait has no marginal association, it is a conservative valid p-value if each relevant input p-value is valid for its null. It does not require independent trait tests. Multiplying by 15, capped at 1, is a Bonferroni calculation across these 15 tests and does not require independent variants. Same estimated sign is reported descriptively; the two-sided joint test is **not** a test that effects have a shared direction.

The numerical thresholds are exploratory here, not confirmatory adjusted significance. The 15 markers were chosen after inspecting prior cross-trait findings, the new screening was not preregistered, and independence of selection from both follow-up sources has not been established. These issues may invalidate nominal conditional error control. Exact overlap, ancestry and phenotype definitions also require further verification. FFS and FI are not the same outcome as HFRS; HFRS can include diagnostic codes reflecting drowsiness or vital exhaustion. The sleep GWAS and frailty GWAS can be associated at the same marker due to different correlated causal variants; fine-mapping and colocalization are needed to investigate shared causal signals. Nine of the 15 candidate rsIDs have multiallelic representations in broader reference data, requiring careful allele-level reconciliation.

A non-significant result is not proof of no association. Effects use different phenotype scales; they should not be summed or interpreted as an individual's disease probability. This screen does not estimate prospective risk, assess score performance, or establish prediction across ancestries.

## Follow-up

Freeze an independent analysis plan and candidate set, verify source cohort overlap and variant representations, and obtain an independent sleep GWAS and independent, sleep-excluded frailty phenotype. Compare signed effects with uncertainty and report all tested candidates. If both traits show credible association at a region, assess ancestry-appropriate LD, fine-map each trait and perform colocalization before claiming shared causal variation; then investigate gene and cell-type evidence. Evaluate any 15–30-marker score against observed outcomes in a genuinely held-out genotyped cohort.

## Reproducibility

The calculation reads the existing `integrated_15_candidate_review.csv` (SHA-256 `ee670368bd74993a40bade4bc497df7337fb80de2886b073e33d6cd54e6bf6bd`). It takes the larger of the already harmonized two-sided p-values and multiplies by 15; it neither modifies original source files nor estimates new GWAS effects. The output CSV has SHA-256 `34eaae508866b95b08162fc56bedeaa96027d7b83aab979fee312340178bd2e3`. P-values and allele caveats trace to `integrated_15_candidate_review.md` and `candidate_hfrs_full_gwas_review.md`.
