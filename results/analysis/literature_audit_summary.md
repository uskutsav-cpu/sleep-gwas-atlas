# Phase-1 literature novelty audit

All **153** locked primary discoveries were reviewed through four fixed Europe PMC query lenses as of **2026-08-28**, then cross-checked against reviewed sleep-GWAS supplements and selected full-text or main-table sources.

- Direct prior genome-wide rg evidence: **97 pairs**
- No direct prior rg located: **56 pairs**
- Cautiously classified apparently unreported: **6 pairs**
- Search queries / screened candidate records: **612 / 1195**

## Classification counts

| Classification | Pairs |
|---|---:|
| APPARENTLY_NOVEL | 6 |
| DIRECT_RG_PREVIOUSLY_REPORTED | 77 |
| DIRECT_RG_REPLICATION_DIFFERENT_DATASET | 20 |
| MR_ONLY | 30 |
| NO_DIRECT_RG_FOUND | 7 |
| OBSERVATIONAL_ONLY | 7 |
| RELATED_GENETIC_EVIDENCE_ONLY | 5 |
| UNCERTAIN | 1 |

## Replication summary

Among 86 pairs with numeric exact-concept comparators, the pair-level correlation between atlas rg and the median published rg was **0.822471832819**; mean absolute difference was **0.0697215116279** and direction concordance was **0.883720930233**.

These comparisons are positive controls, not independent replication: most reuse or overlap public GWAS inputs.

## Apparently unreported candidates

- `sleep_apnea__healthspan`: rg=+0.4005, FDR=3.08e-20, MODERATE_CONFIDENCE
- `shortsleep__parental_lifespan`: rg=-0.2976, FDR=2.02e-16, HIGH_CONFIDENCE
- `snoring__healthspan`: rg=+0.2927, FDR=1.9e-14, HIGH_CONFIDENCE
- `sleep_apnea__parental_lifespan`: rg=-0.2617, FDR=5.58e-12, MODERATE_CONFIDENCE
- `longsleep__parental_lifespan`: rg=-0.2395, FDR=1.33e-08, HIGH_CONFIDENCE
- `sleep_efficiency__frailty`: rg=-0.1520, FDR=0.00012, HIGH_CONFIDENCE

## Interpretation limits

Search hits were screened rather than accepted automatically. Co-occurrence is not direct rg evidence. `APPARENTLY_NOVEL` means apparently unreported under this dated, recorded search—not first-ever. MR, observational association, and global genetic correlation are kept distinct, and contradictions are retained in the audit.
