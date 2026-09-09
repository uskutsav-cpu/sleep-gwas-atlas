# Track B Pair B selection rationale

## Frozen decision

Pair B is **insomnia ↔ ADHD**. This identity was frozen before any Track B local, pleiotropy, fine-mapping, colocalization, tissue, cell, or molecular-QTL result existed.

- Pair B identity SHA-256 (`B\tinsomnia\tadhd\n`): `951158c99371ccec539fd4842ec1fb095b18c10ce1b66521c6bd50e3bad0a623`
- Complete 72-row candidate table SHA-256: `03caed93ad94eccc483ef9d86caede91a57ccbeca8eef74276ef346313aa0bde`
- Three-role pair manifest SHA-256: `b3c520424415895df18b87141d9b99721d3bf49330844524731dd563da12dd4f`

## Eligible family

The locked core contains 12 sleep traits and six eligible diagnosed brain/psychiatric disorders (Alzheimer's disease, Parkinson's disease, MDD, schizophrenia, bipolar disorder, and ADHD), yielding 72 candidates. The existing extension contains psychosocial traits but no additional diagnosed brain-disorder phenotype meeting this definition, so it contributes no eligible Pair B row.
Non-significant candidates were retained in the complete table and explicitly marked `NOT_AUDITED_NONDISCOVERY`; the pre-existing literature audit was discovery-only and cannot support pair-specific novelty conclusions for those rows.

## Prospective selection procedure

Candidates were filtered lexicographically, not assigned an opaque score:

1. clean primary Phase-1 status and PASS h² QC for both traits;
2. locked core FDR < 0.05;
3. dense signed GWAS available for both traits;
4. exact-concept phenotype match in the existing evidence map;
5. concordant direct rg evidence from a changed/expanded dataset with at least partial independence;
6. among survivors, descending |rg|, then ascending P, then stable trait identifiers.

Generic human single-cell and scATAC resources were recorded as feasibility evidence only. They were not treated as pair-specific cell/QTL validation, and no downstream result was inspected.

## Ordered survivors

| Rank | Sleep trait | Brain trait | rg | FDR | Existing replication class |
|---:|---|---|---:|---:|---|
| 1 | insomnia | adhd | 0.4033 | 2.568e-35 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |
| 2 | sleep_apnea | adhd | 0.3278 | 5.493e-22 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |
| 3 | sleep_apnea | bipolar | 0.1649 | 1.676e-08 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |
| 4 | chronotype | scz | -0.1299 | 2.52e-11 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |
| 5 | insomnia | bipolar | 0.1113 | 7.512e-05 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |
| 6 | chronotype | bipolar | -0.0529 | 0.03832 | DIRECTIONAL_REPLICATION_PARTIAL_INDEPENDENCE |

## Why Pair B won

Insomnia–ADHD was the leading survivor because it combined clean h²/QC, an extremely strong locked global signal, dense inputs for both traits, exact-concept prior evidence, and the largest absolute rg among candidates meeting all preceding gates. This is not a novelty claim: direct rg has already been reported, and the prior comparison is only partially independent/overlapping.

## Required next gate

The existing evidence is not a fully independent replication. Pair B cannot advance beyond a global finding until a fresh cohort-overlap audit and, if feasible, standardized LDSC against a genuinely independent external ADHD GWAS are completed. A materially incompatible, well-powered opposite-direction result is a no-go for mechanistic escalation.

## Non-selection caveats

MDD candidates were not eligible for the dense-input gate because the local file is HapMap3-prefiltered. Parkinson candidates likewise lack dense downstream-ready input. Literature saturation was retained as a caveat rather than used post hoc to override the prospective gates. No candidate may replace Pair B because a later analysis is easier or more positive.
