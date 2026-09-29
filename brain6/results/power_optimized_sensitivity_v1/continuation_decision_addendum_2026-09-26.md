# Continuation decision addendum — 2026-09-26

This addendum reconciles the initial candidate screen with the user's full Phase 5 promotion rule. The earlier `README.md` and `candidate_replacement_screen_complete_v2.tsv` are preserved as hash-bound historical screen outputs; the completed candidate aggregation receipt binds that README, so it remains unchanged.

## Technical advancement and family QC are separate decisions

The user's Phase 5 examples are disjunctive: a candidate advances for further technical sensitivity evaluation if it improves finite local-h² processability by at least 10 percentage points **or** reduces low-local-h² `NOT_RUN` by at least 25%. Dashti 2019 continuous duration meets both: +16.63 percentage points and a 31.08% reduction. Its strict univariate gate-pass count rises by only six, below the older +250 strict-pass example; that does not negate the independently specified processability criteria.

Candidate-specific LDSC overlap preparation was completed and passes integrity QC. The substituted seven-trait family has 3,305/17,465 `NOT_RUN` cells, exceeding the frozen 873-cell limit. Therefore no bivariate LAVA result was generated or promoted, and no downstream locus is eligible. The canonical v3 decision remains immutable and failed QC.

The source-bound current decision table is `candidate_replacement_screen_complete_v3.tsv` with provenance in its adjacent JSON. It records one technical advancement candidate and zero candidates passing the family-wide QC gate. The source-bound lower-bound audit is `sensitivity_family_feasibility_v1.json`: even perfect eligibility for both sleep traits leaves 1,758 `NOT_RUN` cells across the five disorders, 885 over the family ceiling. Further sleep-only replacements cannot make this frozen family pass.

## Current recommendation

- Insomnia: **KEEP** the canonical UK Biobank phenotype; no accessible equivalent with demonstrated material local-h² improvement is available.
- Long sleep: **SUPPLEMENT WITH SENSITIVITY** using the continuous-duration trait for the completed trait-only screen and overlap preparation. Retain the canonical ≥9-hour result; do not call the sensitivity a replacement or independent replication.
- Full bivariate sensitivity and downstream biology: **BLOCKED BY FROZEN FAMILY QC**. Do not relax thresholds or launch gated downstream analyses.

No downstream sleep–frailty association results were consulted for these decisions.
