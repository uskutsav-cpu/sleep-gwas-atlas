# Continuous-duration sensitivity: overlap and readiness

Audit date: 2026-09-26  
Status: `PASS_OVERLAP_AUDIT_FAMILY_FAILS_FROZEN_QC`

## Scope

This is a candidate-specific LDSC overlap-input audit for the already completed continuous sleep-duration local-h² screen. It is separate from canonical LAVA v3. No sleep–disorder association results were inspected to select the candidate or estimate these overlap inputs.

The source-bound LDSC run used the continuous-duration candidate and the frozen SCZ, bipolar-disorder, and Parkinson-disease summary statistics. Its three generated covariance matrices and LDSC logs pass the integrity validator. The DOI-registered Zenodo copy of the European LDSC reference archive passed its published MD5 check and the recorded SHA-256 check. Per-pair intercepts and artifact hashes are recorded in `overlap_integrity_audit_v1.json` and the external `overlap.provenance.json` receipt.

## Family readiness

The candidate improved finite local-h² processability over canonical long sleep by 415/2,495 loci (16.63 percentage points) and reduced low-h² `NOT_RUN` loci from 1,271 to 876 (31.08%). It gained six loci through the unchanged strict univariate gate. These outcome-blinded technical gains satisfy the user's predeclared example thresholds for investigating a power-optimized sensitivity; they do not change the canonical decision or constitute evidence of genetic sharing.

Under the frozen seven-trait substitution-family accounting, the candidate yields 3,305 `NOT_RUN` cells among 17,465, or 18.92%, against a ceiling of 873 cells (5%). The family fails the unchanged QC rule. Accordingly, the bivariate LAVA stage and downstream association or biological follow-up were not started. There are two candidate/disorder strict-gate locus intersections for SCZ (loci 266 and 1719), and none for bipolar disorder or Parkinson disease; these are eligibility intersections only, not local-rg findings.

## Interpretation

Canonical v3 remains `FAILED_QC_NOT_PROMOTED`, immutable. The continuous-duration screen is a technically justified sensitivity candidate, but it does not pass family QC and cannot replace the canonical trait or support downstream conclusions. The new LDSC matrices document candidate-specific sample-overlap covariance inputs; they do not authorize or imply bivariate results.

## Reproduction

- Integrity audit: `python3 brain6/scripts/validate_power_optimized_lava_overlap_v1.py`
- Outcome-blinded power comparison: `promotion_rule_reassessment_2026-09-26.json`
- Matrix and source receipts: external run `overlap_v1_retry/overlap.provenance.json`
- Integrity decision and source hashes: `overlap_integrity_audit_v1.json`
