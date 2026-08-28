# Discovery extension adversarial review

Generated: 2026-08-28T06:20:19Z

Verdict: **PRE-RESULT DESIGN REVIEW COMPLETE; FINDINGS-LEVEL REVIEW BLOCKED BY ACQUISITION**

## Scope and non-result boundary

The immutable core verifier reports: `CORE_CHECKPOINT_OK commit=22df92c8d6894f93c50561ae8a86acf2f48f50cb traits=45 pairs=396 result_sha256=161756ac61775ad3393572cf5fbdf51690eb47ff010b8a50dbcc5fabf79bc077`

The extension panel remains the pre-result lock of 100 traits and 1200 planned sleep-by-extension tests. This review challenges the extension protocol, software readiness, provenance, and current blocking conditions. It does not review biological findings that have not been generated.

Present real downstream finding classes: none.

Absent finding classes: extension h2, global rg, pair novelty, replication, local architecture, pleiotropic loci, fine-mapping/colocalization, mechanistic synthesis.

## Blocking finding

The exact compressed source family is 214.705 GiB and the locked 1.15 safety factor requires 246.91 GiB free. The recorded preflight found 5.42 GiB free, a 241.49 GiB shortfall. No bulk download was started. The current volume cannot support the required retained full-resolution inputs, so substituting HapMap3-only data or a result-selected subset would violate the fine-mapping/colocalization and selection contracts.

No threshold has been weakened to make this blocker disappear.

## Challenge ledger

| Risk | Required control | Current status | Evidence |
|---|---|---|---|
| Core drift | Re-run byte-level core checkpoint before and after every extension stage. | PASS | CORE_CHECKPOINT_OK commit=22df92c8d6894f93c50561ae8a86acf2f48f50cb traits=45 pairs=396 result_sha256=161756ac61775ad3393572cf5fbdf51690eb47ff010b8a50dbcc5fabf79bc077 |
| Post-result panel selection | Panel membership/order/hash must predate extension rg. | PASS | 2026-08-28T04:13:05Z |
| Multiple-testing leakage | Never combine the extension BH family with the immutable 396 core tests. | PASS_CONTRACT | extension_harmonization_policy.json |
| Weak h2 | Exclude rerun h2 Z<4 from primary rg without replacement. | PENDING_REAL_H2 | absent |
| LDSC intercept inflation | Exclude rerun intercept>1.2 from primary rg and retain sensitivity status. | PENDING_REAL_H2 | absent |
| UK Biobank sample overlap | Inspect cross-trait intercepts and disclose overlap; do not equate LDSC adjustment with independent replication. | PENDING_REAL_RG | absent |
| Sparse or proxy phenotypes | Retain exact phenotype definitions and distinguish medication/proxy traits from diagnoses. | PASS_METADATA | candidate_traits.tsv |
| Panel-level novelty inflation | Require pair-level direct/same-phenotype/same-direction/same-sleep-context audit. | PENDING_PAIR_AUDIT | absent |
| Replication non-independence | Require non-overlapping participants and separately sourced summary statistics. | PENDING_REPLICATION | absent |
| Global-to-local overreach | Do not call global rg evidence of a shared locus; retain a prespecified globally-null secondary local set. | CODE_READY_INPUTS_BLOCKED | absent;present:sha256=37888988af5fdf6148d701477c388b4f07857d97565f0db78d8257e7adff6b49 |
| Pleiotropy or mediated effects | Evaluate horizontal, vertical/mediated, shared-factor, and sample-overlap alternatives; run PLACO+ on genome-wide data only. | CODE_READY_INPUTS_BLOCKED | absent;present:sha256=5bed771d4e8579d14785ac59d7798b815270f0da252d4238b13f6d54f04b61d7 |
| Colocalization overclaim | Report H0-H4, priors, sensitivity, and claim guards; colocalization is not causality. | CODE_READY_INPUTS_BLOCKED | absent;present:sha256=10a13727ffbe9c7f9ca07a2d1a53fe2489ae17d36f86a6354002d58cbfad4328 |
| Synthetic/real result contamination | Synthetic tests stay under synthetic or temporary paths and carry explicit markers. | PASS | seven isolated synthetic workflows |
| Storage-driven partial acquisition | Do not silently analyze a result-selected subset of the locked panel. | PASS_BLOCKED | BLOCKED_INSUFFICIENT_STORAGE |
| Local LD-reference mismatch | Require checksum-locked ancestry-matched LAVA/HDL-L references; do not fall back silently to a smaller panel. | PASS_BLOCKED | present:sha256=37888988af5fdf6148d701477c388b4f07857d97565f0db78d8257e7adff6b49 |
| Local multiplicity or h2-gate leakage | Freeze the pair-by-locus family, LAVA local-h2 Bonferroni gate, and local-rg BH family before result access. | PASS_CONTRACT | config/local_architecture_contract.json |
| Replication candidate attrition | Lock the entire Tier A/B candidate family before source curation and preserve NO_INDEPENDENT_DATASET outcomes. | PASS_CONTRACT | scripts/21_prepare_replication_queue.py+22_lock_replication_manifest.py+23_collate_replication.py |
| LAVA simulation instability | Freeze a deterministic nonzero simulation seed per selected pair and retain it in the manifest. | PASS_CONTRACT | config/local_architecture_contract.json |
| Fine-mapping LD or allele mismatch | Require identical dense SNP order/alleles, checksum-locked signed LD, PSD/symmetry checks, RSS-LD s, and kriging outlier diagnostics. | PASS_CONTRACT_INPUTS_BLOCKED | config/fine_mapping_colocalization_contract.json |
| Molecular-QTL source attrition | Complete eQTL/sQTL/pQTL searches per locked locus and preserve evidence-backed NO_SUITABLE_DATASET rows. | PASS_CONTRACT_INPUTS_BLOCKED | scripts/33_prepare_finemapping_queue.py+34_lock_finemapping_manifest.py |
| Coloc prior cherry-picking | Lock p1/p2/p12 and the full p12 sensitivity grid before result access; label prior robustness separately. | PASS_CONTRACT_INPUTS_BLOCKED | config/fine_mapping_colocalization_contract.json |
| Mechanistic source or release drift | Require exact release/accession, local query/result snapshot, SHA-256, access date, license/terms, and primary citation for every supported evidence row. | PASS_CONTRACT_INPUTS_BLOCKED | config/mechanistic_annotation_contract.json+config/mechanistic_sources.tsv |
| Mechanistic chain gap concealment | Emit every required chain edge, mark unsupported edges MISSING, and forbid a complete causal narrative from annotation-only evidence. | PASS_CONTRACT_INPUTS_BLOCKED | scripts/39_validate_mechanistic_evidence.py+40_synthesize_mechanisms.py |

## Findings-level questions that remain blocked

- Heritability: do rerun h2 Z scores, intercepts, liability conventions, effective sample sizes, and phenotype definitions support each retained extension trait?
- Global correlation: are cross-trait intercepts, valid-allele SNP overlaps, direction stability, duplicate phenotypes, and the isolated 1,200-test family acceptable?
- Novelty: does every FDR hit survive direct pair-by-pair literature searches, publication-duplicate checks, and same-phenotype/same-direction/same-sleep-context review?
- Replication: are cohorts truly non-overlapping, ancestry/build/effect definitions compatible, effects concordant, and heterogeneity acceptable?
- Local and pleiotropic analyses: are source-verified LD references compatible, local h2 estimable, simulation seeds locked, and horizontal/vertical/shared-factor/sample-overlap alternatives still plausible?
- Fine-mapping and colocalization: are dense alleles and signed LD exactly aligned, LD-summary diagnostics clean, credible sets pure, models converged, and H0-H4 conclusions robust over the locked priors?
- Mechanism: does every proposed edge have an exact release/accession, local snapshot, checksum, primary citation, compatible tissue/cell context, and an explicit noncausal claim boundary?

## Current adversarial conclusion

The extension design has explicit controls for the main foreseeable failure modes, including pre-result panel and family locks, unavailable-dataset retention, ancestry/build/effect-allele checks, independent replication, local LD validation, prior sensitivity, and edge-level mechanistic provenance. Those controls have synthetic execution evidence only. They cannot validate the absent real analyses.

Publication-grade discovery claims remain disallowed until acquisition and every downstream findings-level challenge above are completed without weakening thresholds or changing the locked families.
