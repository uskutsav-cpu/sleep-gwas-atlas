# Exploratory downstream method readiness — 2026-09-27

All seven mounted dense GWAS files match the locked byte sizes and SHA-256 values.
The source cards retain different `N` semantics (`total` or `effective`).
This does not itself validate a trait-specific SuSiE/coloc sample-size model.

| Method | Class | Current state | Evidence/limitation |
|---|:---:|---|---|
| five_track_candidate_reconciliation | A | COMPLETE_EXPLORATORY | All 25 pair candidates and 20 geographic groups are receipt bound. |
| independent_ld_and_signed_direction | A/B | COMPLETE_EXPLORATORY_WITH_EXCEPTIONS | Five cross-pair genotype LD edges and 2,686 signed-LD candidate statuses are versioned; six raw R values fail range QC. |
| locus_specific_replication | A | NOT_ESTABLISHED | Existing FinnGen ADHD result is global pair-level directional replication only. |
| susie_fine_mapping | B | NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED | Pinned R lacks susieR; reviewed trait-specific N/case-fraction and full regional LD contract absent; UKB raw LD has documented out-of-range entries. |
| trait_trait_coloc | B | NOT_RUN_METHOD_PREREQUISITES_UNREVIEWED | Pinned R lacks coloc; reviewed regional dataset, sample-size, prior and overlap settings absent; multi-signal path also requires valid fine-mapping. |
| eqtl_sqtl_coloc | B | NOT_RUN_NO_VERIFIED_QTL_INPUT | No source-verified full-cis QTL dataset or reviewed molecular follow-up configuration is locally bound. |
| regulatory_cell_type_pathway | B | DESCRIPTIVE_COORDINATE_CONTEXT_ONLY | B Ensembl lead-level regulatory lookup was cached; no weighted functional input, tissue panel/background or multiplicity family is reviewed. |
| lava_local_rg_and_final_tiers | C | BLOCKED_LAVA | Canonical v3 failed frozen family QC: 3,720 NOT_RUN versus maximum 873. |

Pinned native runtime: `susieR=FALSE coloc=FALSE`; PLINK: `None`; PLINK2: `None`.
Locally bound QTL directories with files: 0.
Missing prerequisites are method blocks, not evidence of a negative biological result.
No protected LAVA decision or frozen promotion rule was modified.
