# Track B mechanism follow-up: pre-result local readiness — revision 3

Revision 3 supersedes `results/track_b/mechanism_followup/v2/PRE_RESULT_READINESS.lock.json` (SHA-256 `6643561657f69d13e6e5c8e0352d9c095a401943e4028585562f81e3de3ac530`). Revisions 1 and 2 remain byte-exact audit trails and must not be used for execution. Revision 3 is the authoritative pre-result readiness contract.

This contract used only local resources and immutable upstream metadata. It did not read fine-mapping, PLACO/pleiotropy, or future mechanism results. Missing, blocked, skipped and unavailable resources are never null findings.

- Policy SHA-256: `30f6cd67b1227421664f55b3cb20f5c2139f9f0293365d687c42fa03309de0a4`
- Complete analysis family: 31 records
- Exact evidence family: 818 filesystem records
- Resource readiness: BLOCKED_BY_DATA=18, BLOCKED_BY_SOFTWARE=6, READY=7
- Current pre-result state: BLOCKED_BY_DATA=18, BLOCKED_BY_SOFTWARE=6, NOT_APPLICABLE_UNTIL_UPSTREAM=7

## Readiness matrix

| Phase | Analysis ID | Resource status | Current status | Exact boundary |
|---:|---|---|---|---|
| 7 | P07_CONJFDR | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no checksum-pinned MATLAB-compatible conjFDR runtime/reference/executor manifest exists |
| 9 | P09_HDL_L | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no checksum-pinned ancestry/build-matched HDL-L reference manifest exists |
| 9 | P09_RHO_HESS | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no separate sealed Track B rho-HESS reference/contract exists |
| 13 | P13_MAGMA_BROAD_TISSUE | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no complete human broad-tissue MAGMA expression matrix is local; FUMA scRNA is not a broad-tissue substitute / no checksum-pinned complete broad-tissue MAGMA adapter is local |
| 13 | P13_SLDSC_GTEX | BLOCKED_BY_DATA | BLOCKED_BY_DATA | the exact local S-LDSC family contains 16 selected GTEx tissues plus control; the required chromosome-wide Lung index-36 and Muscle_Skeletal index-38 annotations and an 18-tissue checksum manifest are absent / the baselineLD v2.2/weights/HapMap3/static EUR reference provenance is absent |
| 13 | P13_GTEX_TISSUE_EQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | all 49 GTEx v8 ge payloads and their 49 tabix indexes are absent / no repository-local checksum-pinned tabix entrypoint exists |
| 13 | P13_GTEX_TISSUE_SQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no preregistered GTEx v8 tissue-sQTL manifest/payload/index family is local / no repository-local checksum-pinned tabix entrypoint exists |
| 14 | P14_BRAIN_CELL_CLASS_DISCOVERY | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no frozen Allen-MTG 75-subtype-to-seven-class mapping/aggregation executor is local |
| 14 | P14_HUMAN_BRAIN_PERICYTE_DISCOVERY | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no checksum-pinned human brain expression atlas containing pericytes/vascular-associated cells is local |
| 15 | P15_BRAIN_CELL_SUBTYPE_DISCOVERY | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 16 | P16_HUMAN_CELL_CLASS_REPLICATION | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no frozen Allen-MTG 75-subtype-to-seven-class mapping/aggregation executor is local / no checksum-pinned cross-atlas class/subtype alignment and replication-classification executor is local |
| 16 | P16_HUMAN_CELL_SUBTYPE_REPLICATION | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no independent human subtype-resolved replication atlas is local; GSE67835 has only seven broad classes and mouse is ineligible / no checksum-pinned cross-atlas class/subtype alignment and replication-classification executor is local |
| 17 | P17_CELL_EQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000559/QTD000569 eQTL payloads and matching .tbi indexes are absent / no repository-local checksum-pinned tabix entrypoint exists |
| 18 | P18_CELL_SQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000563/QTD000573 leafcutter sQTL payloads and matching .tbi indexes are absent / no repository-local checksum-pinned tabix entrypoint exists |
| 19 | P19_THREE_WAY_EQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000559/QTD000569 eQTL payloads and matching .tbi indexes are absent / no repository-local checksum-pinned tabix entrypoint exists / no checksum-pinned multi-signal three-way coloc executor exists; pairwise coloc cannot substitute |
| 19 | P19_THREE_WAY_SQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000563/QTD000573 leafcutter sQTL payloads and matching .tbi indexes are absent / no repository-local checksum-pinned tabix entrypoint exists / no checksum-pinned multi-signal three-way coloc executor exists; pairwise coloc cannot substitute |
| 20 | P20_CATLAS_SCATAC | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 20 | P20_SCREEN_CHROMATIN | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 20 | P20_HOCOMOCO_MOTIF | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no local checksum-pinned hg38 FASTA/sequence family exists; network sequence calls are forbidden for this local contract |
| 21 | P21_ABC_ENHANCER_GENE | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 21 | P21_PCHIC_ENHANCER_GENE | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 21 | P21_TWAS | BLOCKED_BY_DATA | BLOCKED_BY_DATA | locked result-blind inventory is valid (98 files/49 contexts/3135665776 bytes) but model payloads are absent / MetaXcan source/entrypoint are pinned, but the required exact .molecular-env/python runtime is absent or unvalidated |
| 22 | P22_REGULATORY_CHAIN_EQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000559/QTD000569 eQTL payloads and matching .tbi indexes are absent / no checksum-pinned same-locus/same-context regulatory-chain integrator exists |
| 22 | P22_REGULATORY_CHAIN_SQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | QTD000563/QTD000573 leafcutter sQTL payloads and matching .tbi indexes are absent / no checksum-pinned same-locus/same-context regulatory-chain integrator exists |
| 23 | P23_SPATIAL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no authorized local human spatial payload/metadata manifest exists |
| 24 | P24_LOCKED_ATLAS_CROSS_SLEEP | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no checksum-pinned complete 12-trait cross-sleep integration executor exists |
| 25 | P25_PAIR_A_VS_B | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no checksum-pinned deterministic Pair-A-versus-Pair-B comparison executor is local |
| 26 | P26_POSITIVE_CONTROL | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 27 | P27_PATHWAYS | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 28 | P28_GRN | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no matched human expression-accessibility dataset/metadata manifest exists for GRN inference / no validated executor matched to a frozen multimodal GRN design exists |
| 29 | P29_BIDIRECTIONAL_MR | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | MR packages and EUR reference are locally valid, but no validated automatic Track B bidirectional MR science executor exists |

## Corrected scientific boundaries

- GTEx tissue eQTL, GTEx tissue sQTL, cell eQTL and leafcutter cell sQTL are four distinct families. GTEx eQTL metadata has 49 frozen `ge` contexts, but the 49 payloads and indexes are absent. No GTEx tissue-sQTL family is locally preregistered.
- Cell eQTL evidence is limited to QTD000559 microglia and QTD000569 neuron. Cell sQTL evidence is limited to QTD000563 microglia and QTD000573 neuron leafcutter data. Payloads/indexes and tabix are absent. PsychENCODE remains controlled supplementary access, not human replication or local evidence.
- Allen Human MTG contains 75 subtype columns spanning seven named class prefixes, but no frozen subtype-to-broad-class mapping/aggregation executor exists. Broad-class discovery is therefore software-blocked; the full subtype unit remains resource-ready but awaits its upstream class gate. Neither Allen MTG nor GSE67835 contains pericytes, so the human pericyte unit is data-blocked.
- GSE67835 is an independent human seven-class resource, but no frozen cross-atlas class alignment/replication classifier exists; two unadjusted enrichment tables are not replication. The two vascular FUMA matrices are mouse mapped to human gene identifiers and cannot be used as human replication. No airway or vascular human expression reference was invented.
- The complete 30-file interpretation code family, including task preparation, runners, recorder and family-wide BH collator, is checksum-pinned. S-LDSC has exactly its locally frozen 16 selected GTEx tissues plus the shared control (375 files) and exact runtime; this is not the complete intended 18-tissue family. Chromosome-wide GTEx Lung index 36 and Muscle_Skeletal index 38 annotations are absent, and the static EUR baselineLD/weights/HapMap3/genotype reference is also absent. Lung must not be described as an airway-specific reference. Broad-tissue MAGMA lacks both its complete tissue matrix and a dedicated adapter.
- CATlas adult human scATAC, SCREEN, ABC, immune-only PCHiC, HOCOMOCO motifs and four pathway resources validate locally. HOCOMOCO lacks a pinned hg38 sequence family. The Pair-A-versus-Pair-B executor is absent.
- The 98-file/49-context TWAS inventory is frozen, but model payloads and the exact MetaXcan runtime are absent. Spatial data, a matched human multimodal GRN resource/runtime, HDL-L, rho-HESS and conjFDR remain blocked.
- MR package/reference resources validate, but the automatic Track B bidirectional executor is absent. A causal claim remains prohibited unless every predeclared robustness gate passes.

## RAM-equivalent decomposition

Locus-dependent QTL, colocalization, chromatin linking and local-correlation analyses may run one complete locus per fresh process; an LD-dependent locus must never be split. Trait, tissue, cell and pathway work may run one frozen unit at a time while preserving the complete preregistered family and correction denominator. Global nuisance estimation must remain genome-wide; chunking is allowed only after global parameters are frozen and mathematical identity is demonstrated. GRN batching remains blocked until method-specific equivalence is proved.

## Claim boundary

`READY` means only that local resources/software passed this result-blind contract. Every such unit remains `NOT_APPLICABLE_UNTIL_UPSTREAM` until its scientific gate is satisfied. A blocker, unavailable resource, skipped unit or failure must be preserved and must never be reported as a negative biological result.
