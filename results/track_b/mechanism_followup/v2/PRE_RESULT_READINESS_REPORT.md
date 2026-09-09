# Track B mechanism follow-up: pre-result local readiness — revision 2

This contract supersedes `results/track_b/mechanism_followup/PRE_RESULT_READINESS.lock.json` (SHA-256 `4ab3c004d109fa856936a5411f4495975b176333916bbe3f70b7974aac9de40f`). Revision 1 is retained byte-exact as an audit trail and must not be used for execution. Revision 2 corrects the cell-sQTL family to the predeclared leafcutter datasets QTD000563 and QTD000573 only.

This is a result-blind resource/software contract. It did not read fine-mapping, PLACO/pleiotropy, or future mechanism results and it makes no post-result scientific claim. Missing, blocked, skipped, and unavailable resources are never interpreted as null evidence.

- Policy SHA-256: `d6adcd5ec9863ddef2cd7c36d6fbff95cd243ebe973bb0e0fa726c1817e70e81`
- Complete analysis family: 27 records
- Exact evidence family: 696 filesystem records
- Resource readiness: BLOCKED_BY_DATA=14, BLOCKED_BY_SOFTWARE=3, READY=10
- Current pre-result state: BLOCKED_BY_DATA=14, BLOCKED_BY_SOFTWARE=3, NOT_APPLICABLE_UNTIL_UPSTREAM=10

## Readiness matrix

| Phase | Analysis ID | Resource status | Current status | Exact boundary |
|---:|---|---|---|---|
| 7 | P07_CONJFDR | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no checksum-pinned MATLAB-compatible conjFDR runtime/reference/executor manifest exists |
| 9 | P09_HDL_L | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no checksum-pinned ancestry/build-matched HDL-L reference manifest exists |
| 9 | P09_RHO_HESS | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no separate sealed Track B rho-HESS reference/contract exists |
| 13 | P13_MAGMA_BROAD_TISSUE | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no complete human broad-tissue MAGMA expression matrix is local; FUMA scRNA is not a broad-tissue substitute |
| 13 | P13_SLDSC_GTEX | BLOCKED_BY_DATA | BLOCKED_BY_DATA | the baselineLD v2.2/weights/HapMap3/static EUR reference provenance is absent |
| 13 | P13_GTEX_TISSUE_QTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | complete selected QTL payload plus .tbi family is absent / no repository-local checksum-pinned tabix entrypoint exists |
| 14 | P14_BRAIN_CELL_CLASS_DISCOVERY | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 15 | P15_BRAIN_CELL_SUBTYPE_DISCOVERY | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 16 | P16_HUMAN_CELL_CLASS_REPLICATION | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 16 | P16_HUMAN_CELL_SUBTYPE_REPLICATION | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no independent human subtype-resolved replication atlas is local; GSE67835 has only seven broad classes and mouse is ineligible |
| 17 | P17_CELL_EQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | complete selected QTL payload plus .tbi family is absent / no repository-local checksum-pinned tabix entrypoint exists |
| 18 | P18_CELL_SQTL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | complete selected QTL payload plus .tbi family is absent / no repository-local checksum-pinned tabix entrypoint exists |
| 19 | P19_THREE_WAY_COLOC | BLOCKED_BY_DATA | BLOCKED_BY_DATA | complete selected QTL payload plus .tbi family is absent / no repository-local checksum-pinned tabix entrypoint exists / no checksum-pinned multi-signal three-way coloc executor exists; pairwise coloc cannot substitute |
| 20 | P20_CATLAS_SCATAC | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 20 | P20_SCREEN_CHROMATIN | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 20 | P20_HOCOMOCO_MOTIF | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no local checksum-pinned hg38 FASTA/sequence family exists; network sequence calls are forbidden for this local contract |
| 21 | P21_ABC_ENHANCER_GENE | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 21 | P21_PCHIC_ENHANCER_GENE | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 21 | P21_TWAS | BLOCKED_BY_DATA | BLOCKED_BY_DATA | locked result-blind inventory is valid (98 files/49 contexts/3135665776 bytes) but model payloads are absent / MetaXcan source/entrypoint are pinned, but the required exact .molecular-env/python runtime is absent or unvalidated |
| 22 | P22_REGULATORY_CHAIN | BLOCKED_BY_DATA | BLOCKED_BY_DATA | complete selected QTL payload plus .tbi family is absent / no checksum-pinned same-locus/same-context regulatory-chain integrator exists |
| 23 | P23_SPATIAL | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no authorized local human spatial payload/metadata manifest exists |
| 24 | P24_LOCKED_ATLAS_CROSS_SLEEP | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | no checksum-pinned complete 12-trait cross-sleep integration executor exists |
| 25 | P25_PAIR_A_VS_B | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 26 | P26_POSITIVE_CONTROL | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 27 | P27_PATHWAYS | READY | NOT_APPLICABLE_UNTIL_UPSTREAM | local resource/software contract is READY; scientific execution awaits its predeclared upstream gate |
| 28 | P28_GRN | BLOCKED_BY_DATA | BLOCKED_BY_DATA | no matched human expression-accessibility dataset/metadata manifest exists for GRN inference / no validated executor matched to a frozen multimodal GRN design exists |
| 29 | P29_BIDIRECTIONAL_MR | BLOCKED_BY_SOFTWARE | BLOCKED_BY_SOFTWARE | MR packages and EUR reference are locally valid, but no validated automatic Track B bidirectional MR science executor exists |

## Scientifically important boundaries

- Human brain discovery is available in Allen Human MTG level 2; independent human replication is available only at seven broad classes in GSE67835. No independent human subtype atlas is local. The two vascular FUMA matrices are mouse mapped to human gene IDs and cannot count as human replication.
- The frozen S-LDSC annotation source family contains 16 GTEx tissues in four domains, but its static genotype/baselineLD/weights/HapMap3 reference family is absent. Lung/airway and skeletal-muscle coverage are not present in that 16-tissue selection.
- eQTL Catalogue and GTEx metadata are present, but the molecular-QTL payloads and matching tabix indexes are not. PsychENCODE remains controlled supplementary access, not a local resource.
- SCREEN, CATlas, ABC and immune-only promoter-capture Hi-C sources are local. HOCOMOCO motifs are local but no immutable local hg38 sequence family is available. No airway-specific or vascular human expression reference was invented.
- The 98-file/49-context phi-enabled TWAS inventory is locked before results, but all model payloads and the exact MetaXcan Python runtime are absent.
- No authorized spatial dataset or matched multimodal human GRN dataset/runtime exists locally. Those absences are blockers, never negative findings.
- MR packages and the 503-sample European GRCh37 LD reference validate, but a production executor enforcing the full bidirectional estimator/overlap/Steiger/colocalization/FDR contract is absent.
- HDL-L, rho-HESS and conjFDR remain blocked by their method-specific runtime/reference contracts.

## RAM-equivalent decomposition

Locus-dependent QTL, colocalization, regulatory linking and local-correlation work may run one complete locus per fresh process; an LD-dependent locus must never be split. Trait/cell/pathway analyses may run one frozen trait/dataset/resource unit at a time while preserving every predeclared family and correction denominator. Global nuisance estimation (for example conjFDR or any method that estimates genome-wide parameters) must remain global; chunked tests are permitted only after those parameters are frozen and mathematical identity is demonstrated.

## Claim boundary

A READY resource status means only that the local inputs and entrypoints passed this pre-result contract. It is not authorization to bypass an upstream gate, not evidence that an association exists, and not a biological or causal conclusion.
