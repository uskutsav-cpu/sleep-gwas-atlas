# Reviewer 2 audit of current sleep–frailty evidence

Audit date: 2026-09-26 21:26 UTC. This audit covers four distinct evidence sets now
reported in the provisional manuscript: (1) the 12 frozen-atlas sleep × FI
correlations, (2) the completed secondary 12 × 7 sleep × latent-frailty
correlation family, (3) a 72-pair read-only aging-context extract from the
same frozen atlas, and (4) the locked FI × 12-sleep LAVA sensitivity
family, active and incomplete at the 21:25 UTC snapshot (17,887/29,940 receipts validated; inventory advanced by five during scanning; runner state reached 17,891 and post-scan inventory 17,892, leaving 12,048 slots). The locked 1% failure lower
bound is already exceeded for all 12 traits. It does not imply that Fried,
HFRS, component, family-level local,
fine-mapping, colocalization, or molecular analyses are complete. The first
15 questions below remain specifically adjudicated for the FI result; the
following section addresses the latent-factor and aging-context evidence.
The FI table is `analysis/frozen_atlas_frailty_global_rg.tsv`, copied from
the frozen atlas with its all-396-pair q-values unchanged. The latent table is
`results/frailty_v1/rg_sleep_latent_factors.tsv`, with its fixed 84-pair BH
family and run manifest. The aging table is
`analysis/frozen_atlas_aging_context_global_rg.tsv`, retaining original
all-396-pair q-values. Input/output hashes and reproduction limits are
recorded in the corresponding manifests and provenance files.

## Reviewer questions and adjudication

| Reviewer question | Evidence and current verdict | Claim consequence / required follow-up |
|---|---|---|
| 1. Is the phenotype definition appropriate? | **Partly supported.** The primary endpoint is the registered European-descent FI GWAS, GCST90020053, described locally as a UK Biobank/TwinGene meta-analysis. It is the prespecified central endpoint. The distinct Fried phenotype is now represented by an acquired full file and is not interchangeable with FI; exact-release build/effect conventions and participant intersections remain unresolved, so it has not been analyzed. The exact operational FI item set and source-level phenotype construction are not reproduced from raw inputs in this checkout. | Call this an FI association only. Do not generalize to all frailty constructs or claim construct replication. Verify the exact Fried source metadata and FI phenotype definition before harmonization or manuscript-level phenotype comparison. |
| 2. Could sample overlap explain it? | **Yes, potentially; unresolved.** Cohort-level metadata indicate expected UK Biobank overlap between the FI and 11 of 12 sleep sources. Exact participant intersections are unknown. Sleep apnea is registered to FinnGen R9, with no shared cohort identified, but this does not prove zero overlap. | Do not call these independent replications. Preserve LDSC cross-trait intercept/overlap diagnostics where available; exact overlap sensitivity remains open. |
| 3. Is heritability adequate? | **Pass; h² independently rerun on the current versioned FI input.** The frozen summary reports h²=0.1093 (SE=0.0050; Z=21.86), intercept=1.02, ratio=0.0509, and PASS. A separate rerun from the current versioned munged FI input reproduced h²=0.1093 (SE=0.0050; Z=21.86), intercept=1.0203 (SE=0.0086), ratio=0.0509, and 1,167,050 SNPs after reference merges. | This corroborates the reported h² gate at displayed precision. It is not an exact end-to-end replay of atlas-v1.0: its legacy h² summary and exact harmonized-input hash remain unavailable, and the current versioned harmonized FI file has a different hash. It also does not establish adequate power for every sleep phenotype/local test. See `analysis/frailty_index_h2_reproduction_manifest_2026-09-23.json`. |
| 4. Does it replicate? | **No frailty-specific independent replication established.** The frozen atlas contains one FI source and no separate FI replication. Fried-score statistics are acquired, but exact-release build/effect conventions are unresolved; the locked plan treats Fried as construct replication and UK Biobank overlap is expected. | Do not use “replicated” or “independently confirmed.” Reassess only after phenotype direction, ancestry, build, effect conventions, and participant overlap are documented and a separately justified analysis is completed. |
| 5. Is the direction concordant? | **Descriptive only; cross-source concordance unavailable.** The frozen family has positive estimates for insomnia, short sleep, long sleep, sleepiness, napping, snoring, sleep apnea and sleep timing; negative estimates for sleep duration, chronotype, sleep efficiency and accelerometer duration. Three of these 12 do not pass all-396 FDR. No independent frailty result is available for direction comparison. | Report signs only with trait coding and source definitions. Do not call directions replicated or biologically concordant across datasets. |
| 6. Is the result robust to alternate frailty definitions? | **Partly explored, not independently validated.** The latent-factor family provides distinct dimensional phenotypes, but it is sensitivity-only because exact participant overlap is unknown. No eligible, verified Fried-score or component GWAS comparison has been run. | Limit FI conclusions to the registered FI. Describe latent factors separately as sensitivity analyses; do not call them construct replication or claim component-level robustness. Fried and component comparisons remain future work in separate, source-validated families. |
| 7. Is the local signal credible? | **No local inference is admissible.** The D13 sensitivity-only 12-pair × 2,495-region family remains incomplete. The 2026-09-26 21:24:40–21:25:11 UTC full audit validated 17,887 receipts and found all 12 trait-level failure lower bounds exceed the frozen 1% allowance. Exact participant intersections remain unknown, and the overlap parameter is LDSC-intercept-derived. | Do not report or interpret partial regional estimates. Continue receipt completion and final collation for provenance, applying unchanged gates; current failures already prohibit local inference. These results establish neither local sharing nor its absence and cannot establish independent replication. |
| 8. Is the shared locus supported by more than one method? | **Not assessed.** No frailty-specific locus evidence exists. | Make no shared-locus or convergent-method claim. A blocked method must not be described as a negative result. |
| 9. Does colocalization support H4 rather than H3? | **Not assessed.** No frailty-specific coloc analysis has run; no eligible loci or matched inputs exist. | Make no colocalization or shared-causal-signal claim. If reached, report PP.H4 and PP.H3 under the frozen priors and diagnostics. |
| 10. Is the gene assignment actually supported? | **Not applicable to current results.** Global rg does not nominate genes; there are no frailty-specific fine-mapped loci. | Do not name candidate causal genes from these correlations. |
| 11. Is the tissue/cell type biologically plausible and statistically supported? | **Not assessed.** Tissue/cell resource indexes are contextual inventories only; no frailty lead loci have been frozen and no locus-matched QTL or cell analysis exists. | Do not make tissue, cell-type, or mechanism claims from resource availability or biological plausibility alone. |
| 12. Are claims stronger than the evidence? | **High risk unless carefully bounded.** The support is for global genetic correlation in a frozen atlas, not causality, mediation, local sharing, shared variants, or clinical prediction. Full source-input reprocessing is unavailable in this checkout. | Use “global genetic correlation” or “genome-wide genetic correlation” for the supported results. Avoid the broader phrase “shared genetic architecture,” which could imply locus-level or functional convergence not demonstrated here. State the frozen-output/source-reprocessing limitation and cohort overlap. Avoid causal or mechanistic language. |
| 13. Could population structure or LD differences explain it? | **Partly mitigated, not eliminated.** The FI row is European ancestry and the plan requires ancestry-matched LD/build. The frozen results are global LDSC estimates; exact source inputs and ancestry-specific sensitivity cannot be reprocessed here. | Keep claims to the reported ancestry/source context. Do not generalize across ancestries. Recheck ancestry, build, allele alignment, and LD reference if source files become available. |
| 14. Is the finding already known? | **Yes, multiple prior genetic analyses exist.** In addition to Song and Deng et al. (2024), the frozen PubMed genetic query contains Li et al. (2024), which reports FI genetic correlations for sleep duration and napping; Lu et al. (2024), Che et al. (2025), Zhang et al. (2025), and Gao et al. (2025) report MR analyses across sleep and FI/Fried frailty; Zhang et al. (2026) includes FI in a multidimensional sleep-aging MR; Huang et al. (2026) reports OSA–frailty MR and is present in the frozen PubMed snapshot as one deduplicated queue record from three query-source occurrences. A 2026 CLSA ordinal Fried-phenotype GWAS adds a small phenotype-side discovery with one lead association, but it does not test sleep or provide reported effect estimates, h², or a validation sample (`analysis/clsa_ordinal_frailty_gwas_resource_audit_2026-09-23.md`). MR does not replicate LDSC rg. Che et al.’s FI and sleep-duration meta-analyses both include UKB participants; exact intersections are unknown. Source overlap for other papers remains incompletely reconciled. Queue records remain unscreened. | Do not claim first sleep–frailty genetic evidence, first sleep-duration/napping global rg with FI, or first sleep–frailty MR. Frame any contribution as phenotype-specific expansion and only claim new estimates after source/QC and independent-replication gates pass. Complete the registered screening. |
| 15. Is there a simpler explanation? | **Plausible alternatives remain.** Shared UK Biobank participants, broad/self-reported sleep phenotypes, FI construct differences, correlated health and behavioral liabilities, and residual population structure could contribute to genome-wide correlations. These alternatives are not separated by the current results. | Interpret results as noncausal genome-wide covariance. Avoid mechanism, direction-of-effect, or intervention conclusions. |

## Retraction follow-up (2026-09-24)

The frozen PubMed queue contains a 2024 sleep–frailty Comment (PMID 38372839) and its 2026 retraction notice (PMID 41811627). The publisher marks the Comment retracted on 11 March 2026. Neither record has been screened; the Comment is not treated as empirical evidence. Its separately cited CKD cohort article (PMID 38289547) is also in the frozen queue and remains unscreened. No screening or PRISMA counts were changed. The retraction notice's reason was not available in the publisher preview and is not inferred. See `analysis/retracted_sleep_frailty_comment_audit_2026-09-24.md`.

## Additional adjudication: latent factors and aging context

| Reviewer question | Evidence and current verdict | Claim consequence / required follow-up |
|---|---|---|
| 16. Are the latent-factor analyses complete and multiplicity controlled? | **Yes, for the frozen secondary family.** The table contains 84 unique sleep × factor pairs; every row records `family_denominator=84`, and 48 have BH q≤0.05. The run manifest records input and log checksums and frozen family/execution hashes. | Report estimates and q-values from the frozen table. Do not recalculate q-values on selected factors or traits. This completion does not make the secondary family primary. |
| 17. Do the latent factors establish independent replication or resolve overlap? | **No.** The factor GWAS contain UK Biobank data; exact participant overlap with the sleep GWAS remains unknown. All 84 rows are marked sensitivity-only. | Do not call any latent-factor association replication. Retain the cohort-overlap caveat for each interpretation. |
| 18. Is the insomnia–general-factor correlation independent of phenotype construction? | **No.** The Foote et al. model includes insomnia among its 30 input deficit GWAS, and all 30 load onto the general factor. The project’s Jansen insomnia GWAS is a different release, but it shares the UK Biobank cohort with the factor study; exact participant intersection is unknown. This is a construct-level part-whole dependent comparison. | State this directly beside the estimate. Do not present the large insomnia–general-factor rg as independent corroboration or replication. The daytime-sleepiness trait is also related, but not identical, to the model’s tiredness/lethargy indicator. |
| 19. Do QC gates support the latent-factor rg results? | **The recorded h² gates pass, with a boundary case.** All seven factor rows are `PASS`; Factor 3 has LDSC intercept 1.195 against the locked ceiling of 1.20. The results do not by themselves establish local-test power or remove possible overlap bias. | State the inherited gate and Factor 3 value; avoid implying broad robustness. Keep downstream regional analyses gated by their own criteria. |
| 20. Do differences among factors show which frailty dimension drives sleep–frailty overlap? | **No direct comparison was performed.** Factor-specific rg estimates and q-significance counts are descriptive; factors and estimates may be correlated, and the analysis does not test their differences. | Do not say a factor drives, explains more of, or is stronger than another. Use the forest plot only descriptively, with its distinct q-value families disclosed. |
| 21. Do the aging-context estimates independently confirm the frailty findings? | **No.** These are 72 read-only correlations already present in the frozen atlas; 26 retain q≤0.05 under its original all-396 correction. The source-level crosswalk classifies 44 pairs as shared-cohort-expected and 28 as possible/unknown; exact participant intersection is unknown for all 72. | Present only as contextual estimates, not replication or a newly run frailty family. Preserve original q-values and state the frozen-atlas reuse. |
| 22. Is a physical frailty component identified as the source of the associations? | **No.** No eligible, source-validated physical-component rg family is available. The acquired Fried file does not resolve component-specific results, and the exact HFRS endpoint remains unresolved. | Do not infer the component from latent-factor labels. The decomposition question remains open pending valid phenotype-specific inputs. |
| 23. Does the expanded evidence support loci, mechanisms, or causal direction? | **No.** The local LAVA family remains active and incomplete (17,887/29,940 receipts validated at 21:25:11 UTC); all 12 trait-level failure lower bounds exceed the locked 1% allowance. No frailty-specific PLACO, fine-mapping, trait–trait colocalization, or locus-matched QTL/cell-type result exists. | Keep claims at genome-wide correlation. Do not infer loci, genes, tissues, cellular pathways, mediation, or causality. Continue the frozen receipt run for provenance and apply the final family audit; current failures already preclude local inference. |

## Frozen result family: allowed statement and downgrades

The frozen FI table contains 12 global rg estimates; 9 have the original
all-396-family BH q-value at or below 0.05. The completed secondary latent
factor table contains 84 estimates; 48 have BH q≤0.05 within that locked
family. The read-only aging-context table contains 72 estimates; 26 have the
original all-396-family q≤0.05. These three counts come from different
families and do not support comparisons across families. The FI and aging
tables reuse frozen atlas estimates; neither is a newly run frailty analysis.
The latent-factor estimates are sensitivity-only, not independent
replication. Quote estimates and q-values from their source tables without
recalculating on selected subsets.

**Current claim strength: descriptive/exploratory for genome-wide correlation
patterns; unsupported for independent replication, component decomposition,
local sharing, pleiotropic loci, shared causal variants, molecular mechanism,
causal direction, and novelty.** The FI h² row passes the repository's
recorded inclusion gate as reported, but its original source inputs cannot be
fully reprocessed in this checkout. The latent-factor h² gates pass as
recorded, although exact overlap remains unresolved. The paper remains
provisional while source/QC, review-screening, and data-specific analysis
gates are incomplete.

## Evidence references

- Frozen estimates and original family q-values: `analysis/frozen_atlas_frailty_global_rg.tsv`.
- Latent-factor estimates, family denominator, and sensitivity-only labels: `results/frailty_v1/rg_sleep_latent_factors.tsv`; input/log checksums: `results/frailty_v1/rg_sleep_latent_factors_run_manifest.tsv`.
- Latent-factor h² gates: `results/frailty_v1/h2_latent_factors.tsv`.
- Frozen aging-context estimates and cohort crosswalk: `analysis/frozen_atlas_aging_context_global_rg.tsv` and `analysis/frozen_atlas_aging_context_overlap.tsv`.
- Frozen atlas hash and source-input limitation: `analysis/frozen_atlas_frailty_audit.md` and `manifests/frozen_atlas_frailty_manifest.json`.
- Cohort-level overlap and exact-overlap unknowns: `analysis/sample_overlap_assessment.tsv`.
- Locked phenotype hierarchy and downstream gates: `config/analysis_plan_v1.yaml` and `config/analysis_plan_decisions.tsv`; LAVA sensitivity-family lock and runtime details: `config/lava_frailty_sensitivity_v1.yaml` (D13).
- Current LAVA full-family receipt audit: `analysis/lava_full_receipt_integrity_2026-09-26_2125_checkpoint.json`; per-trait gate summary: `analysis/lava_full_receipt_failure_gate_audit_2026-09-26_2125_checkpoint.md`. The run remains incomplete, all 12 locked trait gates fail, and no local inference is permitted. Earlier coordinator and interim checkpoints remain historical.
- Review progress and manual-search limitation: `review/README.md`, `review/reporting_checklist_status.tsv`, and `review/manual_search_instructions.md`.
