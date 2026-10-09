# Adversarial statistical-validity review v1

Client audit date: 2026-10-08. Reviewer scope: independent read-only scientific audit requested by Phase 10. No manuscript prose, source modification, historical result modification, native LDSC run, or test-suite-based scientific validation was performed by this reviewer.

**Verdict: frozen arithmetic and available printed-log concordance pass; fully independent replication and calibrated discovery/validation heterogeneity do not pass. The archived multivariate correlation-uncertainty export has an additional actionable estimand-labeling error.**

## Executed independent checks

`statistical_validity_v1.py` uses Python standard-library CSV, SHA-256, normal-tail arithmetic and an independently written BH suffix-minimum calculation. It imports no repository analysis implementation. Its receipt is `statistical_validity_v1.json`; `statistical_validity_v1.tsv` contains 89 compact numerical diagnostics with exact input paths and one-based file lines. The receipt hashes every numerical/log input read. The script runs from any working directory and writes only these reviewer outputs.

The archival root used here is `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas`. Core tables were read from `sleep_unified_research_v1/sources/recovered/results/tables/`; extension and replication tables were read from `discovery_extension/results/`.

| Check | Independently verified result | Boundary |
|---|---|---|
| Core Cartesian family | 12 × 33 = 396 unique rows, no omissions or duplicates | Family retained separately from extension |
| Core BH, all 396 original P values | Adjusted values match exactly; 161 numerical positives = 153 primary + 8 QC-failed sensitivity | 153 is the primary scientific count under the frozen all-396 correction |
| Original standalone h2/inclusion | 45 rows; 43 pass; T2D fails intercept; melanoma fails h2 Z | No classification mismatch with original logs |
| Core printed logs | All 396 rg/SE/Z/P/intercept records and 45 standalone h2 records match | Historical four-significant-digit serialization reproduced, not full unrounded estimates |
| Extension family and BH | 12 × 100 = 1,200; exact BH agreement; 603 positives | Discovery h2/rg logs absent at paths named in these extension tables during this audit |
| Extension recorded QC | 100/100 standalone and 1,200/1,200 pairwise outcome h2/intercept gates pass | Table arithmetic only for extension discovery |
| Replication family | 217 unique candidates, original manifest order retained; alpha = 0.05/217 = 0.0002304147465437788 | Denominator includes unavailable and QC-excluded candidates |
| Replication available estimates | 41 records match the eight retained rg logs; all 13 source h2 records match their logs | Printed-log concordance, not native rerun |
| Replication statistical classes | 23 positive, 18 concordant below threshold, 17 QC excluded, 159 no qualifying source | Zero two-trait-independent positives in the audited historical design |
| Effect-difference arithmetic | 41/41 historical zero-covariance Z/Q/P calculations agree | Covariance-zero assumption remains unverified |

For both global families, replacing original P values with normal tails of printed Z or printed rg/SE changes no BH labels. These substitutions are rounding diagnostics and do not recover omitted precision. There are no zero P values or nonpositive SE values in either frozen global table. All 18 recovered core artifacts named in the original checkpoint match their expected SHA-256, including the stale report discussed below. Exact recovery is distinct from scientific reproduction.

## Findings requiring action

### S1 — Multivariate export labels conditional standardized-covariance SE as rg SE

Priority: high for any proposed covariance-aware rg contrast or measurement-group experiment. The original 396 Python LDSC family is unaffected by this finding.

`scripts/25_genomicsem_covariance.R:154–170` exports `sqrt(diag(V_Stand))` as `genetic_correlation_se`, then computes correlation Z/P from it. In the pinned [GenomicSEM ldsc implementation](https://raw.githubusercontent.com/GenomicSEM/GenomicSEM/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R), lines 477–489 create `V_Stand` by rescaling covariance uncertainty using estimated diagonals as fixed constants. This does not propagate uncertainty in the estimated h2 denominators.

Independent arithmetic on archival `results/tables/ldsc_covariance_pairs.tsv` and `ldsc_sampling_covariance_1035x1035.tsv.gz` establishes:

- For all 1,035 entries, exported SE equals `SE(cov_ij)/sqrt(h2_i*h2_j)` to maximum relative error 6.33×10⁻¹⁵.
- Exported correlation Z equals covariance Z to maximum absolute error 1.03×10⁻¹⁴.
- All 45 diagonal correlations are exactly 1 but have positive exported correlation SE. A ratio-standardized self-correlation is identically 1; these positive SE values quantify uncertainty of the scaled variance parameter instead.

For off-diagonal `r = c/sqrt(a*b)`, the ratio delta-method gradient for `(a,c,b)` is `(-r/(2a), 1/sqrt(a*b), -r/(2b))`. Applying `gradient' V gradient` to the corresponding archived covariance elements changes the SE by more than 20% for 11 of 990 off-diagonal entries. Examples:

| Archived pair-table line | Pair | Exported fixed-scale SE | Ratio delta-method SE, conditional on archived V |
|---:|---|---:|---:|
| 48 | sleepdur–shortsleep | 0.0353847 | 0.00869510 |
| 372 | sleep_efficiency–accel_sleep_duration | 0.0541714 | 0.0256572 |
| 660 | ibd–crohn | 0.0745765 | 0.0132370 |
| 960 | longevity–parental_lifespan | 0.0669151 | 0.0897080 |

The diagnostic is not a calibrated replacement estimate: it inherits the archived 1,082-block V and its assumptions. Preserve the original export; label its uncertainty as fixed-scale standardized covariance. Before using rg SE/covariance for new inference, produce a separate justified ratio transformation with independent implementation agreement or aligned-block ratio jackknife, and freeze that experiment. Diagonal correlation P values should not be interpreted as tests of correlation being nonzero.

### S2 — Historical replication is qualified outcome-side validation

Priority: high for scientific claims. `discovery_extension/config/replication_manifest.tsv:2` explicitly records the same discovery sleep source for insomnia. The validation command at `discovery_extension/logs/replication/rg/rg_replication_insomnia.log:11` starts with `data/munged/insomnia.sumstats.gz`; the discovery script `14_rg_extension.sh` uses the same sleep-input convention. Equivalent source reuse is documented for all 23 positives in the prior evidence table `tables/replicated_23.tsv` under `sleep_GWAS_reused=True` and `two_trait_independence=NOT_FULLY_INDEPENDENT_SLEEP_REUSED`.

The 23 statistical positives represent eight external outcomes and seven sleep phenotypes, with 18 exact and five explicitly comparable phenotype matches. They are not 23 independent diseases or 23 independent two-trait replications. Preserve historical `REPLICATED` labels in originals but use qualified external outcome-side validation in the unified ledger. Cohort identity and cross-trait intercepts do not certify zero individual overlap. No independent sleep-source admissibility claim was tested by this reviewer.

### S3 — The seven heterogeneity flags are nominal selected-subset diagnostics

Priority: high for claims of effect differences. `discovery_extension/scripts/23_collate_replication.py:174–176` uses `Var(delta)=SE_discovery²+SE_validation²`. Correct general variance includes `−2Cov(rg_discovery,rg_validation)`. Shared sleep inputs do not establish covariance zero. A cross-trait intercept is a different estimand and cannot be substituted for this covariance.

`replication_results.tsv:2` gives insomnia–abdominal pain: discovery 0.5591 (0.0359), validation 0.4098 (0.0332), delta −0.1493, nominal covariance-zero P=0.00226361. The 41-record independent receipt reproduces all values. Seven of 23 selected positives have nominal P<0.05; across all 41 estimated pairs there are 16 nominal, 10 exploratory BH, and four exploratory Bonferroni positives. These latter counts are diagnostics, not new frozen confirmatory findings.

Selection of discovery associations and subsequent selection of validation positives also affect the distribution of effect differences. Attenuation can be consistent with winner's curse, measurement/ascertainment differences or unknown covariance; no particular cause is identified. Absence of nominal heterogeneity does not establish equivalence.

The archival 45-trait covariance matrix does not contain the complete extension and validation outcome sets and cannot furnish their missing shared-estimator covariance. A new joint jackknife must align actual genomic deletion boundaries across estimates; default block vectors based on different SNP-overlap sets cannot be stacked merely because each has 200 entries. Arbitrary rho grids remain sensitivity diagnostics and are not estimated corrections.

### S4 — Three lipid traits have standalone/pairwise intercept disagreement

Priority: medium for preserved inference; high for new confirmatory models using lipids. The standalone core gate is applied faithfully, but all 36 sleep–LDL/HDL/triglyceride primary correlations have pairwise outcome intercept >1.2. Twenty of these 36 are significant under the original all-396 BH.

| Trait | Recovered h2-summary line | Standalone intercept | Pairwise intercept range | Original significant pair count |
|---|---:|---:|---:|---:|
| LDL | 39 | 1.141 | 1.222–1.247 | 4 |
| HDL | 35 | 1.196 | 1.408–1.423 | 8 |
| Triglycerides | 29 | 1.165 | 1.205–1.214 | 8 |

Examples are recovered `rg_matrix.tsv:16–18`, insomnia with LDL, HDL and triglycerides. All three standalone SSD h2 logs explicitly invoke the two-step estimator. `docs/full_covariance_results.md` already qualifies HDL, and archival multivariate QC retains LDL and triglycerides with intercepts 1.1939199 and 1.1994159. Thus an automatic historical relabeling from pairwise intercepts is not justified. Extend the estimator/SNP-set sensitivity qualification to all three; resolve differences before admitting any new lipid-conditioned analysis. Passing a single historical gate is insufficient to describe them as uniformly clean across estimators.

### S5 — Nonphysical liability h2 remains a source/scale diagnostic

Priority: medium for rg interpretation, high for liability-h2 or genetic-covariance claims. Recovered `h2_summary.tsv:44` retains MS as PASS with h2≈5.982, SE≈1.196; its original SSD log prints 5.9823 (1.1958). Melanoma at line 46 prints 3.731 (1.0456) and is excluded by h2 Z<4. Values above one cannot be reported as literal fractions of liability variance. MS has no significant sleep pair in the frozen family, so this finding does not change the 153 primary-positive count.

Review source identity, per-SNP N/effective N, effect conventions and sample/population prevalence for both traits. A warning alone does not establish validity for new multivariate interpretation. Conversely, liability conversion does not itself change rg; preserve historical classifications rather than treating h2>1 as a silently added exclusion rule. [Official LDSC scale documentation](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) distinguishes scale conversion of h2/covariance from scale-invariant rg and documents the role of intercepts for overlap.

### S6 — Eight of the 17 historical UNDERPOWERED candidates are intercept failures

Priority: medium for classifications. Nine candidates fail h2 Z; eight fail the intercept gate despite h2 Z>22. For example `replication_results.tsv:32` is insomnia–knee arthrosis, source h2 Z=22.71875 and intercept=1.2236. General arthrosis has h2 Z=22.47727 and intercept=1.22 (`replication_source_h2.tsv:11`). Calling all 17 underpowered misstates eight failures. Retain original labels, add explicit `QC_INELIGIBLE_INTERCEPT` versus `QC_INELIGIBLE_H2_Z` reasons, and keep all 17 in the 217 denominator.

### S7 — Exact recovery does not repair stale report content or change correction denominators

Priority: medium for evidence consistency. The checkpoint-matching recovered `phase1_summary_report.md:18–24` still describes three pending sources and only 27 parsed h2 traits, while lines 45–46 state 136 pairs and 51 FDR positives. This is an earlier snapshot, despite its exact expected SHA-256. Preserve it and mark it stale; the recovered numerical tables and checkpoint provide the final family evidence.

The corrected-count display must distinguish 161 all-row numerical BH positives from 153 primary positives. Eight T2D sensitivity positives are excluded from primary inference; melanoma contributes zero. The existing `fdr_primary_phase1` column uses 372 tests and yields 155 positives, which is not the frozen all-396 primary count. Switching correction columns adds insomnia–longevity (`rg_matrix.tsv:23`, q396=0.0505325, q372=0.04992) and accelerometer sleep duration–Crohn disease (`:339`, q396=0.0503556, q372=0.0497610). Do not switch columns, combine the 396/1,200 families, or shrink replication to 41 estimated tests.

## Remaining human/statistical review questions

1. Resolve original source and scale suitability for MS/melanoma and the three lipid estimator disagreements; neither printed logs nor artifact hashes answer these questions.
2. Review cohort overlap and phenotype compatibility for any proposed two-trait replication sources before examining outcomes.
3. Validate any new shared-estimator covariance using a method with aligned genomic blocks, denominator uncertainty and numerical calibration; retain the zero-covariance results as diagnostics until then.
4. For new SEM, inspect archived indefiniteness and conditioning: `ldsc_covariance_diagnostics.tsv:2–4` records three negative genetic-matrix eigenvalues and sampling-V condition number≈1.457×10⁹. The pinned GenomicSEM implementation warns about >1,000 jackknife blocks. Structural finiteness is not calibration; use a justified smaller admitted trait set and document any smoothing without treating it as evidence recovery.
5. Retain correlated phenotypes, shared UKB cohorts and discovery-selected candidates as interpretation constraints. Counts of significant results cannot establish objective-versus-self-report biological differences, distinct disease components, mediation or causality.

This reviewer found no arithmetic implementation error in the original 396/1,200 BH families or 217-candidate threshold classification. Native replay is handled separately by the parent task; this report does not certify those reruns, recover unrounded jackknife statistics, establish independence, identify novel biology, or approve human manuscript authoring.
