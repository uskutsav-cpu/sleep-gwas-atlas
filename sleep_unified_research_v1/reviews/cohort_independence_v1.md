# Separate adversarial cohort-independence review, v1

Review date: 2026-10-08, America/Chicago. Scope: source and cohort structure of the locked 217-candidate historical validation family, its 41 estimated and 23 positive subsets, the result-free MVP clinical-insomnia × ten FinnGen R13 proposal, and the original FinnGen R9 sleep-apnea boundary. This is technical review documentation, not manuscript prose. No new pair outcome was inspected, no GWAS body downloaded, no individual-level data accessed, and no original or SSD input changed.

**Verdict: zero completed historical results qualify as independent both-trait replication. All 41 estimated correlations reuse the discovery sleep GWAS; 23 are threshold-positive qualified external outcome-side validation. The proposed MVP × FinnGen design removes the known shared sleep source relative to original discovery, but remains conditionally admitted clinical-phenotype transport with unknown individual cross-enrollment and incomplete execution gates.**

The review derives counts from the original source queue, results and sleep metadata, independently checks locked row order, recomputes direction and the 0.05/217 success rule, and checks all 41 outcome/sleep paths against eight retained LDSC commands. Existing review conclusions were not used to generate these counts. The machine-readable 217-row evidence matrix is `cohort_independence_v1.tsv`; the checks, input hashes, counts and primary-source access limitations are in `cohort_independence_v1.json`.

## Independently counted original subsets

| Subset | Rows | UKB sleep × FinnGen R13 outcome | UKB sleep × MVP outcome | FinnGen R9 apnea × MVP outcome | No external source |
|---|---:|---:|---:|---:|---:|
| Entire locked candidate family | 217 | 48 | 8 | 2 | 159 |
| Outcome source available | 58 | 48 | 8 | 2 | 0 |
| Correlation estimated | 41 | 36 | 4 | 1 | 0 |
| Threshold-positive | 23 | 19 | 3 | 1 | 0 |

Across all 217 rows, 184 use UKB sleep GWAS and 33 use FinnGen R9 sleep apnea. The 159 source-unavailable rows comprise 128 UKB-sleep and 31 apnea rows. The original statistical categories independently reconcile to 23 threshold-positive, 18 directionally concordant below threshold, 17 QC-ineligible and 159 unavailable. All 41 estimates are directionally concordant. The 17 historically labeled `UNDERPOWERED` comprise nine h2-Z failures and eight intercept failures; all have absent pair estimates. They are not null genetic-correlation results and the intercept failures should not be described purely as low power.

There are 13 selected external source files among the 58 source-available rows. The 23 positives cover eight external source phenotypes and seven sleep traits: 19 FinnGen rows across seven endpoints and four MVP abdominal-pain rows. Phenotype-match metadata are 37 exact /21 comparable among 58 available, 28 exact /13 comparable among 41 estimated, and 18 exact /five comparable among 23 positives. “Exact” is the historical coding/name curation category; it does not establish identical ascertainment, population composition or genetic estimands.

## All four discovery–validation intersections

Let Ds and Do denote discovery sleep and outcome sources, and Vs and Vo validation sleep and outcome sources. The four cross-stage intersections must all be considered; comparing only Do with Vo misses direct sleep reuse and UKB cross-trait overlap.

| Cross-stage intersection | All 217 declared/planned rows | 41 executed estimates | 23 positive estimates |
|---|---|---|---|
| Ds ∩ Vs | Same planned sleep source; 41 actually executed, 176 pair estimates absent | Same sleep input in all 41 | Same sleep input in all 23 |
| Ds ∩ Vo | Cohort-source distinct in 58; no Vo in 159 | Distinct named cohort in all 41 | Distinct named cohort in all 23 |
| Do ∩ Vs | Same UKB cohort in 184, FinnGen R9 versus UKB in 33; Vs is only planned where no pair ran | Same UKB cohort in 40; distinct source for one apnea row | Same UKB cohort in 22; distinct source for one apnea row |
| Do ∩ Vo | UKB versus FinnGen/MVP in 58; no Vo in 159 | Distinct named cohort in all 41 | Distinct named cohort in all 23 |

For Ds ∩ Vs, shared source identity is established by source IDs, discovery script `14_rg_extension.sh:36–40`, validation script `49_rg_replication.sh:16–29`, and each retained command in `logs/replication/rg/`. This is historical code/run-path evidence; this reviewer did not compare missing historical source bodies byte for byte. Known reuse is enough to reject independent both-trait replication. The one sleep-apnea positive is FinnGen R9 sleep apnea × MVP abdominal pain, with rg 0.3433 in discovery and 0.3920 in validation; it still reuses sleep apnea. It is not a FinnGen R9 × R13 positive.

Within discovery, UKB sleep and Pan-UKB outcomes share the named cohort for 184 candidates (40 estimated /22 positive); actual pairwise overlapping participant counts are unknown. FinnGen R9 sleep apnea versus Pan-UKB is source-distinct for 33 candidates (one estimated/positive). Within historical validation, all 58 source-available pair designs use different named cohorts: UKB–FinnGen, UKB–MVP or FinnGen R9–MVP. Actual cross-enrollment is not measured. No validation correlation exists for the other 176 rows, so planned cohort relations must not be represented as executed tests.

## Cohort-source distinction is not an individual-ID audit

The original source queue labels all 58 available sources `NON_OVERLAPPING_CONFIRMED`, but its evidence text describes only the independent organization/cohort relative to Pan-UKB. It does not cover every cross-stage intersection and supplies no individual enrollment intersection or numerical upper bound. Therefore the supported replacement description in a new ledger is **cohort-source distinct; individual intersection unknown**. Preserve the original label as historical data, with its qualification visible alongside it. A public cohort name, different countries, a later release, European ancestry, or a nonsignificant cross-trait intercept cannot certify zero shared people.

The [official FinnGen R9 documentation](https://finngen.gitbook.io/documentation/r9) identifies Finnish biobanks and health registries; the [primary MVP resource](https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/) identifies MVP participants and VA EHR/enrollment phenotypes. These support distinct named-cohort designs, not proof against individual cross-enrollment. The original source queue uses one UKB-only Jansen insomnia source, not the full Jansen UKB/23andMe meta-analysis; mixing those versions would change the overlap assessment.

LDSC may model overlap-associated covariance under its assumptions. The [authors’ documentation](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) and [primary method](https://www.nature.com/articles/ng.3406) support estimating the cross-trait intercept. That does not establish independent source provenance, identify shared people, or estimate covariance between two different rg estimates. Source reuse does not by itself invalidate each historical global rg estimate; it invalidates the strongest independence claim and changes the uncertainty needed for comparing estimates.

## FinnGen R9 sleep-apnea boundary

The exact original sleep-apnea source is `finngen_r9_sleep_apnoea`, release `FinnGen_R9_G6_SLEEPAPNO_2023-05-04`, 38,998 cases/336,659 controls, N=375,657. The queue admits only two MVP outcomes for apnea: abdominal pain (`GCST90479148`, estimated and threshold-positive) and drug-allergy history (`GCST90479330`, h2-ineligible and unestimated). Ten otherwise selected FinnGen R13 outcomes were explicitly rejected for sleep apnea because the existing sleep source already uses FinnGen R9. Those ten rows and their exact rejection evidence are flagged in the matrix. The other 21 apnea-unavailable rows lack another admitted source.

This exclusion is appropriate. A new FinnGen release is not a cohort-disjoint validation cohort. [FinnGen’s release FAQ](https://docs.finngen.fi/faq/about-finngen-data/why-are-there-differences-in-the-gwas-results-between-data-freezes-releases) documents changing study subjects, legacy composition, endpoints and pipelines; it does not supply a numerical R9/R13 membership intersection. Do not claim a measured overlap fraction or assume strict nesting. Treat the shared FinnGen source as an independence failure until an explicitly disjoint sample and compatible phenotype are documented. Smaller P or larger N in R13 cannot repair that design.

The core panel also contains FinnGen R9 sources for MS, asthma, T2D, CAD and melanoma. Core sleep-apnea correlations with those endpoints share FinnGen source cohorts. They can remain overlap-aware discovery analyses, but cannot function as independent confirmations merely because their phenotype names differ. This observation is source-metadata-only; no new core pair outcomes were examined here.

## Result-free MVP insomnia × ten FinnGen outcomes

The proposed source contract selects EUR `GCST90475826`, MVP insomnia PheCode 327.4, 78,566 cases/329,572 controls, N=408,138. The [primary dbGaP GIA phenotype manifest](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?phd=8759&study_id=phs002453.v1.p1) independently supports those exact counts. Its separate HARE row, all-population row and organic/persistent-insomnia endpoint have different counts and must not be substituted. The accession itself could not be live recertified through the web open tool in this review; exact file identity therefore relies on the pinned source contract and remains subject to full-body verification.

The original discovery pair design is UKB questionnaire insomnia × Pan-UKB outcome. Proposed validation is MVP clinical insomnia × FinnGen R13 outcome. All four cross-stage source comparisons are distinct by named cohort: Ds–Vs=UKB–MVP; Ds–Vo=UKB–FinnGen; Do–Vs=UKB–MVP; Do–Vo=UKB–FinnGen. Within validation, MVP–FinnGen is also cohort-source distinct. Every individual-level intersection remains unknown. Within original discovery, both traits still share UKB.

The proposed ten rows use six historical exact outcome matches and four comparable matches. They span six historical positives, one directionally concordant estimate and three historically QC-ineligible outcomes. Including all ten source-selected insomnia rows rather than only the six positives preserves the declared result-free scope. The two MVP historical outcomes are properly excluded from this pilot: pairing new MVP insomnia with an MVP outcome would create within-validation shared-cohort sampling and would need a separate protocol. It is not scientifically interchangeable with the proposed cross-cohort design.

The FinnGen R13 outcomes are the same outcomes used in the original outcome-side validation. Thus a successful future pilot would be independently sourced on both traits relative to **original discovery**, while sharing the outcome data with **earlier validation**. It cannot be counted as two mutually independent validation waves or pooled as if the old and new estimates were statistically unrelated.

MVP clinically coded insomnia and UKB frequent self-reported sleep complaints are related, different phenotypes. The frozen protocol correctly limits success to **cohort-distinct clinical-insomnia transport validation**. It does not authorize exact questionnaire-insomnia replication, causal interpretation, a new primary biological question or local/molecular follow-up. Unknown cross-enrollment and phenotype transport must remain in any result label. The stricter original contract asks for no participant overlap or quantitatively demonstrated negligible overlap; its individual-level condition has not been certified here. Source distinction alone must not be described as having satisfied that stronger condition.

At the source-contract snapshot, full compressed-body hash verification, source-specific schema/effect/SE/N/build checks, current source harmonization, sleep/outcome h2 and intercept gates, and prospective power checks are incomplete. This review did not inspect any newly estimated MVP pair result. This is a design-admissibility verdict, not a positive replication verdict. Keep all 217 rows and 0.05/217=0.0002304147465437788. Preserve QC failures and nulls, and do not remove the three formerly QC-ineligible outcomes or replace a source after seeing a new pair result.

## Shared-sleep effect differences and winner’s curse

The original difference variance is `SE_discovery² + SE_validation²`. The variance actually needed is `SE_discovery² + SE_validation² − 2 Cov(rg_discovery,rg_validation)`. Sleep-source reuse supplies common noisy GWAS information and common estimated sleep-h2 denominators; UKB discovery outcome–validation sleep cohort overlap adds another dependence route in 40 of 41 pairs. The sign and magnitude of total covariance have not been established. A within-pair cross-trait intercept is not that between-estimate covariance, and an arbitrary correlation grid is not a covariance estimator.

The independently counted zero-covariance nominal P<0.05 flags are **16 among all 41 estimates**, of which **seven among the 23 selected positives**. Referring to “seven heterogeneity findings” without specifying the positive subset obscures the other nine. All flags remain uncalibrated selected-subset diagnostics until a valid covariance estimator and selection-aware interpretation are supplied. The absence of a flag does not establish effect equivalence.

Absolute validation rg is smaller in 36/41 estimates and 19/23 positives. Those are descriptive counts, not evidence that a particular mechanism, population difference or winner’s curse explains attenuation. The family was selected after discovery FDR, magnitude, QC and novelty screening; the 23-positive subset is additionally selected on validation results. Winner’s curse is therefore a plausible threat to discovery magnitudes and naive difference tests. Reusing sleep data means it is not automatically removed by swapping the outcome. Clinical ascertainment, phenotype differences, uncertainty, power and sampling covariance offer additional explanations. No shrinkage correction, selection correction or calibrated causal attribution was estimated here.

Joint genomic deletion boundaries, compatible SNP intersections, native ratio-rg jackknife marginal SE agreement and independent implementation agreement are needed before covariance-aware effect comparisons can pass. Differently partitioned 200-delete vectors cannot merely be paired by their row index. A corrected covariance alone also would not undo discovery selection or make questionnaire and clinical insomnia identical. h2 Z≥4 is an eligibility rule, not demonstrated power to reject a small rg; threshold-negative results require the frozen detectable-effect/power report.

## Exact human scientific review requirements

1. **Cohort/genetic epidemiology reviewer:** sign the Ds/Do/Vs/Vo membership matrix for the exact UKB-only sleep sources, Pan-UKB EUR, two historical MVP accessions, eleven selected FinnGen R13 endpoints and the proposed MVP insomnia accession. Separate absence of a shared named cohort from absence of individually cross-enrolled or related participants. If a strong no/negligible-overlap claim is needed, use authorized membership/linkage evidence or a justified quantitative upper bound and sensitivity assessment; otherwise retain unknown individual overlap and the cohort-distinct label. No investigator contact or restricted-data access is authorized by this review.
2. **FinnGen release expert:** verify the R9/R13 shared-sample risk using release/sample-inclusion documentation; retain the ten sleep-apnea exclusions absent an explicitly disjoint design. Review the five FinnGen R9 core-outcome relations as within-cohort discovery analyses, not independent confirmations. Do not invent a R9/R13 overlap percentage from release totals.
3. **Sleep clinician/phenotyper:** approve the clinical-to-questionnaire transport estimand, exact PheCode/ICD mapping, case/control exclusions, and each of the four outcome-definition differences in the ten-row pilot. Approve only related-phenotype transport language; frequent complaints are not automatically insomnia disorder. Review MVP veteran demographics and Finnish registry ascertainment as generalizability constraints.
4. **Statistical geneticist:** decide whether covariance estimation is sufficient for discovery–validation differences; review aligned deletion blocks, ratio-h2 denominator uncertainty, original selection, correlated sleep/outcome traits and reusing FinnGen outcomes across validation waves. Treat all 16 nominal cov0 flags and the seven-positive subset as diagnostics until the relevant calibration is established. Determine whether any effect-consistency claim between clinical and questionnaire insomnia is a scientifically valid estimand at all.
5. **Source/QC reviewer:** verify whole-file identity and current harmonization before testing; require identical inherited source gates for the ten outcomes, including three previously QC-ineligible rows, and prospective power/detectable-effect reporting. Confirm the 217-family denominator, result-independent stop and complete negative/unavailable ledger. A successful software run or passing h2 threshold cannot substitute for these checks.

The historical evidence supports qualified outcome-side validation and source-aware discovery inference. It does not yet support a claim that any completed result is an independent replication of both traits. A future MVP pilot can improve source independence relative to discovery only after the frozen gates succeed and the human review accepts the qualified phenotype/cohort interpretation.

## Ten-row pilot source list (no new pair outcomes)

| Original pair ID | Existing FinnGen outcome | Historical class | Outcome match |
|---|---|---|---|
| `insomnia__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | REPLICATED | EXACT |
| `insomnia__panukbb_phecode__550_2__both_sexes__na__na` | `finngen_r13_K11_DIAHER` | REPLICATED | EXACT |
| `insomnia__panukbb_phecode__558__both_sexes__na__na` | `finngen_r13_K11_OTHENTERCOL` | REPLICATED | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| `insomnia__panukbb_phecode__496_2__both_sexes__na__na` | `finngen_r13_J10_BRONCHNAS` | UNDERPOWERED | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| `insomnia__panukbb_phecode__318__both_sexes__na__na` | `finngen_r13_SMOKING_DEPEND` | REPLICATED | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| `insomnia__panukbb_icd10__j44__both_sexes__na__na` | `finngen_r13_J10_COPD` | REPLICATED | EXACT |
| `insomnia__panukbb_icd10__m17__both_sexes__na__na` | `finngen_r13_M13_ARTHROSIS_KNEE` | UNDERPOWERED | EXACT |
| `insomnia__panukbb_phecode__740__both_sexes__na__na` | `finngen_r13_M13_ARTHROSIS` | UNDERPOWERED | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| `insomnia__panukbb_icd10__m51__both_sexes__na__na` | `finngen_r13_M13_INTERVERTEB` | REPLICATED | EXACT |
| `insomnia__panukbb_icd10__k80__both_sexes__na__na` | `finngen_r13_K11_CHOLELITH` | DIRECTIONALLY_CONCORDANT | EXACT |

## Twenty-three historical positives with row-level evidence

All rows below reuse sleep and have unknown individual cross-enrollment between distinct named cohorts. These are qualified outcome-side validations. The TSV retains each source-queue and historical command line, all four cross-stage relations, and the entire absent/QC-failed/nonsignificant family.

| Locked row | Pair ID | External source | Outcome match |
|---|---|---|---|
| 1 | `insomnia__panukbb_phecode__785__both_sexes__na__na` | `gwas_catalog_GCST90479148` | EXACT |
| 4 | `insomnia__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | EXACT |
| 10 | `insomnia__panukbb_phecode__550_2__both_sexes__na__na` | `finngen_r13_K11_DIAHER` | EXACT |
| 12 | `shortsleep__panukbb_phecode__785__both_sexes__na__na` | `gwas_catalog_GCST90479148` | EXACT |
| 15 | `insomnia__panukbb_phecode__558__both_sexes__na__na` | `finngen_r13_K11_OTHENTERCOL` | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| 23 | `insomnia__panukbb_phecode__318__both_sexes__na__na` | `finngen_r13_SMOKING_DEPEND` | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| 27 | `insomnia__panukbb_icd10__j44__both_sexes__na__na` | `finngen_r13_J10_COPD` | EXACT |
| 30 | `shortsleep__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | EXACT |
| 34 | `shortsleep__panukbb_phecode__318__both_sexes__na__na` | `finngen_r13_SMOKING_DEPEND` | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| 40 | `shortsleep__panukbb_phecode__558__both_sexes__na__na` | `finngen_r13_K11_OTHENTERCOL` | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| 42 | `sleep_apnea__panukbb_phecode__785__both_sexes__na__na` | `gwas_catalog_GCST90479148` | EXACT |
| 51 | `shortsleep__panukbb_icd10__j44__both_sexes__na__na` | `finngen_r13_J10_COPD` | EXACT |
| 55 | `shortsleep__panukbb_icd10__m51__both_sexes__na__na` | `finngen_r13_M13_INTERVERTEB` | EXACT |
| 60 | `longsleep__panukbb_phecode__318__both_sexes__na__na` | `finngen_r13_SMOKING_DEPEND` | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| 68 | `insomnia__panukbb_icd10__m51__both_sexes__na__na` | `finngen_r13_M13_INTERVERTEB` | EXACT |
| 74 | `longsleep__panukbb_phecode__785__both_sexes__na__na` | `gwas_catalog_GCST90479148` | EXACT |
| 119 | `snoring__panukbb_icd10__k80__both_sexes__na__na` | `finngen_r13_K11_CHOLELITH` | EXACT |
| 130 | `snoring__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | EXACT |
| 141 | `longsleep__panukbb_icd10__j44__both_sexes__na__na` | `finngen_r13_J10_COPD` | EXACT |
| 172 | `napping__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | EXACT |
| 178 | `sleepdur__panukbb_icd10__m51__both_sexes__na__na` | `finngen_r13_M13_INTERVERTEB` | EXACT |
| 184 | `longsleep__panukbb_icd10__k21__both_sexes__na__na` | `finngen_r13_K11_REFLUX` | EXACT |
| 204 | `longsleep__panukbb_icd10__m51__both_sexes__na__na` | `finngen_r13_M13_INTERVERTEB` | EXACT |
