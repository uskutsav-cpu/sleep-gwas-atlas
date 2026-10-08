# Adapted STROBE/STREGA reporting evidence map

Execution/access date: **2026-10-08**. This is an evidence-completion checklist for a human author, not a completed journal manuscript checklist. **All manuscript page/line pointers are PENDING_AUTHOR.** No manuscript exists in this task and no reporting-compliance certification is issued.

Primary references: [STROBE official checklists](https://www.strobe-statement.org/checklists/), [combined cohort/case-control/cross-sectional checklist PDF](https://www.strobe-statement.org/download/strobe-checklist-cohort-case-control-and-cross-sectional-studies-combined), [STROBE explanation article DOI10.1371/journal.pmed.0040296](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0040296), [STREGA statement DOI10.1371/journal.pmed.1000022](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000022), and [primary STREGA article PDF, Table1](https://stacks.cdc.gov/view/cdc/50845/cdc_50845_DS1.pdf). Original STROBE checklist and STREGA full article were successfully retained and inspected locally; a separate PMC challenge is a failed fetch. Hash/access ledgers retain this distinction.

STROBE was developed for observational cohort, case-control and cross-sectional reports; STREGA extends genetic-association reporting. This work analyzes aggregate GWAS summary statistics and does not recruit participants or perform local genotyping. STROBE/STREGA therefore require an explicit adaptation, not an unsupported claim of conventional participant-level compliance. Items concerning source participant selection, genotype QC, medication treatment, population stratification and relatedness remain applicable to upstream provenance even when no local test was performed. Missing upstream information is **UNRESOLVED**, never automatically not applicable. MR-specific checklists are inapplicable because no MR analysis was executed.

`APPLICABLE`/`ADAPTED` means an evidence requirement, not a pass. Every row is **PARTIAL_OR_AUTHOR_PENDING** unless explicitly local-not-applicable; links to package paths are evidence pointers that must be checked against the final authored report. These are paraphrased audit criteria; use the original official checklists for journal-upload forms.

| Item | Audit criterion | Applicability | Evidence pointer | Required qualification/completion |
|---|---|---|---|---|
| 1a | Study design in title/abstract | APPLICABLE | `ANALYSIS_PROTOCOL.md` | Author identify summary-GWAS association design |
| 1b | Balanced abstract | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md;06_NOVELTY_AND_PRIOR_ART.md` | Human abstract pending |
| 2 | Scientific rationale | APPLICABLE | `PRIOR_STUDY_COMPARISON.md` | Author rationale without inflated novelty |
| 3 | Objectives/hypotheses | APPLICABLE | `ANALYSIS_PROTOCOL.md` | Separate historical exploratory family from audits |
| 4 | Design early in report | APPLICABLE | `ANALYSIS_PROTOCOL.md` | Human Methods pending |
| 5 | Sources/settings/dates | APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Resolve original dates/release/source subsets |
| 6a | Eligibility and sampling | ADAPTED | `02_SOURCE_AND_PROVENANCE_AUDIT.md;tables/source_integrity_ledger.tsv` | Summarize upstream selection and audit trait selection |
| 6b | Matching/follow-up cohorts | NOT_APPLICABLE_LOCAL | `Aggregate GWAS inputs only` | No local matched recruitment or participant follow-up |
| 7 | Variables/outcomes/confounders | APPLICABLE | `tables/source_integrity_ledger.tsv;tables/replicated_23.tsv` | Source coding and covariates incomplete |
| 8 | Measurement/data sources | APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md;04_REPLICATION_INDEPENDENCE_AUDIT.md` | Exact/comparable phenotype distinctions required |
| 9 | Bias safeguards | APPLICABLE | `04_REPLICATION_INDEPENDENCE_AUDIT.md;05_HETEROGENEITY_VALIDATION.md` | UKB overlap,selection,measurement and shared sleep GWAS |
| 10 | Study-size rationale | ADAPTED | `ANALYSIS_PROTOCOL.md;tables/global_1200.tsv` | 12×100 tests are not participant N |
| 11 | Quantitative transformations/cutoffs | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md;tables/source_integrity_ledger.tsv` | Retain source scales,thresholds and missingness |
| 12a | Statistical methods | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md;05_HETEROGENEITY_VALIDATION.md` | Native LDSC replay blocked; distinguish numerical audit |
| 12b | Subgroups/interactions | ADAPTED | `03_GLOBAL_RESULTS_VERIFICATION.md` | Exploratory strata only |
| 12c | Missing-data handling | APPLICABLE | `tables/replication_family_217.tsv` | Keep unavailable/QC-ineligible outcomes and missing h2 |
| 12d | Follow-up/matching/sampling weights | NOT_APPLICABLE_LOCAL | `No individual-level cohort or survey reanalysis` | Record source study treatment if relevant |
| 12e | Sensitivity analyses | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md;05_HETEROGENEITY_VALIDATION.md` | Precision/covariance grids are audits,not native GWAS robustness |
| 13a | Numbers at each stage | ADAPTED | `figures/source_data/fig3_flow.tsv;tables/replication_family_217.tsv` | Tests,phenotypes and participants have different denominators |
| 13b | Reasons for exclusion | APPLICABLE | `tables/replication_family_217.tsv;02_SOURCE_AND_PROVENANCE_AUDIT.md` | Retain failed/unavailable stages |
| 13c | Flow diagram | ADAPTED | `figures/fig3_flow.pdf` | Trait/test flow does not substitute for upstream cohort flow |
| 14a | Participant characteristics | ADAPTED | `tables/source_integrity_ledger.tsv` | Upstream metadata incomplete; no reconstructed participant demographics |
| 14b | Missingness by variable | ADAPTED | `tables/source_integrity_ledger.tsv;tables/replication_family_217.tsv` | Local metadata missingness plus source individual-level missingness unresolved |
| 14c | Participant follow-up time | NOT_APPLICABLE_LOCAL | `Aggregate GWAS association reanalysis` | No participant follow-up performed |
| 15 | Outcome numbers/summary | ADAPTED | `tables/source_integrity_ledger.tsv` | Source case/control,N and prevalence require verified upstream records |
| 16a | Estimates and uncertainty | APPLICABLE | `tables/global_1200.tsv;tables/replicated_23.tsv` | rg/SE/P/FDR/CI; do not invent adjusted participant models |
| 16b | Category boundaries | APPLICABLE | `tables/source_integrity_ledger.tsv` | Different sleep and diagnosis thresholds explicit |
| 16c | Absolute-risk translation | NOT_APPLICABLE | `Global rg does not estimate risk or treatment effect` | Do not invent clinical risk benefit |
| 17 | Additional analyses | APPLICABLE | `05_HETEROGENEITY_VALIDATION.md;tables/HETEROGENEITY_MASTER.tsv` | Qualify nominal heterogeneity and unknown covariance |
| 18 | Principal results | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md` | Author Results pending |
| 19 | Limitations/direction of bias | APPLICABLE | `04_REPLICATION_INDEPENDENCE_AUDIT.md;06_NOVELTY_AND_PRIOR_ART.md` | No human discussion supplied |
| 20 | Cautious interpretation | APPLICABLE | `SLEEP_EDITORIAL_NOVELTY_ASSESSMENT.md` | No causal,molecular,clinical or first-priority conclusion |
| 21 | Generalisability | APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | EUR-focused sources and selected outcomes limit transfer |
| 22 | Funding and funder roles | APPLICABLE | `15_ETHICS_AUTHOR_APPROVAL_CHECKLIST.md` | Investigator input; no invented None |
| STREGA3 | First report or replication | APPLICABLE | `06_NOVELTY_AND_PRIOR_ART.md` | Outcome-side validation and substantial prior overlap |
| STREGA6a | Genetic subset selection | ADAPTED | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Upstream subset selection remains partially unresolved |
| STREGA7b | Genetic nomenclature/stratification variables | ADAPTED | `tables/source_integrity_ledger.tsv` | GWAS/build/ancestry/allele/source identity |
| STREGA8b | DNA/genotyping/platform/call/error rates | NOT_PERFORMED_LOCAL;UPSTREAM_APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Cite upstream genotype QC; unknown rates remain missing |
| STREGA9b | Treatment-related quantitative bias | UPSTREAM_APPLICABLE | `tables/source_integrity_ledger.tsv` | Source medications/covariates unresolved where relevant |
| STREGA11 | Treatment effect transformations | UPSTREAM_APPLICABLE | `tables/source_integrity_ledger.tsv` | No local participant medication adjustment |
| STREGA12a | Software/version/settings | APPLICABLE | `REPRODUCE.md;requirements/;scripts/` | Pinned numerical environment differs from unrecovered native runtime |
| STREGA12f | Hardy-Weinberg handling | NOT_PERFORMED_LOCAL;UPSTREAM_APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | No genotype-level test; source QC descriptions required |
| STREGA12g | Genotype/haplotype inference | NOT_PERFORMED_LOCAL;UPSTREAM_APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | No local inference; source imputation/platform provenance required |
| STREGA12h | Population stratification | APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Source PCs/ancestry,EUR LD references,intercepts/limitations |
| STREGA12i | Multiple testing | APPLICABLE | `03_GLOBAL_RESULTS_VERIFICATION.md` | Preserve1200 BH,217 Bonferroni,heterogeneity multiplicity |
| STREGA12j | Relatedness handling | UPSTREAM_APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Source exclusion/mixed models unresolved; participant overlap audit |
| STREGA13a | Genotyping attempted/successful | NOT_PERFORMED_LOCAL;UPSTREAM_APPLICABLE | `02_SOURCE_AND_PROVENANCE_AUDIT.md` | Source counts required; no invented local flow |
| STREGA14a | Genetic descriptive information | ADAPTED | `tables/source_integrity_ledger.tsv` | N/case-control/h2/nSNP and metadata gaps |
| STREGA15 | Outcome data by genotype | NOT_APPLICABLE_LOCAL | `No individual genotype/outcome records` | Not estimable from global rg |
| STREGA16d | Multiplicity-adjusted estimates | APPLICABLE | `tables/global_1200.tsv;tables/replication_family_217.tsv` | Adjusted probabilities; rg not shrunk by BH |
| STREGA17b | Complete genetic results | APPLICABLE | `tables/global_1200.tsv;tables/replication_family_217.tsv` | Include all1200 and all217,not positives only |
| STREGA17c | Detailed data access | APPLICABLE | `13_DATA_AND_CODE_AVAILABILITY_INVENTORY.md;tables/source_integrity_ledger.tsv` | Access/licensing and missing raw inputs qualified |

Upstream-genotyping items also require appropriate source citations, rather than claiming the package reproduced genotyping. SNP-level allele frequencies, candidate-variant counts, haplotype distributions and genotype-stratified risk are not fabricated. Source restrictions and ethics decisions require investigators; checklist completion does not create data-use permission.
