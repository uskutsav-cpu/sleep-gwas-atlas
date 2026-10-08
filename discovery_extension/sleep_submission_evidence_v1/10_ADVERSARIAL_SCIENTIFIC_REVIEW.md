# Automated adversarial scientific review

Four separate automated perspectives were used: root integration/claim review, numerical review with an independent incomplete-gamma tail implementation, provenance/source review, and primary literature/full-supplement review. Cross-review identified and corrected a clipped forest interval, incomplete source-to-claim mapping, ambiguous phenotype labels, a real-data dependency in the synthetic test, and inadequate exact-universe assertions. These are automated checks, not human specialist approvals.

| Challenge | Evidence and finding | Disposition |
|---|---|---|
| Are 23 validations independent for both traits | Both original rg scripts reuse data/munged sleep inputs | Fully independent two-trait claim rejected; outcome-side classification retained |
| Are definitions equivalent | 18 exact and five documented comparable definitions; full 217-row source comparison | Comparable definitions remain marked; phenotype sensitivity is not exact replication |
| Are correlations causal | No causal design or verified local/molecular evidence was executed | Causal, gene, mechanism and clinical-benefit claims prohibited |
| Does UKB reuse affect discovery | Most sleep inputs and Pan-UKB outcomes share a cohort; FinnGen sleep apnea has a different cohort | Intercepts retained; no intercept-based proof of independence |
| Were correction families changed | Exact 1,200 BH and 217 Bonferroni universes and original values checked | PASS arithmetic; no favorable exclusions |
| Was the 100-phenotype panel prospective | Historical lock/receipts and Git temporal evidence reviewed | Evidence supports before-rg selection; retrospective evidence cannot prove undocumented analyst behavior |
| Was replication reporting selective | All 217 include 159 unavailable, 17 QC/power-ineligible, 18 directional and 23 positives | Complete selected family retained; selected family does not represent all 603 hits |
| Do seven heterogeneous pairs weaken consistency | Seven nominal flags survive zero-covariance recalculation; covariance grid changes inferential counts | Qualified; seven is not a family-corrected discovery count |
| Do correlated diseases inflate breadth | Eight external phenotypes; 28 profile similarities 0.775–0.985 | Comparisons are not independent diseases/loci; similarity is not phenotype rg |
| Has SLEEP or other literature already reported these | Morrison full supplement plus 2026 full tables overturn several historical no-direct-rg labels | Novelty weakened; use current pair table, not historical labels |
| Are original inputs missing | Nine checkpoint originals missing; original sleep h2 incomplete; no local GWAS/LD replay inputs | Full source reproduction BLOCKED |
| Do headline positives use weak measurement | 20 questionnaire self-report positives, two partner-report snoring positives and one registry phenotype; zero actigraphy validations | Clinical depth and generalizability remain substantial limitations |
| Does every claim map to data | 1,200-row claim ledger; full 217 sources; exact figure source hashes | Machine consistency gate required |
| Are positives sensitive to unsupported assumptions | Rounding cannot change the primary counts; heterogeneity depends on unknown covariance | Primary arithmetic robust, inferential consistency qualified |
| Were unavailable local analyses called null | Local, PLACO, fine-mapping and coloc results remain missing upstream | No zero/null substitution |
| Were new analyses confirmatory | Sensitivities explicitly retrospective; no new confirmatory family | No confirmatory promotion |
| Can a specialist reproduce central conclusions | Numerical/figure replay possible; current full literature cache acquisition and original GWAS replay have separate boundaries | CONDITIONAL handoff; native science reproduction not passed |

The central first-ever novelty and fully independent two-trait replication interpretations were falsified. The narrower findings remain an exploratory global association screen with outcome-side validation and explicitly limited measurement depth. No new causal or molecular finding emerged. Specialist human statistical-genetics and sleep-science review remains necessary.

Independent review receipts are in `logs/provenance_independent_review.json`, `logs/numerical_independent.json`, `logs/numerical_figure_source_audit.json`, and `logs/literature_audit_receipt.json`; exact receipt names are listed in the final hash inventory. Historical test failures are retained separately and cannot be converted to scientific passes.
