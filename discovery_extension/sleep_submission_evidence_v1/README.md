# Sleep GWAS research and submission support

Overall SLEEP assessment: **NO_GO_FOR_SLEEP_AS_CURRENTLY_SUPPORTED**. Human manuscript handoff: **CONDITIONAL**. This is a completed non-manuscript evidence package with explicit external-data, scientific and human blockers; no manuscript or submission was created.

Start with [18 Manuscript author handoff](18_MANUSCRIPT_AUTHOR_HANDOFF.md), [17 Readiness](17_SLEEP_EDITORIAL_READINESS.md), [Final blockers](FINAL_BLOCKERS_AND_RESTART_ACTIONS.md), and [Final execution receipt](FINAL_EXECUTION_RECEIPT.json). The 18 numbered audits cover science, prior art, sources, ethics, reporting and journal rules. The main-source baseline is recorded in the receipt. Historical outputs remain unchanged.

All 1,200 discovery comparisons and all 217 validation candidates were checked numerically at archived precision. The 603 BH discoveries remain. All 23 original positives remain Bonferroni positive but are **external-outcome-side replication**, with the same sleep GWAS; none is fully independent two-trait replication. Seven nominal heterogeneity flags assume zero estimator covariance. Nine original core checkpoint files remain missing. Native GWAS/LDSC reruns: zero.

The `tables/` directory contains canonical TSVs, a nine-sheet editable `supplementary_data.xlsx`, and six editable numerical tables under `tables/word/`. The `figures/` directory contains five research figure types and one preliminary graphical abstract, seven pages total, with PDF/vector/300-dpi PNG, source data and alt-text metadata. These are separate support assets, not a manuscript. Final figure selection, print scaling and journal portal supplement rules need human review.

Independent receipts: `logs/numerical_independent.json`, `logs/numerical_figure_source_audit.json`, `logs/provenance_independent_review.json`, `logs/workbook_independent_validation.json`, `logs/visual_qa.json`, `logs/clean_checkout_receipt.json`. All 25 new tests pass. The historical full test suite has 12 failures and 37 errors, retained explicitly; passing source-free tests cannot certify native reproduction.

Use [REPRODUCE](REPRODUCE.md) and `scripts/finalize_package.py --verify` to check the released hash inventory from the repository root. A standalone output copy is a review packet; reproduction uses the dedicated Git branch because canonical historical inputs remain in their original repository paths. New full publisher/supplement/search caches and raw GWAS are excluded from public release; exact acquisition URLs and hashes remain in `sources/`. No new project-code license, public archive or DOI has been asserted.

Machine gates: `tables/final_acceptance_gates.tsv`. PASS concerns the named evidence level only. FAIL, BLOCKED and NOT_APPLICABLE rows are deliberate scientific boundaries and cannot be promoted by software test success.
