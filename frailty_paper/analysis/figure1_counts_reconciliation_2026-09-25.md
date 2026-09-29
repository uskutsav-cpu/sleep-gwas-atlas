# Figure 1 review-count reconciliation and visual QA

Review date: 2026-09-25.

## Discrepancy and correction

The saved Figure 1 provenance and source-data table described an older review snapshot: 60,989 PubMed records, 4,897 duplicate occurrences, and a 56,092-record queue. The current `review/prisma_flow.tsv` and `review/source_record_counts.tsv` agree on 61,009, 4,892, and 56,117, respectively. The difference was in the figure generator's hard-coded counts and rendered artifacts; no screening or inclusion status had changed.

Updated `scripts/37_plot_evidence_framework.py` to derive the counts from both ledgers and fail if the ledgers disagree, the source-minus-duplicates arithmetic does not reconcile, manual-database records appear without a revised PubMed-only panel, or screening/inclusion status changes. The screening note now says the current deduplicated PubMed queue is ready while screening and eligibility counts remain pending.

Regenerated the source-data TSV, SVG, vector PDF, 4,000 × 3,000 PNG, and provenance JSON. The new flow counts are 61,009 identified, zero manual-database imports, 4,892 duplicate occurrences linked, 56,117 queued, zero screened, zero reports assessed, and zero included. All recorded output SHA-256 values match the provenance JSON. Visual inspection of the PNG confirms the corrected figures are legible and the four-panel layout remains unclipped.

The separate 17-occurrence PubMed increment dated 2026-09-23 through 2026-09-25 remains outside this figure: the frozen review cutoff is 2026-09-22, and the increment remains in its separate cache. Query hits are not screening decisions.

## Validation

- Focused Figure 1 count-reconciliation tests: 7 passed.
- Existing Figure 1 analysis summary preserved: 12 FI pairs (9 significant under inherited all-396 BH), 84 latent-factor sensitivity pairs (48 significant), and 72 aging-context pairs (26 significant).
- No review queue, reviewer decision, analysis result, or scientific threshold was changed.
- Figure 1 remains a provisional PubMed-only acquisition/design framework, not a completed PRISMA flow or final manuscript figure.

Provenance and generated file hashes: `figure1_evidence_framework.provenance.json`.
