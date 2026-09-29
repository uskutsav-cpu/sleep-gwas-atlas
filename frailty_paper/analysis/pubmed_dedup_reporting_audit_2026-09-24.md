# PubMed deduplication and reporting audit — 2026-09-24

- Audit source HEAD: `a11985d2b120bc17cc5154f177932e57311e3b4b` on `frailty-paper-v1`.
- Review-builder manifest SHA-256: `9585af71dc1cd50e3a798e971721bd0c381e175320f122fbabcaee37f0e04bbb`.
- Output-manifest SHA-256: `7301215e3d9fa40cfd479b0afbc4cc73b3cd3e8c1c4c227bda22e6055480379c` (matches `record_build_manifest.json`).
- Builder SHA-256: `a21916abc8c47ac434a24d3ee5e7b0623082358363b9e026f8abcc5e03cc7f25` (matches the recorded build).
- All 10 outputs in `review/record_build_outputs.tsv` reverified by declared byte size and SHA-256.

The verified PubMed-only build parsed 311 XML batches and 60,989 source records, retained 56,092 records, linked 4,897 duplicate occurrences, and imported zero manual-database records. The builder uses exact PMID, then normalized DOI when title evidence agrees or a unique title candidate exists, then normalized-title matching only when identifiers do not conflict; otherwise records remain separate. Duplicate memberships and source/query IDs are retained in the audit files, and source-priority selection is recorded in the manifest.

The current title/abstract screening validator passes with 56,092 records and zero decisions; full-text records, exclusions and included studies remain zero. No review-completion claim is supported. Embase, Scopus, Web of Science and PsycINFO exports remain unsupplied, so the full multi-database deduplication and screening remain pending.

The manuscript Methods now describes the completed PubMed-only build and explicitly limits those counts to PubMed. PRISMA-S item 16 is `PARTIAL_PUBMED_SNAPSHOT_ONLY` because the manuscript reports only this partial record-management step and licensed-source integration remains outstanding. The focused post-edit checks passed: 28/28 quantitative claims (`review_dedup_claim_audit_2026-09-24_1422.log`, SHA-256 `61dcc61d68246e2094870173ab9ec3123fe19a07c60c09cf571412da743df356`) and 92 checklist rows / 93 line references / zero errors (`review_dedup_locator_audit_2026-09-24_1422.log`, SHA-256 `ff6ffe71b62cd2d678dfbd16e2aaf1ccde8fee446b12c1ab9a65a9992b41c7da`). Machine-readable outputs: `manuscript_quantitative_claims_audit_2026-09-24.json` and `reporting_checklist_locator_audit_2026-09-23.json`.
