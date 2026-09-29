# Risk-of-bias appraisal packet

The review protocol requires design-appropriate JBI appraisal for observational studies, two independent reviewers, item-level judgments, evidence notes and source locators. No study has been screened in or appraised yet; both TSV files are header-only templates and are not evidence of appraisal completion.

## Tool selection and version control

Use `tool_registry.tsv` to route analytical cross-sectional, cohort/longitudinal and case-control studies to the current JBI tool for that design. JBI revised its analytical cross-sectional tool in 2026 and its cohort tool in 2025. Before appraisals begin, retrieve the official file, record its exact version/date and SHA-256 in `study_register.tsv`, and keep that file with the review audit materials. Recheck JBI for updates at the appraisal freeze. Never mix item sets from different versions within a study or silently migrate an in-progress appraisal; document and assess any change as a protocol amendment.

For GWAS or other genetic association studies, do not force an observational JBI checklist onto the design. Record `NOT_APPLICABLE_GENETIC_ASSOCIATION` for the JBI tool field and separately assess applicable STREGA/STROBE reporting and the source-specific genetic QC, ancestry, overlap and replication evidence. A reporting checklist is not itself a risk-of-bias score.

## Current checklist file archive status (2026-09-24)

The current [official JBI Critical Appraisal Tools registry](https://jbi.global/critical-appraisal-tools) links to the three archived, design-specific files recorded in `tool_registry.tsv`:

- Analytical cross-sectional (2026): `tools/jbi_analytical_cross_sectional_2026.docx` (612,527 bytes; SHA-256 `a37aaaacc15fe4bb128d9cbcfb16f978e4ced69213867973c0837852ddc0a2f9`).
- Cohort (2025): `tools/jbi_cohort_checklist_2025.docx` (1,304,359 bytes; SHA-256 `21a30409b5b87f4dd208686bd1d869484bc173502c96bdca4a6ee44c94491e6c`).
- Case-control (current registry file; no revision year assigned by JBI): `tools/jbi_case_control_current.docx` (701,047 bytes; SHA-256 `57dc568ed410385d188f62d568ffc4753c9be005819b808f4ae87b33aac68378`).

All three files were downloaded from their official JBI-hosted URLs, passed DOCX/OOXML integrity validation and were rendered and visually inspected page-by-page. Details and printed citations are in `../../analysis/jbi_observational_checklists_archive_audit_2026-09-24.md`; cohort-specific verification is also in `../../analysis/jbi_cohort_checklist_archive_audit_2026-09-24.md`. Prior retrieval attempts remain in `tool_download_log.tsv`; the ordinary browser/web reader failed to display DOCX responses, but the exact public JBI HTTPS downloads succeeded. The archived source files were not edited or reconstructed.

No studies have yet been screened in or appraised. At appraisal freeze, recheck the registry and record the exact file and SHA-256 in each eligible study's `study_register.tsv` row. Do not begin appraisal judgments before eligible studies and full texts enter the dual-review workflow.

## Entry rules

- Add one row per study × tool item × reviewer to `appraisal_items.tsv`.
- The two reviewers work independently before comparing decisions. Use only `Yes`, `No`, `Unclear`, or `Not applicable` as `judgment`, following the frozen protocol.
- Missing or unreported information is `Unclear`; `Not applicable` requires a brief reason.
- Every judgment must include a concise evidence note and a page, section, table, figure, supplement or other stable source locator. Do not infer a favorable judgment from absent reporting.
- Keep disagreements visible in both original reviewer rows. Record resolution in `adjudicated_judgment`, `adjudicator_id`, and `adjudication_note`; do not overwrite the independent decisions.
- Do not calculate percent scores, collapse item judgments to “low/moderate/high,” or exclude a study based on an invented cutoff. The protocol does not define such a threshold.
- Preserve all included studies in the appraisal register, including studies for which full text or enough information remains unavailable; label the status and limitation instead of silently omitting them.

## Official sources

- [JBI Critical Appraisal Tools page](https://jbi.global/critical-appraisal-tools) (checked 2026-09-22; authoritative list of current design-specific tools and citations).
- [Revised JBI analytical cross-sectional tool](https://jbi.global/sites/default/files/2026-05/Assessment%20of%20Risk%20of%20Bias%20for%20Analytical%20Cross-sectional%20Studies%202026.docx).
- [Revised JBI cohort tool](https://jbi.global/sites/default/files/2026-09/Checklist_for_Cohort_Studies%202025.docx).
- [JBI case-control checklist](https://jbi.global/sites/default/files/2026-05/Checklist_for_Case_Control_Studies.docx).

The case-control tool is listed on the current JBI page, but the page does not identify a revision year in its citation block; therefore the exact downloaded file and checksum must be frozen before appraisal rather than assigning it an unsupported version date.
