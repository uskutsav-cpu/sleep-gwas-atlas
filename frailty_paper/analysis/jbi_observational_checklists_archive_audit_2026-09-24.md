# JBI observational checklist archive audit

Audit time: 2026-09-24 02:56 UTC.

## Official source and archived tool files

The current [JBI Critical Appraisal Tools registry](https://jbi.global/critical-appraisal-tools) links the three design-specific files below. Each DOCX was fetched directly from its exact JBI-hosted URL through an approved public HTTPS route on 2026-09-24. No JBI account, form, alternate host, restricted content or access-control challenge was used.

| Tool | Official source | Archived file | Size (bytes) | SHA-256 | Identity check |
|---|---|---|---:|---|---|
| Analytical cross-sectional (2026) | [JBI DOCX](https://jbi.global/sites/default/files/2026-05/Assessment%20of%20Risk%20of%20Bias%20for%20Analytical%20Cross-sectional%20Studies%202026.docx) | `review/risk_of_bias/tools/jbi_analytical_cross_sectional_2026.docx` | 612,527 | `a37aaaacc15fe4bb128d9cbcfb16f978e4ced69213867973c0837852ddc0a2f9` | File title and printed citation identify Barker et al., *JBI Evidence Synthesis* 2026;24(3):401–408. |
| Cohort (2025) | [JBI DOCX](https://jbi.global/sites/default/files/2026-09/Checklist_for_Cohort_Studies%202025.docx) | `review/risk_of_bias/tools/jbi_cohort_checklist_2025.docx` | 1,304,359 | `21a30409b5b87f4dd208686bd1d869484bc173502c96bdca4a6ee44c94491e6c` | File title and printed citation identify Barker et al., *JBI Evidence Synthesis* 2025;23(3), Supplemental Digital Content 1, DOI 10.11124/JBIES-24-00103. See the earlier cohort-only audit for its detailed checks. |
| Case-control (current JBI registry link) | [JBI DOCX](https://jbi.global/sites/default/files/2026-05/Checklist_for_Case_Control_Studies.docx) | `review/risk_of_bias/tools/jbi_case_control_current.docx` | 701,047 | `57dc568ed410385d188f62d568ffc4753c9be005819b808f4ae87b33aac68378` | File title is “JBI Critical Appraisal Checklist for Case Control Studies”; its explanatory section cites the JBI Manual for Evidence Synthesis, Chapter 7 (2020). The registry does not assign a revision year, so none is inferred. |

All three archived copies match the downloaded-source SHA-256 values and pass DOCX/OOXML ZIP integrity checks. Text extraction confirmed the expected tool titles and item sets (8 analytical cross-sectional items, 11 cohort items, and 10 case-control items). The detailed cohort verification is recorded in `jbi_cohort_checklist_archive_audit_2026-09-24.md`.

## Render verification

The analytical cross-sectional file rendered to 10 pages and the case-control file to 6 pages using the bundled document renderer. Every page was visually inspected. The checklist forms, answer choices, explanatory questions and footer content are legible; no clipping, overlap, missing glyphs or broken checklist content was observed. The cohort DOCX rendered to 15 pages and was inspected in the preceding audit. Rendered PDFs and PNGs were temporary QA outputs; only the original DOCX source files are archived.

## Scope

The exact current JBI checklist files are now archived for all three observational designs routed by the review protocol. This closes the checklist-file acquisition blocker. No studies have yet been screened into the review or appraised; no risk-of-bias judgments were created. Recheck JBI at appraisal freeze and record the exact tool path and checksum in each eligible study's register row.
