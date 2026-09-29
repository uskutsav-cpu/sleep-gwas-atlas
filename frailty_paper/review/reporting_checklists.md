# Reporting checklist status

This is a readiness scaffold, not a declaration of reporting-guideline compliance. The concise cross-framework status is recorded in `reporting_checklist_status.tsv`. Full item-level tables, including planned manuscript section and page/line fields, are generated under `../paper/reporting_checklists/` by running `python frailty_paper/scripts/14_build_reporting_checklists.py` from the repository root. Page/line remains `NOT_DRAFTED` until there is a manuscript.

PRISMA 2020 and PRISMA-S are applicable to the planned systematic review. PRISMA-S contains 16 search-reporting items; the current PubMed search snapshot is auditable, but database exports, full review search coverage and screening are incomplete. PRISMA 2020 contains 27 items; most result-dependent items remain pending. The protocol was frozen after the 2026-09-22 PubMed search snapshot, before screening; disclose this timing and do not describe the protocol as prospectively registered.

STREGA is applicable to a paper reporting genetic associations. Its 22 STROBE parent items, with genetic-association reporting requirements added where applicable, are scaffolded in `../paper/reporting_checklists/strega.tsv`; complete the checklist once the GWAS results and manuscript exist. MOOSE is currently conditional because no meta-analysis of observational effect estimates is planned. Reassess if that scope changes.

For study-level risk of bias, a design-specific JBI appraisal packet is now available in `risk_of_bias/`. It routes analytical cross-sectional studies to the revised 2026 JBI tool and cohort studies to the revised 2025 tool. The packet remains empty because no studies have been screened in; actual tool files must be versioned and checksummed before two independent appraisals. It defines no numeric score or overall low/moderate/high threshold.

Primary sources reviewed:

- PRISMA 2020 statement and 27-item checklist: https://www.bmj.com/content/372/bmj.n71
- PRISMA-S 16-item checklist and explanation: https://pmc.ncbi.nlm.nih.gov/articles/PMC8270366/
- STREGA statement: https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000022
- MOOSE statement: https://pubmed.ncbi.nlm.nih.gov/10789670/
- Current JBI design-specific appraisal tools and citations: https://jbi.global/critical-appraisal-tools (checked 2026-09-22)
