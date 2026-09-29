# Systematic review and evidence map protocol

**Protocol version:** 1.0  
**Frozen:** 2026-09-22, before title/abstract screening  
**Registration:** None recorded. Do not imply prospective registration.  
**Search history:** PubMed XML acquisition and the three searches described in `pubmed/search_log.tsv` were completed on 2026-09-22 before this protocol file was frozen. No screening or study inclusion had been performed at protocol freeze. This timing will be reported transparently.

## Questions

**Primary review question:** Among adults, what observational and genetic evidence links sleep or circadian phenotypes with frailty?

**Secondary evidence-map question:** What evidence connects sleep/circadian phenotypes with physical frailty components, sarcopenia, strength, mobility, physical activity, cognition, disability, falls, mortality, healthspan or longevity?

## Eligibility criteria

- **Population:** Human adults aged 18 years or older. Mixed-age studies are eligible when adult results are separable or the study population is predominantly adult and the age information permits a defensible adult interpretation; otherwise classify as unclear and seek clarification/full text.
- **Exposure:** Sleep duration, insomnia, sleep quality/fragmentation/efficiency, chronotype/circadian timing, daytime sleepiness, napping, sleep-disordered breathing/apnea, snoring, or objectively measured sleep. Sleep medication alone is not a sleep phenotype unless the study also measures sleep.
- **Primary outcome:** A named or otherwise operationalized frailty construct, including Frailty Index, Fried phenotype/score, Clinical Frailty Scale, validated frailty scale, or a clearly described frailty classification.
- **Secondary outcomes:** Frailty components and the aging-related outcomes listed in the secondary question. These are evidence-map outcomes and are not interchangeable with frailty.
- **Study designs:** Human observational cross-sectional, cohort, case-control and longitudinal studies; GWAS and other human genetic association studies; genetic correlation, Mendelian randomization, fine-mapping or molecular-QTL studies when sleep and frailty/aging relevance is identifiable. Reviews may be used to find primary studies but will not be counted as primary evidence.
- **Exclusions:** Non-human or in-vitro-only studies; editorials, protocols and abstracts without usable results (retain as background leads where useful); studies without a sleep/circadian exposure or eligible outcome; pediatric-only studies; studies where relevant results cannot be separated and adult/frailty relevance cannot be established.
- **Language:** No language exclusion at search. If translation is unavailable, record as awaiting classification rather than excluding for results.

## Outcomes and extraction

For the primary synthesis, extract the sleep measurement, frailty definition, design, population/cohort, sample size, age, sex, ancestry, effect estimate and uncertainty, adjustment set, follow-up, direction, and authors' conclusion. For genetic studies additionally extract accession/provenance, genome build, ancestry, sample size/cases/controls, method, sample overlap, LD reference, genomic result, QTL tissue/cell context, and the level of evidence actually supported. Use `extraction_template.tsv`; leave unavailable values blank and document why in notes. Never infer missing sample sizes, accessions or results.

## Screening

1. Merge database exports with the PubMed record set while retaining each source database and source-record identifier.
2. Deduplicate using exact PMID, then exact normalized DOI, then normalized title; preserve all duplicate sources and the selected retained record in the audit trail.
3. Two reviewers should independently screen title/abstract and full text. Resolve disagreement by discussion or a third reviewer; record both decisions and the resolution.
4. Record full-text exclusion reasons using a mutually exclusive primary reason, with a free-text note if needed.
5. Do not mark the review complete until all supplied database exports have been deduplicated and screened and outstanding full texts are resolved or explicitly listed as unavailable.

## Risk of bias

Use the current JBI critical-appraisal checklist appropriate to each observational design. Record item-level judgments as Yes / No / Unclear / Not applicable, with a short evidence note and study/page locator. Missing information is Unclear, not low risk. For genetic association reporting, separately audit relevant STREGA/STROBE items; this is not a substitute for observational risk-of-bias assessment. Preserve two independent appraisals and adjudication.

## Synthesis

Describe observational and genetic evidence separately. Stratify narrative synthesis by sleep measurement mode/phenotype and frailty definition/design. Present secondary aging phenotypes as adjacent context, not as frailty equivalents. Do not pool estimates unless design, phenotype, outcome, effect scale and population are sufficiently comparable; any meta-analysis proposal must be documented before examining pooled results. Report heterogeneity and inconsistent/null findings. Genetic correlation is not causation; MR is interpreted only under its assumptions; colocalization, fine-mapping and QTL overlap are not proof of a causal gene or mechanism.

## Planned subgroup and sensitivity summaries

Where the evidence supports comparison, distinguish self-report from objective sleep, diagnosed from questionnaire phenotypes, FI from Fried/other frailty definitions, cross-sectional from longitudinal design, and exact from comparable phenotype replication. For genetic evidence, report ancestry, cohort overlap, independent replication availability, and method-specific QC. Do not run underpowered or unplanned subgroup tests as confirmatory analyses.

## Search sources and dates

PubMed searches and exact strings are in `pubmed/search_log.tsv`; raw XML and the generated source/checksum inventory are retained in `pubmed/`. Prepare manual exports for Embase, Scopus, Web of Science, and PsycINFO if available using `manual_search_instructions.md`. Search dates, interfaces, full strategies, result counts, deduplication and export dates must be recorded per database.

## Amendments

Operational and search-history amendments are recorded in `protocol_amendments.tsv`. As of 2026-09-24, the two recorded amendments do not change eligibility, outcomes, appraisal or synthesis methods. AMEND-002 registers one publisher-found post-snapshot candidate outside the frozen PubMed queue; it is unscreened and is not counted as included evidence. Record any later change with date, reason, screening status and affected files. Do not silently change eligibility, outcomes, appraisal or synthesis methods.
