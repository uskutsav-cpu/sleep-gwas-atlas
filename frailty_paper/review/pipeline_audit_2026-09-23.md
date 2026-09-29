# Systematic-review pipeline audit — 2026-09-23

## Frozen protocol and search

`protocol.md` is version 1.0, frozen on 2026-09-22 before title/abstract
screening. It defines adult population, sleep/circadian exposures, primary
frailty outcomes, secondary aging-map outcomes, study designs, exclusions,
extraction, two-reviewer screening, adjudication, JBI appraisal and narrative
synthesis. The protocol explicitly records that the search preceded protocol
freeze and that no prospective registration is recorded. Do not present it as
prospectively registered.

`pubmed/search_log.tsv` and the query files retain the three executed PubMed
searches and exact strategies. The acquisition reconciliation documents 60,989
records from PubMed only. The 311 retained XML batches freshly parsed; no
licensed manual database records are present. Embase, Scopus, Web of Science
and optional PsycINFO instructions and empty search-register/log templates
exist. No search was executed in those licensed services, and no export was
imported.

## Record and queue audit

| File | Data rows | Observed state |
|---|---:|---|
| `all_records.tsv` | 60,989 | Every PubMed record is represented; current file has a header and all source/dedup audit fields. |
| `deduplicated_records.tsv` | 56,092 | Each retained record includes duplicate count and source audit fields; 4,897 duplicate occurrences are linked. |
| `screening/title_abstract_queue.tsv` | 56,092 | Includes two independent reviewer decisions, reasons and adjudication fields; all decisions are blank. |
| `screening/fulltext_queue.tsv` | 0 | Header only; no title/abstract decisions have routed records to full text. |
| `screening/fulltext_exclusions.tsv` | 0 | Header only; no full-text exclusions. |
| `included_studies.tsv` | 0 | Header only; no studies included or extracted. |

`extraction_template.tsv` contains the planned observational and genetic
fields, including cohort/design, sample and age, ancestry, sleep/frailty
definitions, effect/uncertainty, adjustment, follow-up, genetic accession/build,
locus/variant/gene, tissue/cell type and risk-of-bias judgments. The empty
`included_studies.tsv` uses the corresponding extraction fields.

The JBI register lists tools by observational design, while genetic-association
reporting is routed to STREGA/STROBE as specified by the protocol. Item-level
appraisal and study-register files are header-only. No study has been appraised;
no study is classified as low risk. Before appraisal, archive the exact
design-specific JBI checklist file/version and checksum, as required by the
tool registry.

The repository review validator passed on 2026-09-23 and reported 56,092
title/abstract records, zero decisions, zero full-text records/decisions, zero
full-text exclusions and zero included studies. This validates queue
consistency; it does not constitute screening, extraction or review completion.

## Required next steps and limits

1. Run the authorized licensed searches, preserve each untouched export and
   exact search history, and populate the search log with platform/date/count
   and source checksums.
2. Import each export, validate source identifiers and cross-database duplicate
   links, rebuild queues without overwriting any entered reviewer decisions,
   and rerun the review validator.
3. Conduct independent title/abstract and full-text screening with documented
   adjudication; record full-text locations and mutually exclusive primary
   exclusion reasons.
4. Extract included-study fields and independently appraise applicable JBI
   items using archived current tools. Report missing information as unclear.
5. Update PRISMA values only from the reconciled imported-and-screened sources.

The review cannot be called complete while licensed exports are missing and
all 56,092 records remain unscreened. Search access is institution/user-side;
this audit makes no network query and does not scrape restricted services.

## Rebuild and provenance validation (2026-09-23 21:14 UTC)

Rebuilt the review artifacts from the retained 311 PubMed XML batches with the current `scripts/08_build_review_records.py` (SHA-256 `a21916abc8c47ac434a24d3ee5e7b0623082358363b9e026f8abcc5e03cc7f25`). Counts were unchanged: 60,989 source records, 56,092 retained records and 4,897 duplicate occurrences. All 10 generated artifacts matched the byte counts and SHA-256 values in `record_build_outputs.tsv`; that manifest hash matched `record_build_manifest.json`. The review validator passed and still reports zero screening decisions, zero full-text records, zero exclusions and zero included studies. This updates builder provenance only; it is not screening or appraisal.
