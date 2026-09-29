# Title/abstract queue metadata audit

Audit date: 2026-09-24. This is a read-only completeness check of the generated title/abstract screening queue; no screening decisions were entered.

## Current queue

The queue contains 56,092 retained records: 54,940 from `frailty_neighborhood` and 1,152 from `primary_sleep_frailty`. The title field is populated for all 56,092 records. Abstract text is absent for 935 records (1.67%): 918 from `frailty_neighborhood` and 17 from `primary_sleep_frailty`. Both reviewer decision columns remain blank for all records. The queue has zero full-text-routed records.

The current queue is 112,139,734 bytes with SHA-256 `d29d6ec4a2d6b14e1207c5a6ab81f19fc108cc164cfc90e0e2b77db1e03b3ab1`, matching the byte count and digest in `review/record_build_outputs.tsv`. This confirms identity with the recorded review-build output; it does not validate record eligibility or make screening decisions.

## Handling implication

Absence of an abstract is missing information, not evidence that eligibility criteria are unmet. Reviewers should screen the available title, use the protocol's unclear/awaiting-classification route when relevance cannot be resolved, and seek the full text where indicated. Do not exclude a record solely because its abstract is absent. Preserve both independent decisions and adjudication in the registered queue.

## Source-XML abstract availability audit

The 935 empty-abstract PMIDs were checked against every occurrence in the 311 retained PubMed XML source files. All 935 occur in the source XML; 53 PMIDs have duplicate source occurrences, and none of those duplicates contains abstract text. Thus the empty fields were not caused by the record builder retaining an abstract-free occurrence over a duplicate that had an abstract. This establishes only that PubMed's retained source snapshot has no abstract text for these records; it says nothing about abstracts or full text in other databases or publisher pages.

Reproduction: `make -C frailty_paper audit-missing-abstract-sources`. The read-only audit verifies every XML source file's byte count and SHA-256 against `review/pubmed_source_files.tsv`, parses all 60,989 source records, and compares abstract availability for the 935 PMIDs. Machine-readable result: `analysis/missing_abstract_source_audit_2026-09-24.json` (queue SHA-256 `d29d6ec4a2d6b14e1207c5a6ab81f19fc108cc164cfc90e0e2b77db1e03b3ab1`; source-manifest SHA-256 `7191f09ab527eb1d8cf2a8a89a6bad5e4b7209daeca3955bc832f4fa40df45ed`). This audit did not change the queue or enter screening decisions.

## Reproduction

Read `review/screening/title_abstract_queue.tsv` as UTF-8 tab-separated data. Count an abstract as missing when its field is empty after whitespace trimming; count title presence and reviewer decisions by the corresponding non-empty fields. Source queue identity is recorded in `review/record_build_outputs.tsv`.
