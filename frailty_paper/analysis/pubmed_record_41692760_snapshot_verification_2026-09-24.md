# PMID 41692760 frozen PubMed snapshot verification

**Audit date:** 2026-09-24

**Conclusion:** PMID 41692760 was present in the frozen 2026-09-22 PubMed acquisition. It is already counted once in the deduplicated screening queue. No supplemental PRISMA record or count change is required.

## Source evidence

The PubMed source XML contains the record in three acquired query batches:

| Query | Batch | Bytes | SHA-256 | Record count |
|---|---|---:|---|---:|
| `frailty_neighborhood` | `pubmed/frailty_neighborhood/segment_2025-09-27_2026-09-22/batch_004200.xml` | 5,477,406 | `04d599ac0e746e8a33715eeb281cb4365038bed333df047b8f656997a2ad9966` | 200 |
| `genetic_sleep_frailty` | `pubmed/genetic_sleep_frailty/segment_1900-01-01_2026-09-22/batch_000000.xml` | 1,136,382 | `26acada96cc2bde31d848eb570234fcb4c721bff806f5afeb65917cf87bf9f30` | 36 |
| `primary_sleep_frailty` | `pubmed/primary_sleep_frailty/segment_1900-01-01_2026-09-22/batch_000000.xml` | 5,531,949 | `32c8c7f5906ea20e9faa7d2b2bef7776b035c862e4d224a04719b1a846fc40a0` | 200 |

The batch byte counts and SHA-256 values match `review/pubmed_source_files.tsv`. The registered query log (`review/pubmed/search_log.tsv`) gives the acquisition date and publication-date ceiling as 2026-09-22 for all three queries. Each XML occurrence carries the same PubMed history: Entrez 2026-02-15, PubMed 2026-02-16, and publication date 2026-02-16. Therefore, the paper was in the retained snapshot before the cutoff; it is not a post-snapshot index addition.

For this audit, each listed XML was parsed and its PubMedArticle with PMID
41692760 was checked for PubMed history dates. The saved source manifest was
then used to verify each file's byte count and SHA-256. The all-record,
deduplicated-record and screening TSVs were counted by exact PMID; the queue
row's two reviewer-decision fields and adjudication field were confirmed blank.

## Record-build and screening evidence

`review/all_records.tsv` has three rows for PMID 41692760, one for each query source. `review/deduplicated_records.tsv` has one row with `duplicate_count=3` and the three source paths above. `review/screening/title_abstract_queue.tsv` has one row with screening ID `primary_sleep_frailty:41692760:11f58fecd8d0`. Its two reviewer decisions and adjudication fields are blank, so the study remains unscreened. The priority TSV is an ordering aid and is not a screening decision.

## Correction to earlier notes

Prior versions of `BLOCKERS.md`, `review/pubmed_reconciliation.md`,
`analysis/prior_genetic_sleep_frailty_audit.md`, the manuscript, and the
reviewer audit incorrectly said the Huang et al. record was absent from the
frozen PubMed snapshot. Those statements are corrected in the current files.
The targeted follow-up found and verified a record that had already been
retrieved three times; it did not add a record or change any PRISMA count.
