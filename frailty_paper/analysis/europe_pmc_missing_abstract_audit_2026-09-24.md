# Europe PMC missing-abstract metadata audit

Audit date: 2026-09-24. This supplemental lookup checked public Europe PMC records for each PMID whose abstract field is blank in the frozen title/abstract screening queue. It did not modify that queue, screening decisions, or the earlier PubMed-source audit.

## Result

All 935 PMIDs matched Europe PMC using `EXT_ID` queries. None had abstract body text. One Europe PMC record (PMID 35187946) had a 21-character structured field containing only the heading `Introduction`; it had no abstract body and was classified as unavailable. No abstract text was persisted from API responses. These findings do not imply ineligibility; reviewers must continue to use the protocol's unclear/awaiting-classification route where appropriate and seek full text when indicated.

## Reproduction and provenance

Run from the repository root:

```bash
frailty_paper/.venv/bin/python frailty_paper/scripts/audit_europe_pmc_missing_abstracts.py
```

The script reads `review/screening/title_abstract_queue.tsv`, queries the [Europe PMC REST search API](https://www.ebi.ac.uk/europepmc/webservices/rest/search) in batches of 50 using `EXT_ID:<PMID>`, and requests `resultType=core`. It counts non-heading text only; HTML heading elements are excluded so a record containing a heading with no body is not treated as an available abstract. Abstract content is processed transiently and never written. The TSV retains only PMID, match status, field/body character counts, body availability, and Europe PMC source. The JSON records the source queue checksum, query parameters, timestamp, counts, and TSV checksum.

Machine-readable outputs:

- `review/europe_pmc_missing_abstract_audit.tsv` (935 rows; SHA-256 `364a87d7f068384eea33f1a0d581b9d1de0575f724cda3fc1ffb075bc682117c`)
- `analysis/europe_pmc_missing_abstract_audit_2026-09-24.json`

The input screening queue remained byte-identical to its build manifest (SHA-256 `d29d6ec4a2d6b14e1207c5a6ab81f19fc108cc164cfc90e0e2b77db1e03b3ab1`). This is an availability check only; it does not screen or classify records and it does not replace the required manual database searches.
