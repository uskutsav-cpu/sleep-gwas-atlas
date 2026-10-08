# Literature evidence acquisition and release boundary

Execution date:2026-10-08. This document describes local evidence recovery; it does not grant redistribution rights. Retained third-party article bodies, publisher HTML, full supplements, extracted complete tables, PDF text, and raw search/web responses are **local review cache**, excluded from the new public evidence release. Selected factual prior-estimate crosswalk rows, bibliographic source URLs/hashes, technical reports and reproducible audit code are the intended tracked outputs. Historical source files already in repository history remain historical artifacts; this execution does not authorize relicensing them or republish new copies.

`literature_acquisition_manifest.tsv` lists26 original artifacts with path, URL/recovery entry point, SHA256, byte count, acquisition status and required-input status. `literature_access_ledger.tsv` duplicates those verified acquisition dispositions for readable review. `literature_source_hashes.tsv` also hashes local derived/cache artifacts; a hash authenticates bytes, not permission or scientific validity. Failed OUP/PMC/author-site responses are labeled failed and cannot be used as primary full text. Two unpreserved exact failed/widget request URLs are explicitly recovery entry points, not claimed original fetch URLs.

## Recreate a lawful local cache

Run from repository root, with curl and Python/openpyxl/pypdf available. **Do not bypass an access control or accept a changed hash silently.** Required current publisher XLSX/PDF URLs use official article supplement links, not third-party mirrors. Historical Morrison/Goodman/Dashti workbooks are obtained via their official landing pages because exact binary request URLs were not retained. Download the named original supplements using authorized access and put them at the manifest paths; verify the SHA256 before using them. Missing manual inputs cause a clear error.

```sh
python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py --acquire
python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py
```

`--acquire` verifies existing required bytes or downloads the current required files, fails on changed content, and regenerates Sun2022 text using pypdf with an expected text hash. It does not fetch HTML/challenge pages or unnecessary optional imaging supplements. Full regeneration also needs the historical extracted Morrison TSVs, prior phenotype inventory/prescreen and original1200/603/217 frozen input tables already identified by the script; it does not reconstruct native GWAS. Retrieval may require manual publisher access and new metadata review. Dependencies used here: Python3.12.14, openpyxl3.1.5, pypdf6.10.0 (verify the bundled runtime if replaying; version drift may change PDF text).

For source-free CI/release verification (Python standard library only):

```sh
python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py --verify-canonical
```

This verifies the retained1200/23 canonical universe, selected crosswalk comparator counts and zero certified novelty. It **does not independently reconstruct or verify publisher claims** when primary cache is absent. Human primary-source recovery/review remains necessary for a new literature audit date. Later searches require new dated query receipts; current indexing completeness is not guaranteed.

## Exact new-cache exclusions for repository owner

Add the following rules to repository `.gitignore` before staging. Existing historical `discovery_extension/provenance/prior_screens/` files are not deleted or edited. These rules intentionally exclude raw responses even when named receipts, because their payload includes publisher source text and third-party search snippets.

```gitignore
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.pdf
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.xlsx
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.html
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.xml
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.txt
/discovery_extension/sleep_submission_evidence_v1/sources/literature*.json
/discovery_extension/sleep_submission_evidence_v1/sources/literature*_extracted.tsv
/discovery_extension/sleep_submission_evidence_v1/sources/journal*.html
/discovery_extension/sleep_submission_evidence_v1/sources/journal*.pdf
/discovery_extension/sleep_submission_evidence_v1/sources/journal*.json
```

The access/acquisition/hash TSVs and this technical Markdown remain trackable. The selected115-row `tables/novelty_crosswalk.tsv` is a factual comparator dataset, not redistribution of the full publisher tables (8,418/4,046/1,724/257/232 rows). Nature Sleep Chart states CC BY-NC-ND4.0: no inferred permission to adapt/relicense publisher figures. Other publishers/checklist licenses must be reviewed directly before any future reproduction. No third-party full article/supplement is newly released as part of these intended outputs.
