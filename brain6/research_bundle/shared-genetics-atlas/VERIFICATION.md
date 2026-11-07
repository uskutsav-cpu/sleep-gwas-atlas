# Verification record, 27 October 2026

- Delivery-wide Markdown link audit: 34 local relative links checked, none missing after correcting the chromosome 5 case study's link.
- Existing atlas `manifest.tsv`: 52 recorded files checked for existence and SHA-256 agreement; no mismatches at audit time. This manifest records selected provenance-linked files, not the entire delivery.
- Expanded AlphaGenome prediction check: `python3 outputs/brain6-alphagenome-pilot/test_expanded.py` passed, verifying all 425 completed variant windows and their input binding.
- AlphaGenome pilot check: `BRAIN6_SOURCE_ROOT="$PWD/work/sleep-gwas-atlas" python3 outputs/brain6-alphagenome-pilot/test_pilot.py` passed 11 tests. The test explicitly checks frozen source membership, liftover references, checksums, selected prediction windows, AVI comparison, and local mocked API behavior. A deliberately failed mock query checks that failures remain visible; it is not a real service failure. Before adding the optional source-root setting, one test failed because the delivered copy's location differed from the source checkout. The check depends on the original checkout and external GRCh37-to-GRCh38 chain; it is not a standalone validation of biological predictions or the external service.
- `bundle_integrity.py` creates and verifies a broader snapshot, excluding live scan logs, partial outputs, bytecode caches and the snapshot itself. See `BUNDLE_INTEGRITY.md` for scope and limits.

The PLACO+ scan was still active at this audit. Scientific limitations, outstanding source data and the failed local-QC gates remain documented in `RESEARCH_ATLAS.md`; successful file checks do not resolve them.

## Later reproducibility repair

A subsequent review found three stale hashes in the selected 52-entry `manifest.tsv` after reports were updated. Those generated-report entries were refreshed; archived-source provenance fields were not changed. `python3 outputs/shared-genetics-atlas/verify_selected_manifest.py` now checks each listed file's size and SHA-256; it does not verify the original upstream study data or cover every delivered file.

Five missing Mak FI–physical-frailty calculation receipts, three companion logs and three scripts were restored as exact copies from the earlier workspace. Their [copy receipt](frailty_sensitivity/mak2023_cross_study/restored_source_receipt.json) records and verifies source and delivered SHA-256 hashes. These copies preserve the prior analyses; matching hashes do not independently validate their statistical assumptions or replace external GWAS inputs. The prior 34-link audit above was a bounded historical check. The newer `python3 outputs/shared-genetics-atlas/verify_local_links.py` checks simple inline Markdown local targets throughout this atlas, not anchors, images, web URLs, reference-style links or scientific validity.

The scan described as active in the historical record above has since terminated; its status and integrity incident are in the [completed PLACO audit](placo_diagnostic/completed_audit/RESULTS.md).
