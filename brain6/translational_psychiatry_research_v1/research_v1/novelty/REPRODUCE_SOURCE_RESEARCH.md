# Source research reproducibility

All commands are run from the repository root. Public literature/source raw objects remain in ignored `work/tp-literature-20261007`; the three fixed FinnGen raw JSON objects are in `work/tp-finngen-fixed-20261007`. No source GWAS bulk archive was acquired by this research agent. The largest newly acquired supplementary object is Demontis MOESM6,9,081,487 bytes; individual object limits were frozen before acquisition.

Ordinary Python3 suffices for metadata acquisition, source-report accounting and clinical fixed-point calculations. XLSX derivation uses openpyxl in the bundled Python runtime. Exact openpyxl version and extraction-script hash are in `extraction_environment.json`. Text-only DOCX derivation uses zip/XML; supplemental PDF text was extracted with `pdftotext -layout` for cohort/source review. Native PDF/XLSX/DOCX/HTML objects, URLs, timestamps and SHA256 are retained in manifests/receipts. Cached-object reuse is explicit and preserves prior receipts when available; no inaccessible source is repeatedly requested by default.

```sh
python3 brain6/translational_psychiatry_research_v1/research_v1/novelty/acquire_and_screen.py search
python3 brain6/translational_psychiatry_research_v1/research_v1/novelty/acquire_and_screen.py extras
/Users/swethasunilkumar/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 brain6/translational_psychiatry_research_v1/research_v1/novelty/extract_competitor_signals.py
/Users/swethasunilkumar/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 brain6/translational_psychiatry_research_v1/research_v1/sources/derive_cohort_graph.py
python3 brain6/translational_psychiatry_research_v1/research_v1/sources/fixed_variant_reanalyze_v1.py
python3 brain6/translational_psychiatry_research_v1/research_v1/sources/build_source_reports.py
python3 brain6/translational_psychiatry_research_v1/research_v1/novelty/validate_source_accounting.py
```

The clinical manifest and acquisition code were frozen before outcomes at2026-10-07T23:04:02Z. Do **not** overwrite or rerun freeze: it deliberately refuses an existing manifest. The acquisition script also refuses silently reacquiring existing point outcomes; the initial parser holds are preserved. Deterministic reanalysis uses the explicit chr-schema correction, with identical hypothesis/allele/test choices. Parent independent raw-JSON checks are a distinct validator.

Signed publisher supplemental links may expire. The Jia URL was copied from the actual public article's supplement anchor; if a future uncached acquisition fails, obtain a freshly advertised public link from the publisher and record a new manifest/receipt, rather than editing/removing signature/access parameters. Restricted full text and bulk GWAS forms remain respected. Source hashes attest the analyzed historical object, not guaranteed future HTTP availability.

The wider prior cache `work/literature-20261007` is reused for Zu main Table3, selected official Ensembl37 mappings and prior GWAS Catalog metadata. It is not silently interpreted as a current participant-level cohort manifest. A clean checkout without retained raw sources can inspect derived results/manifests but cannot honestly reproduce empirical extraction; data-dependent execution should report that requirement rather than become a synthetic pass.
