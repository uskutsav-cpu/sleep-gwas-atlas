# Reproduction and execution boundary

From the new branch root, with Python3.11+:

```
python brain6/confirmatory_v2/scripts/audit_evidence.py
python -m pytest -q brain6/confirmatory_v2/tests
python brain6/confirmatory_v2/scripts/run_contract_tests.py --suite source-free
PYTHONPATH=extensions/brain6 python brain6/confirmatory_v2/scripts/run_brain6_tests.py
python brain6/confirmatory_v2/scripts/validate_manuscript.py
python brain6/confirmatory_v2/scripts/build_figures.py
```

Continuation dependencies: numpy,pandas,scipy,matplotlib,pytest,beautifulsoup4; recorded actual versions in qc/environment.json. Production requirements remain separately pinned in original requirements files. No local native R/PLINK production run is claimed. Original archived-data integration: `python brain6/confirmatory_v2/scripts/run_contract_tests.py --suite integration`; full unchanged suite: `python -m unittest discover -s tests -v`. Full Brain6 integration requires restored sources: `PYTHONPATH=extensions/brain6 python -m pytest brain6/scripts/tests extensions/brain6/tests`.

Optional bounded public acquisition:

```
python brain6/confirmatory_v2/scripts/acquire_literature.py
python brain6/confirmatory_v2/scripts/reacquire_manifests.py --manifest brain6/confirmatory_v2/review/focused_acquisition_manifest.tsv
# Explicit opt-in network replay; public metadata only, ≤manifest byte cap/object:
python brain6/confirmatory_v2/scripts/reacquire_manifests.py --execute --manifest brain6/confirmatory_v2/review/reference_acquisition_manifest.tsv
python brain6/confirmatory_v2/scripts/diagnostic_extension.py
```

The coordinate/metadata extension requires local raw public literature cache acquired under its manifests. Without it, use shipped derived tables and receipts; do not fabricate missing responses. Raw responses are excluded from Git and handoff. Checkpoint receipts are atomically updated after each acquisition; successful objects have SHA256. Downloads are bounded30MB orlower; querytruncationdeclared. Inventory hashes stream1MB blocks; numerical audits process small tables and avoid loading retained large GWAS files. Figures are deterministic, seeded nowhere because no stochastic operation is run. Matplotlib usesAgg; fixed metadata date removed fromPDF to permit byte-stable figures. Native analyses remain blocked rather than launching expensive downloads.

PDF: edit the standaloneMANUSCRIPT/brain6_manuscript.tex in the Codex built-in editor; use compile_latex_document for diagnostics. Terminal export via installedTectonic: `tectonic --outdir brain6/confirmatory_v2/MANUSCRIPT brain6/confirmatory_v2/MANUSCRIPT/brain6_manuscript.tex`. Figures are supplied separately as vectorPDF/SVG with numerical source tables. Original source package remains unchanged; manifests bind both old/new source hashes. Do not run inventory.py on unrelated mutable checkouts expecting identicaltimestamps; inventories are session evidence.
