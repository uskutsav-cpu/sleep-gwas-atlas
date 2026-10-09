# Reproduction entry points

Run from repository root. The scripts write only the new research package or this task's new SSD directory. Historical original outputs are never execution destinations.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/00_recover_inventory.py
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/01_verify_native_inputs.py
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/02_recover_archived_extension.py --dry-run
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/02_recover_archived_extension.py
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/04_native_reproduction_runner.py --stage all
```

The archive reader scans a 56,985,877,290-byte compressed research archive sequentially but retains only receipt-pinned inputs/reference/QC/logs (about 1.2 GB). The input/source hash ledger distinguishes current hashes from historical expected matches. The official reference check retrieves the exact 33.4 MB archive from Zenodo and compares extracted member bytes; no alternative ancestry/reference is selected.

The native control used Python 3.9.23 in the recovered original LDSC environment and CBIIT revision `6c673952cee74bd5c57aef1555a03b1c015399a0`. `native_ldsc_capture.py` captures stock return values without changing estimation. The same numerical source files were compared against a fresh checkout of that commit. Dedicated library versions are recorded in the native receipt and are distinct from the workflow NumPy/Pandas pins.

After at least 3 GiB internal free space and all source hashes are admitted, execute one sequential stage:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/04_native_reproduction_runner.py --stage core --execute
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/04_native_reproduction_runner.py --stage extension --execute
PYTHONDONTWRITEBYTECODE=1 python3 sleep_unified_research_v1/scripts/04_native_reproduction_runner.py --stage validation --execute
```

These execute 45+100+13 h2 fits and 396+1,200+41 correlations in 190 jobs. A stock invocation or return code alone is insufficient; receipt cardinality, finite estimates, source/admissibility, exact family correction and original-versus-native comparison must all be audited. Existing sealed successful outputs may be skipped only after argument/output hash verification; failed/unsealed runs are preserved for explicit investigation.

This runner reproduces estimators from recovered processed inputs. Complete dense-source harmonization/munging reproduction remains a separate required chain. It does not secretly reconstruct missing sources, substitute newer data, admit a novel hypothesis, or convert the historic validation into independent sleep replication.

The bounded MVP acquisition is separately frozen and uses `scripts/03_acquire_admitted_mvp_source.py`. It is not retried automatically over existing partial/source bytes. Preprocessing, h2 eligibility and any related-phenotype independent pair tests must pass the frozen protocol before their outcomes are inspected.

Source/count audits are `06_audit_mvp_schema.py` and `09_audit_processed_source_counts.py`. They retain original source bytes and refuse to overwrite a completed diagnostic. The latter excludes macOS AppleDouble auxiliary entries and distinguishes the full HM3 template from rows with finite N/Z. `07_consolidate_inventory.py` combines the preserved initial inventory with separate exact/current/missing hash ledgers.

Figure generation requires Matplotlib 3.10.7 plus the versions recorded in `manifests/figure_runtime_versions_v1.json`. The task-created plotting runtime resides on the new SSD namespace. Set `MPLCONFIGDIR` and `TMPDIR` there and `PYTHONDONTWRITEBYTECODE=1`; do not install or repair dependencies inside an original archive. On exFAT, task-created AppleDouble style sidecars caused the first plot attempt to fail and were removed only from this new plotting environment; the cleanup receipt preserves the exact list. Run `08_make_evidence_figures.py` to produce seven figures. The 100-row extension atlas is intended for its vector PDF/SVG or full-resolution PNG.

Run the portable artifact checks from a clean checkout:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s sleep_unified_research_v1/tests -v
```

Tests use committed numerical tables, reviewer receipts and the pilot delete arrays; the SSD is unnecessary. Ten passing tests validate evidence invariants and guards, not raw-source reconstruction, complete native runs, phenotype equivalence, independence or novel biology. The final file manifest can be checked with `scripts/10_seal_package.py --verify`.

Prospective `FROZEN_NEW_ANALYSIS_PROTOCOL_v2.md` explicitly stops new clinical-insomnia pair tests while accession-specific phenotype coding remains uncertified. Historical core/extension/validation replay may resume independently after resources pass; the new-source phenotype gate does not change historical thresholds or sources.
