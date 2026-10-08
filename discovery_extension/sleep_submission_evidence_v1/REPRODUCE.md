# Reproduction boundary and commands

Run from the repository root using Python 3.11. The source-free numerical audit uses only the Python standard library. Figures use the exact pinned dependencies in `requirements/figures.txt`. Keep raw restricted data outside Git.

```sh
python3.11 -m venv work/evidence-env
task_python=work/evidence-env/bin/python
$task_python -m pip install -r discovery_extension/sleep_submission_evidence_v1/requirements/figures.txt -r discovery_extension/sleep_submission_evidence_v1/requirements/tooling.txt
# Verify delivered bytes before regeneration changes outputs or execution logs.
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/finalize_package.py --verify
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/verify_statistics.py
$task_python -m unittest discover -s discovery_extension/sleep_submission_evidence_v1/tests -v
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/build_figures.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/independent_statistics.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py --verify-canonical
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/build_support.py
$task_python discovery_extension/synthetic/test_extension_ldsc_collation.py
$task_python discovery_extension/sleep_submission_evidence_v1/scripts/validate_package.py
```

Optional source-dependent replay is separate: `audit_provenance.py` needs the recorded local historical worktrees/Git history and acquisition receipts; `audit_literature.py --acquire` reconstructs required official-source cache bytes, followed by `audit_literature.py` with `requirements/exports.txt` dependencies. This does not repeat every current database search. Do not treat `--verify-canonical` as fresh source verification. The provenance audit supports local historical worktree paths; inspect `--help` or its script before rerunning and preserve current acquisition receipts. The literature audit uses saved extracted primary tables and retains explicit unresolved classifications; live searches must be repeated for a later execution date. The workbook builder uses `@oai/artifact-tool` from Codex bundled dependencies and reads canonical TSV tables; editable XLSX is a static research-data export, not a replacement statistical engine.

Hash verification: `shasum -a 256 -c discovery_extension/sleep_submission_evidence_v1/hashes/PACKAGE_SHA256SUMS` from the repository root. Final receipt/manifests intentionally do not hash themselves. Current literature cache acquisition is `audit_literature.py --acquire`; source-free canonical evidence checking is `audit_literature.py --verify-canonical`. These are distinct evidence levels. For the workbook, run `workbook_input.py work/workbook_data.json`, expose the bundled artifact-tool modules in the script directory, then run `build_workbook.mjs work/workbook_data.json discovery_extension/sleep_submission_evidence_v1` with bundled Node. For editable Word tables, run `build_word_tables.py` with bundled Python/python-docx, then the Documents skill's `render_docx.py` for every output before delivery.

Scientific native rerun: BLOCKED. No original sleep/extension munged GWAS or EUR LD/weight reference is local. Acquire identities and hashes from original source manifests and licenses; use the frozen versioned streaming scripts only after adequate compute/network/storage feasibility is approved. Never insert regenerated files under missing frozen artifact names. The old full dense mirror estimates about 246.91 GiB and was not initiated.

Existing historical tests with missing scientific artifacts or absent R tooling are not a native scientific pass. Initial sparse-checkout failures are retained in logs; complete-checkout baseline and new-package results are separately recorded. Historical scientific and manuscript paths are unchanged.
