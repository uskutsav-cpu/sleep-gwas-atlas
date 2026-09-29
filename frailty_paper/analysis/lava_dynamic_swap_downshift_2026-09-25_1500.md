# LAVA dynamic-swap downshift checkpoint — 2026-09-25 15:00 UTC

## Worker transitions

The user-authorized resume began from the independently verified 8,267-receipt checkpoint with four workers. The coordinator safely reduced 4→3 at 14:47:54 UTC when macOS reported 1,748/2,048 MiB swap (≥82%). It reduced 3→2 at 14:52:53 UTC at 4,842/5,120 MiB (≥90%). The host dynamically enlarged swap, so the original fraction-only guard could miss the second pressure rung and a concurrent pause update could be overwritten by stale control state. I requested safe pauses, allowed claimed loci to finish, and corrected both conditions in commit `39acb0b4986fc724b936b5d2ba811c55f0b657da`.

The updated scheduler uses percentage or absolute swap-use triggers (1,536 MiB for 6→4/4→3; 3,584 MiB for 3→2), rereads the latest concurrency control before appending fallback state, and preserves a newer pause request. The focused test suite passes 16/16, including dynamic denominator and pause-preservation regressions.

## Clean pause audit

After all in-flight loci finished, the runner exited with zero workers, zero active jobs, zero launch failures, zero stale recoveries, and no `runner.lock`. The new additive external checkpoint records 8,408/29,940 valid receipts (21,532 remaining). Independent verification confirms:

- All 8,267 carried-forward receipt hashes still match the prior manifest.
- All 141 receipts created during this resume pass identity/schema validation.
- Zero duplicate receipt keys, malformed event rows, failures, or receipt validation errors.
- The frozen analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- The prepared-input manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

Checkpoint: `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/paused_checkpoint_20260925T1500Z_deep_pressure_downshift.json` (SHA-256 `40a357178bf3cd1bc082d0f518d407675bf7e9fe8dd1384f38dff406003e1db1`). Receipt manifest: `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/receipt_manifest_20260925T1500Z_deep_pressure_downshift.tsv` (SHA-256 `4bc061fbca6b589a72fe351603a75220d8e727410e7ed0b05c27e466ef08690f`). The 2-worker resume command is `python3 frailty_paper/scripts/51_run_lava_sensitivity_parallel.py --resume-paused --workers 2`.

The campaign is not complete. Existing first-pair and chronotype QC failures remain under the frozen rules; local sharing and PLACO remain unjustified. The systematic review is unscreened and manual/source-access blockers remain as documented in the readiness report.
