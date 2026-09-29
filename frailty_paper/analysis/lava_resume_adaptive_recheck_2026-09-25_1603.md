# FI×sleep LAVA adaptive-resume checkpoint — 2026-09-25 16:03 UTC

## Safe pause and reconciliation

The live two-worker coordinator received pause request `pause_20260925T1546Z_requested_four_worker_restart` at 15:47:12 UTC. Both active loci finished, both workers emitted `worker_finished`, the coordinator state showed zero active jobs/workers, and `runner.lock` disappeared at 15:47:43 UTC. The campaign paused at 8,483/29,940 receipts.

The exact slot audit found 8,483 receipts and 21,457 pending slots. The initial hold-generation audit exposed an off-by-one pair-order value in the temporary pause-hold metadata; all 21,457 hold records were corrected before resumption. Full verification then passed for all 8,483 receipt identities/statuses and every pending hold’s trait, pair order, locus, and pause ID. Total reconciliation was exactly 29,940; no normal claim, active job, or runner lock remained. The frozen analysis-lock SHA-256 is `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256 is `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

## Adaptive resume

Resumed with the user-requested `--workers 4`. Startup’s swap guard reduced 4→3 at 6,649/7,168 MiB; the coordinator reduced 3→2 at 5,856/7,168 MiB through the absolute-use deep-pressure guard. Worker 3 completed its current locus and emitted `worker_finished` at 16:02:58 UTC. At 16:03:41 UTC the runner reported 8,486/29,940 receipts, two effective/live workers, two active claims, zero launch failures and zero stale recoveries. The three new receipt events since resume have unique keys, valid schema/identity/status and nonempty logs; no errors.

## Host pressure

At approximately 16:03 UTC, load averages were 26.15/25.16/27.12 on eight cores, system-wide free memory was 36%, and swap use was 7,775.88/8,192 MiB (94.9%). Continue at two workers while deep swap pressure persists. No stable completion estimate.
