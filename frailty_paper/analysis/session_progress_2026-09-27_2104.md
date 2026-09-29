# Frailty project checkpoint — 2026-09-27 21:04 UTC

The coordinator had exited at 20:54 while the two registered workers continued recording receipt events. It was safely reattached at 21:03 with `--workers 2`; it reused PIDs 15810 and 15811 and launched no extra workers. Recovery details: `lava_coordinator_reattach_2026-09-27.md`. At 21:05:01 the runner reported 21,436/29,940 receipts, two active claims, and zero stale recoveries.

The independent full audit at 21:04:23–21:04:56 validated 21,416 receipts; inventory reached 21,419 and runner state 21,420 shortly after. Receipt/claim/duplicate checks were clean and both PIDs live. All 12 gates fail (3,726 negative-variance, 252 no-reference-SNP, two other process failures). Figure 4 was refreshed from this audit. The throughput window reset on coordinator restart; do not use the older mixed-worker rate for a new ETA.

The 56,117-record title/abstract queue remains undecided. Two distinct reviewer archives are ready. Licensed-source exports/full texts, human decisions/adjudication, exact HFRS object access and exact cohort intersections remain outstanding.
