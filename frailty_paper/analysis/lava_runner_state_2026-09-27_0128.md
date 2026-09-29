# FI×sleep LAVA coordinator state — 2026-09-27T01:28:29Z

Read-only state snapshot; not a full receipt-integrity audit.

- Runner state: 18,835/29,940 receipts; 11,105 remain.
- Six workers requested; two effective. Coordinator and worker PIDs 72428, 72534 and 72535 were live.
- At the 01:28:29 UTC sample, worker 72534 held `insomnia` locus 1809 and worker 72535 held `snoring` locus 1808. Their active R child command lines matched the trait/locus claims; both claims were distinct.
- No launch failures or stale claims had been recorded in this invocation.
- The frozen lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Latest full receipt audit remains 01:16:19–01:16:29 UTC, which validated 18,735 files; two arrived during that scan, runner state reached 18,736, and receipt/claim/duplicate issues were zero. All 12 locked 1% gates failed.
- The 01:28:17–01:28:29 sample added two receipts in 12 seconds. This short interval is too brief to use as a throughput or completion estimate.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
