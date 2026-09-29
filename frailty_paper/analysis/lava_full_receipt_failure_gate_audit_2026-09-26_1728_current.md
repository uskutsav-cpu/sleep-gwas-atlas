# FI×sleep LAVA full-family audit — 2026-09-26 17:30 UTC

The read-only auditor validated 16,973/29,940 receipts (56.69%) from 2026-09-26T17:29:10+00:00 to 2026-09-26T17:30:14+00:00. Inventory advanced 16,973→16,975; runner state advanced 16,972→16,975. Receipt, claim, duplicate-claim, and duplicate trait/locus identity issue lists were empty. Workers 78548 and 78549 were live on distinct claims: sleep_efficiency locus 1523 and sleep_timing locus 1523. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective because the existing safeguard detected swap use at 2,772/3,072 MiB (90%). Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged; runner lock hash matches. The audit is read-only and did not mutate run state or receipts.

From the 17:23:36 audit (16,934 validated) to this scan end, 39 receipts were added in 6m38s (~353/hour); 12,967 remain. The simple-rate ETA is ~36.7 hours, provisional because short-window throughput varies. Failure categories total 3,260: 3,008 all-phenotype negative-variance and 252 no-reference-SNP. All 12 frozen 1% gates fail; local-sharing inference and PLACO remain inadmissible. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1728_current.json`.
