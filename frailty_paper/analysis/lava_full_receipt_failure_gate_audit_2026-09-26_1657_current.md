# FI×sleep LAVA full-family audit — 2026-09-26 16:57 UTC

The read-only auditor validated 16,841/29,940 receipts (56.25%) from 2026-09-26T16:57:07+00:00 to 2026-09-26T16:57:29+00:00. Inventory advanced 16,841→16,843; runner state advanced 16,841→16,843. Receipt, claim, duplicate-claim, and duplicate trait/locus identity issue lists were empty. Workers 78548 and 78549 were live on distinct claims: shortsleep locus 1510 and sleep_apnea locus 1510. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective because the existing safeguard detected swap use at 2,772/3,072 MiB (90%). Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged; runner lock hash matches. The audit is read-only and did not mutate run state or receipts.

From the 16:52:28 audit (16,800 validated) to this scan end, 41 receipts were added in 5m01s (~490/hour); 13,099 remain. The simple-rate ETA is ~26.7 hours, provisional because this short window is not representative. Failure categories total 3,234: 2,982 all-phenotype negative-variance and 252 no-reference-SNP. All 12 frozen 1% gates fail; local-sharing inference and PLACO remain inadmissible. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1657_current.json`.
