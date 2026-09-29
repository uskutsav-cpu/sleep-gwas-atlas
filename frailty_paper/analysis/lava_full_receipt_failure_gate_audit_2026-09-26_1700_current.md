# FI×sleep LAVA full-family audit — 2026-09-26 17:01 UTC

The read-only auditor validated 16,859/29,940 receipts (56.31%) from 2026-09-26T17:00:31+00:00 to 2026-09-26T17:01:21+00:00. Inventory advanced 16,859→16,862; runner state advanced 16,859→16,862. Receipt, claim, duplicate-claim, and duplicate trait/locus identity issue lists were empty. Workers 78548 and 78549 were live on distinct claims: napping locus 1512 and shortsleep locus 1512. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective because the existing safeguard detected swap use at 2,772/3,072 MiB (90%). Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged; runner lock hash matches. The audit is read-only and did not mutate run state or receipts.

From the 16:57:29 audit (16,841 validated) to this scan end, 18 receipts were added in 3m52s (~279/hour); 13,081 remain. The simple-rate ETA is ~46.9 hours, provisional because this short window is not representative. Failure categories total 3,238: 2,986 all-phenotype negative-variance and 252 no-reference-SNP. All 12 frozen 1% gates fail; local-sharing inference and PLACO remain inadmissible. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1700_current.json`.
