# FI×sleep LAVA full-family audit — 2026-09-26 15:00 UTC

The read-only auditor validated 16,173/29,940 receipts (54.02%) from 2026-09-26T15:00:17+00:00 to 2026-09-26T15:00:24+00:00. Two files arrived during the scan (16,173 at start; 16,175 after); runner-verified state stayed at 16,173. Receipt, claim, duplicate-claim, and duplicate trait/locus identity issue lists were empty. Workers 78548 and 78549 were live on distinct sleep-efficiency locus 1443 and sleep-timing locus 1443 claims. Launch failures and stale-claim recoveries were zero.

Six workers are requested and two effective under the persisted safeguard; swap use reached 2,772/3,072 MiB (90%), which reduces concurrency. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged; the runner lock matches. The audit is read-only and did not mutate run state or receipts.

Failure categories total 3,112: 2,870 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% trait gates fail; local-sharing, PLACO, and downstream inference remain inadmissible. No ETA is reported from this short, variable receipt interval. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1500.json`.
