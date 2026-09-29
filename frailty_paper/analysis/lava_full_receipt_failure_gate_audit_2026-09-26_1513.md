# FI×sleep LAVA full-family audit — 2026-09-26 15:13 UTC

The read-only auditor validated 16,342/29,940 receipts (54.58%) from 2026-09-26T15:13:13+00:00 to 2026-09-26T15:13:27+00:00. Three files arrived during the scan (16,342 at start; 16,345 after); runner-verified state advanced 16,341→16,343. Receipt, claim, duplicate-claim, and duplicate trait/locus identity issue lists were empty. Workers 78548 and 78549 were live on distinct sleep-efficiency locus 1460 and sleep-timing locus 1460 claims. Launch failures and stale-claim recoveries were zero.

Six workers are requested and two effective under the persisted safeguard; swap use reached 2,772/3,072 MiB (90%), which reduces concurrency. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged; the runner lock matches. The audit is read-only and did not mutate run state or receipts.

Failure categories total 3,145: 2,903 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% trait gates fail; local-sharing, PLACO, and downstream inference remain inadmissible. No ETA is reported from short, variable receipt intervals. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1513.json`.
