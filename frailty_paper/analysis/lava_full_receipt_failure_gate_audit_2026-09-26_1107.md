# FI×sleep LAVA full-family audit — 2026-09-26 11:07 UTC

The read-only full-family auditor validated 13,153/29,940 receipt files from 11:07:08 to 11:07:14 UTC. Inventory was stable at 13,153 during the scan; runner state advanced 13,149→13,153. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Workers 78548 and 78549 passed PID liveness checks on distinct claims: sleep_apnea locus 1156 and sleepdur locus 724. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective under the persisted swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged.

From the 10:57:28 audit (13,004 validated), observed throughput was approximately 924 receipts/hour. At the 11:07:14 scan, 16,787 remained, corresponding to about 18.2 hours at that rate; this remains a provisional estimate. Failure categories total 2,530 (2,288 all-phenotype negative-variance; 242 no-reference-SNP). All 12 frozen 1% gates fail, so local-sharing inference and PLACO remain inadmissible. No scientific inputs, settings, thresholds, or receipts were changed by the audit. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1107.json`.
