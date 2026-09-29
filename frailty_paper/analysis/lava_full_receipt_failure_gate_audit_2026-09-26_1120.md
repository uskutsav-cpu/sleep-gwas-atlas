# FI×sleep LAVA full-family audit — 2026-09-26 11:20 UTC

The read-only full-family auditor validated 13,319/29,940 receipt files from 11:20:18 to 11:20:24 UTC. Inventory advanced to 13,320 during scanning while runner state held at 13,319. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Workers 78548 and 78549 passed PID liveness checks on distinct sleep_apnea loci 732 and 1165. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective under the existing swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged.

Compared with the 11:13:20 audit (13,228 validated), the observed rate is approximately 783 receipts/hour. At 11:20:24, 16,620 remained, giving a provisional ~21.2 hours at that rate. Process-failure categories total 2,557 (2,315 all-phenotype negative-variance; 242 no-reference-SNP); all 12 frozen 1% gates fail, so local-sharing inference and PLACO remain inadmissible. No scientific inputs, settings, thresholds, or receipts changed. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1120.json`.
