# FI×sleep LAVA full-family audit — 2026-09-26 10:46 UTC

The full read-only auditor validated 12,857/29,940 receipt files from 10:46:35 to 10:46:40 UTC. Inventory reached 12,859 during the scan while runner state held at 12,857. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Workers 78548 and 78549 passed liveness checks on distinct claims: napping locus 710 and sleep_efficiency locus 1141. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective under the persisted swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged.

Since the 10:42:46 audit (12,804 validated receipts), the observed short-window rate is approximately 831 receipts/hour. About 17,083 remain, corresponding to ~20.6 hours at that rate; this is a provisional short-window estimate. Failure categories total 2,482: 2,240 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen 1% gates fail, so local-sharing inference remains inadmissible. No scientific inputs, settings, thresholds, or receipts were changed by this audit. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1045.json`.
