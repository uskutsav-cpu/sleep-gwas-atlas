# FI×sleep LAVA full-family audit — 2026-09-26 11:13 UTC

The full read-only auditor validated 13,228/29,940 receipt files from 11:13:20 to 11:13:36 UTC. Inventory grew to 13,231 during the scan; runner state advanced 13,225→13,230. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Workers 78548 and 78549 passed liveness checks on distinct claims: insomnia locus 728 and snoring locus 1160. Launch failures and stale-claim recoveries were zero.

Six workers remain requested and two effective under the persisted swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged.

From the 10:57:28 audit (13,004 validated) to this scan, observed throughput was approximately 847 receipts/hour. At the 11:13:36 scan, 16,712 remained, corresponding to about 19.7 hours at the observed rate; this is a provisional estimate. Process failures total 2,545 (2,303 all-phenotype negative-variance; 242 no-reference-SNP). All 12 frozen 1% gates fail, so local-sharing inference and PLACO remain inadmissible. No scientific inputs, settings, thresholds, or receipts changed. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1114.json`.
