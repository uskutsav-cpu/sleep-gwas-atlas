# FI×sleep LAVA full-family audit — 2026-09-26 09:25 UTC

The read-only full-family auditor validated 11,783/29,940 immutable receipt files from 09:25:22 to 09:25:27 UTC. The receipt inventory grew by one during the scan; the coordinator reported 11,784 at 09:25. Remaining at the coordinator snapshot: 18,156.

Receipt-integrity and claim-identity checks report zero issues. Two workers were live on distinct jobs. Launch failures and stale recoveries remain zero. The runner requests six total workers and is held to two effective by the adaptive swap guard at 2,772/3,072 MiB (90%). Frozen lock hash (`74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`) and input-manifest hash (`8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`) are unchanged.

From the 09:06 audit (11,565 validated at 09:06:43) to this scan (11,783 validated at 09:25:22), throughput was approximately 700 receipts/hour. At that observed rate, 18,156 remaining receipts imply about 26 hours. This estimate is provisional and depends on the current resource-constrained two-worker execution. CPU/per-process RAM and disk-I/O telemetry remain unavailable because `ps`, `sysctl`, and `iostat` are blocked by sandbox process-inspection restrictions.

The all-traits 1% failure gate remains failed: 2,014 all-phenotype negative-variance process outcomes and 242 no-reference-SNP outcomes; no local-sharing inference is admissible. The locked plan check passes (12 sleep traits, multiplicity 396), while exact participant-level cohort intersections remain unknown. Review validation and both 113-batch reviewer packet dry runs pass structurally, but no human decisions were supplied or applied. Licensed database exports, adjudication, full-text retrieval, exact custom HFRS summary statistics, and independent replication remain unresolved. Overall project readiness remains **NO-GO**.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0925.json`.
