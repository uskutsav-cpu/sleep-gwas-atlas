# Full-family LAVA receipt integrity and frozen failure-gate audit — 2026-09-26 05:34 UTC

Audit window: 2026-09-26T05:34:39+00:00 to 2026-09-26T05:34:43+00:00. The run continued during this read-only scan.

## Receipt integrity

The scan verified **9,924** receipt files and matching nonempty logs. File count was unchanged during the scan; runner state reported 9,924 receipts at scan start and 9,924 at scan end. Structural/identity/log issues: **0**; duplicate receipt identities: **0**. Expected family size: 29,940 (12 traits × 2,495 loci). Receipt hash stream: `6d20c3e632e148db92a1c832c9674e142962a724be15def63793e40cea62e208`.

Two active claims were process-verified alive (PIDs 78548, 78549). The coordinator requested 6 workers and is running 2 under its current resource fallback. Launch failures: 0; stale recoveries: 0. Analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Frozen 1% locus-failure gate

The complete predicate fails for all 12 traits. `PROCESS_FAILED` receipts total **1,911**: 1,669 all-phenotype negative-variance and 242 no-reference-SNP outcomes. Per-trait lower bounds include these failures and non-tested univariate statuses.

| Sleep trait | PROCESS_FAILED | Failure lower-bound loci | Lower-bound fraction | Passes frozen 1% gate? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1,786 | 71.58% | No |
| chronotype | 304 | 1,326 | 53.15% | No |
| insomnia | 80 | 297 | 11.90% | No |
| longsleep | 118 | 361 | 14.47% | No |
| napping | 86 | 299 | 11.98% | No |
| shortsleep | 106 | 328 | 13.15% | No |
| sleep_apnea | 46 | 254 | 10.18% | No |
| sleep_efficiency | 146 | 382 | 15.31% | No |
| sleep_timing | 147 | 383 | 15.35% | No |
| sleepdur | 85 | 316 | 12.67% | No |
| sleepiness | 111 | 349 | 13.99% | No |
| snoring | 95 | 337 | 13.51% | No |

No threshold was relaxed and no receipt was rerun, reclassified, collated, or overwritten by this audit. Continue the existing campaign, but do not interpret local sharing or launch downstream PLACO/fine-mapping/colocalization under the failed frozen gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0534.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py`.
