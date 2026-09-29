# Full-family LAVA receipt integrity and frozen failure-gate audit — 2026-09-26 05:30 UTC

Audit window: 2026-09-26T05:30:34+00:00 to 2026-09-26T05:30:40+00:00. The run continued during this read-only scan.

## Receipt integrity

The scan verified **9,909** receipt files and matching nonempty logs. File count was unchanged during the scan; runner state reported 9,909 receipts at scan start and 9,909 at scan end. Structural/identity/log issues: **0**; duplicate receipt identities: **0**. Expected family size: 29,940 (12 traits × 2,495 loci). Receipt hash stream: `31c94627322cdb5a1c62dad3a25145ebb1fb2beff0b3472880f013641d7b47a9`.

Two active claims were process-verified alive (PIDs 78548, 78549). The coordinator requested 6 workers and is running 2 under its current resource fallback. Launch failures: 0; stale recoveries: 0. Analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Frozen 1% locus-failure gate

The complete predicate fails for all 12 traits. `PROCESS_FAILED` receipts total **1,910**: 1,668 all-phenotype negative-variance and 242 no-reference-SNP outcomes. Per-trait lower bounds include these failures and non-tested univariate statuses.

| Sleep trait | PROCESS_FAILED | Failure lower-bound loci | Lower-bound fraction | Passes frozen 1% gate? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1,786 | 71.58% | No |
| chronotype | 304 | 1,326 | 53.15% | No |
| insomnia | 80 | 296 | 11.86% | No |
| longsleep | 118 | 359 | 14.39% | No |
| napping | 86 | 298 | 11.94% | No |
| shortsleep | 105 | 326 | 13.07% | No |
| sleep_apnea | 46 | 254 | 10.18% | No |
| sleep_efficiency | 146 | 382 | 15.31% | No |
| sleep_timing | 147 | 382 | 15.31% | No |
| sleepdur | 85 | 315 | 12.63% | No |
| sleepiness | 111 | 348 | 13.95% | No |
| snoring | 95 | 336 | 13.47% | No |

No threshold was relaxed and no receipt was rerun, reclassified, collated, or overwritten by this audit. Continue the existing campaign, but do not interpret local sharing or launch downstream PLACO/fine-mapping/colocalization under the failed frozen gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0530.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py`.
