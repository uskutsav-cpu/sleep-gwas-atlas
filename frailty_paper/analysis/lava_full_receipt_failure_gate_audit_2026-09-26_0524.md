# Full-family LAVA receipt integrity and frozen failure-gate audit — 2026-09-26 05:24 UTC

Audit window: 2026-09-26T05:24:45+00:00 to 2026-09-26T05:24:49+00:00. The run continued during this read-only scan.

## Receipt integrity

The scan verified **9,891** receipt files and matching nonempty logs. The file count was unchanged during the scan; runner state advanced from 9,890 to 9,891. Structural/identity/log issues: **0**; duplicate receipt identities: **0**. Expected family size is 29,940 (12 traits × 2,495 loci). Receipt hash stream: `8ac0aa53df9469ae80a477c69a7aaae2749e6bf93545f535a077a00aa0379c85`.

Two active claims were process-verified alive (PIDs 78548, 78549). The coordinator requested 6 workers and is running 2 effectively under the recorded swap fallback. Launch failures: 0; stale recoveries: 0. Frozen analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Frozen 1% locus-failure gate

The complete frozen predicate fails for all 12 traits. `PROCESS_FAILED` receipts total **1,910**: 1,668 all-phenotype negative-variance and 242 no-reference-SNP outcomes. The per-trait lower bound also includes non-tested univariate statuses.

| Sleep trait | PROCESS_FAILED | Failure lower-bound loci | Lower-bound fraction | Passes frozen 1% gate? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1,786 | 71.58% | No |
| chronotype | 304 | 1,326 | 53.15% | No |
| insomnia | 80 | 296 | 11.86% | No |
| longsleep | 118 | 359 | 14.39% | No |
| napping | 86 | 298 | 11.94% | No |
| shortsleep | 105 | 326 | 13.07% | No |
| sleep_apnea | 46 | 254 | 10.18% | No |
| sleep_efficiency | 146 | 381 | 15.27% | No |
| sleep_timing | 147 | 381 | 15.27% | No |
| sleepdur | 85 | 315 | 12.63% | No |
| sleepiness | 111 | 348 | 13.95% | No |
| snoring | 95 | 336 | 13.47% | No |

No threshold was relaxed and no receipt was rerun, reclassified, collated, or overwritten by this audit. Continue the existing atomic-claim campaign, but do not interpret local sharing or launch downstream PLACO/fine-mapping/colocalization under the failed gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0524.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py`.
