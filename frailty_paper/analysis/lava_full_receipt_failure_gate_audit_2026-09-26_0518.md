# Full-family LAVA receipt integrity and frozen failure-gate audit — 2026-09-26 05:18 UTC

Audit window: 2026-09-26T05:17:58+00:00 to 2026-09-26T05:18:03+00:00. The run continued during this read-only scan.

## Receipt integrity

The scan verified **9,870** receipt files and matching nonempty logs. The file count was unchanged during the scan; runner state reported 9,870 receipts at start and 9,870 at end. Structural/identity/log issues: **0**; duplicate receipt identities: **0**. Expected family size is 29,940 (12 traits × 2,495 loci). The receipt hash stream is `89957bd4db1ded9f62122f5d3d43ebc0f472798b2567becd1cc784687d082c51`.

Two active claims were process-verified alive (PIDs 78548, 78549): longsleep locus 989, sleepdur locus 563. The coordinator requested 6 workers and is running 2 effectively under the recorded swap fallback: `swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two`. Launch failures: 0; stale recoveries: 0. Frozen analysis lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Frozen 1% locus-failure gate

The full predicate already fails for all 12 traits. `PROCESS_FAILED` receipts alone total **1,910**: 1,668 all-phenotype negative-variance outcomes and 242 no-reference-SNP outcomes. The table reports both these receipt failures and the broader per-trait lower bound including non-tested univariate statuses.

| Sleep trait | PROCESS_FAILED | Failure lower-bound loci | Lower-bound fraction | Passes frozen 1% gate? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1,786 | 71.58% | No |
| chronotype | 304 | 1,326 | 53.15% | No |
| insomnia | 80 | 296 | 11.86% | No |
| longsleep | 118 | 359 | 14.39% | No |
| napping | 86 | 298 | 11.94% | No |
| shortsleep | 105 | 326 | 13.07% | No |
| sleep_apnea | 46 | 254 | 10.18% | No |
| sleep_efficiency | 146 | 379 | 15.19% | No |
| sleep_timing | 147 | 380 | 15.23% | No |
| sleepdur | 85 | 314 | 12.59% | No |
| sleepiness | 111 | 347 | 13.91% | No |
| snoring | 95 | 335 | 13.43% | No |

No threshold was relaxed and no receipt was rerun, reclassified, collated, or overwritten by this audit. Keep the current atomic-claim campaign running for planned receipt completion, but do not interpret local sharing or launch downstream PLACO/fine-mapping/colocalization under the failed frozen gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0520.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py` (SHA-256 `d064792f2afba1f8a9e2232fa78d4fae33a346205574382291809942abed4269`).
