# Full-family LAVA receipt integrity and frozen failure-gate audit — 2026-09-26 05:44 UTC

Audit window: 2026-09-26T05:43:59+00:00 to 2026-09-26T05:44:03+00:00. The run continued during this read-only scan.

## Receipt integrity

The scan verified **9,960** receipt files and matching nonempty logs. File count was unchanged during the scan; runner state was 9,959 at both scan start and end, then caught up to 9,960 at 05:44:14 UTC. Structural/identity/log issues: **0**; duplicate receipt identities: **0**. Expected family size: 29,940 (12 traits × 2,495 loci). Receipt hash stream: `32b7088f0095aa94de9123541507bb4fcafd2906419fb80e8581f5e9e90b0205`.

Two active claims were process-verified alive (PIDs 78548, 78549). The coordinator requested 6 workers and is running 2 under its current resource fallback. Launch failures: 0; stale recoveries: 0. Analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Frozen 1% locus-failure gate

The complete predicate fails for all 12 traits. `PROCESS_FAILED` receipts total **1,912**: 1,670 all-phenotype negative-variance and 242 no-reference-SNP outcomes. Per-trait lower bounds include these failures and non-tested univariate statuses.

| Sleep trait | PROCESS_FAILED | Failure lower-bound loci | Lower-bound fraction | Passes frozen 1% gate? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1,786 | 71.58% | No |
| chronotype | 304 | 1,326 | 53.15% | No |
| insomnia | 80 | 298 | 11.94% | No |
| longsleep | 118 | 362 | 14.51% | No |
| napping | 86 | 299 | 11.98% | No |
| shortsleep | 106 | 328 | 13.15% | No |
| sleep_apnea | 46 | 257 | 10.30% | No |
| sleep_efficiency | 147 | 383 | 15.35% | No |
| sleep_timing | 147 | 387 | 15.51% | No |
| sleepdur | 85 | 317 | 12.71% | No |
| sleepiness | 111 | 351 | 14.07% | No |
| snoring | 95 | 338 | 13.55% | No |

No threshold was relaxed and no receipt was rerun, reclassified, collated, or overwritten by this audit. Continue the existing campaign, but do not interpret local sharing or launch downstream PLACO/fine-mapping/colocalization under the failed frozen gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0543.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py`.
