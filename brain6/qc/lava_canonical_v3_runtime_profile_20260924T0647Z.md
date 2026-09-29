# Canonical LAVA v3 runtime and throughput snapshot

Snapshot time: **2026-09-24 06:47 UTC**
Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`
Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. This is a read-only profile; the run locks, receipts, checkpoints, and output root were not modified.

## Receipt-bound progress

The verifier accepted **2,053/2,495 locus receipts** (14,371/17,465 canonical cells), with **0 invalid receipts**. Statuses are **11,487 TESTED**, **2,884 NOT_RUN**, and **0 FAILED**. The frozen untested-cell ceiling is 873, exceeded by 2,011; family QC is guaranteed to fail. Continue the active run for complete accounting only; this profile is not a family decision and does not authorize pairwise promotion.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| ADHD | 1,900 | 153 | 0 |
| Bipolar disorder | 1,952 | 101 | 0 |
| Insomnia | 1,531 | 522 | 0 |
| Long sleep | 987 | 1,066 | 0 |
| Major depressive disorder | 1,598 | 455 | 0 |
| Parkinson's disease | 1,541 | 512 | 0 |
| Schizophrenia | 1,978 | 75 | 0 |

## Measured runtime

Metrics were joined only through accepted receipts to each receipt's worker config and sibling log. All **2,053** accepted loci had exactly one matching compute metric; no metric was missing or duplicated.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 48.56 s |
| Median locus compute wall time | 39.03 s |
| 90th percentile locus compute wall time | 82.65 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.669 s |
| Receipt throughput, 15-minute window | 488/hour |
| Receipt throughput, 60-minute window | 442/hour |
| Receipt throughput, 120-minute window | 430.5/hour |
| Remaining loci | 442 |

The corresponding arithmetic finish projections are 0.91, 1.00, and 1.03 hours. These projections are descriptive and workload-sensitive, not completion guarantees. The frozen 4-worker persistent process/batch setting remains unchanged.

## Reproduction

From the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```
