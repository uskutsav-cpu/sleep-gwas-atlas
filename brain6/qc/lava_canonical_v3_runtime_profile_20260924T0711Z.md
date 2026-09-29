# Canonical LAVA v3 runtime and throughput snapshot

Snapshot time: **2026-09-24 07:11 UTC**

Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`

Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. The frozen execution config and run-local execution lock both specify four workers. This read-only profile did not modify the scheduler, locks, receipts, checkpoints, or output root.

## Receipt-bound progress

The verifier accepted **2,246/2,495 locus receipts** (15,722/17,465 canonical cells), with **0 invalid receipts**. Statuses are **12,535 TESTED**, **3,187 NOT_RUN**, and **0 FAILED**. The frozen untested-cell ceiling is 873, exceeded by 2,314; family QC is guaranteed to fail. Continue the active run for complete accounting only; this profile is not a family decision and does not authorize pairwise promotion.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| ADHD | 2,074 | 172 | 0 |
| Bipolar disorder | 2,132 | 114 | 0 |
| Insomnia | 1,657 | 589 | 0 |
| Long sleep | 1,079 | 1,167 | 0 |
| Major depressive disorder | 1,745 | 501 | 0 |
| Parkinson's disease | 1,689 | 557 | 0 |
| Schizophrenia | 2,159 | 87 | 0 |

## Measured runtime

Metrics were joined only through accepted receipts to each receipt's worker config and sibling log. All **2,246** accepted loci had exactly one matching compute metric; no metric was missing or duplicated.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 46.78 s |
| Median locus compute wall time | 37.41 s |
| 90th percentile locus compute wall time | 79.64 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.663 s |
| Receipt throughput, 15-minute window | 496/hour |
| Receipt throughput, 60-minute window | 477/hour |
| Receipt throughput, 120-minute window | 447/hour |
| Remaining loci | 249 |

The corresponding arithmetic finish projections are 0.50, 0.52, and 0.56 hours. These projections are descriptive and workload-sensitive, not completion guarantees. Four workers/five-locus batches remain the frozen policy.

## Reproduction

From the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```
