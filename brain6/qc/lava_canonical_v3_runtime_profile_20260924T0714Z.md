# Canonical LAVA v3 runtime and throughput snapshot

Snapshot time: **2026-09-24 07:14 UTC**

Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`

Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. The frozen execution config and run-local execution lock both specify four workers. This read-only profile did not modify the scheduler, locks, receipts, checkpoints, or output root.

## Receipt-bound progress

The verifier accepted **2,286/2,495 locus receipts** (16,002/17,465 canonical cells), with **0 invalid receipts**. Statuses are **12,750 TESTED**, **3,252 NOT_RUN**, and **0 FAILED**. The frozen untested-cell ceiling is 873, exceeded by 2,379; family QC is guaranteed to fail. Continue the active run for complete accounting only; this profile is not a family decision and does not authorize pairwise promotion.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| ADHD | 2,108 | 178 | 0 |
| Bipolar disorder | 2,171 | 115 | 0 |
| Insomnia | 1,686 | 600 | 0 |
| Long sleep | 1,101 | 1,185 | 0 |
| Major depressive disorder | 1,773 | 513 | 0 |
| Parkinson's disease | 1,715 | 571 | 0 |
| Schizophrenia | 2,196 | 90 | 0 |

## Measured runtime

Metrics were joined only through accepted receipts to each receipt's worker config and sibling log. All **2,286** accepted loci had exactly one matching compute metric; no metric was missing or duplicated.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 46.41 s |
| Median locus compute wall time | 37.06 s |
| 90th percentile locus compute wall time | 78.46 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.663 s |
| Receipt throughput, 15-minute window | 500/hour |
| Receipt throughput, 60-minute window | 483/hour |
| Receipt throughput, 120-minute window | 452/hour |
| Remaining loci | 209 |

The corresponding arithmetic finish projections are 0.42, 0.43, and 0.46 hours. These projections are descriptive and workload-sensitive, not completion guarantees. Four workers/five-locus batches remain the frozen policy.

## Reproduction

From the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```
