# Canonical LAVA v3 runtime and throughput snapshot

Snapshot time: **2026-09-24 07:00 UTC**
Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`
Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. This is a read-only profile; the run locks, receipts, checkpoints, and output root were not modified.

## Receipt-bound progress

The verifier accepted **2,168/2,495 locus receipts** (15,176/17,465 canonical cells), with **0 invalid receipts**. Statuses are **12,109 TESTED**, **3,067 NOT_RUN**, and **0 FAILED**. The frozen untested-cell ceiling is 873, exceeded by 2,194; family QC is guaranteed to fail. Continue the active run for complete accounting only; this profile is not a family decision and does not authorize pairwise promotion.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| ADHD | 2,001 | 167 | 0 |
| Bipolar disorder | 2,060 | 108 | 0 |
| Insomnia | 1,607 | 561 | 0 |
| Long sleep | 1,041 | 1,127 | 0 |
| Major depressive disorder | 1,684 | 484 | 0 |
| Parkinson's disease | 1,632 | 536 | 0 |
| Schizophrenia | 2,084 | 84 | 0 |

## Measured runtime

Metrics were joined only through accepted receipts to each receipt's worker config and sibling log. All **2,168** accepted loci had exactly one matching compute metric; no metric was missing or duplicated.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 47.42 s |
| Median locus compute wall time | 38.12 s |
| 90th percentile locus compute wall time | 81.14 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.662 s |
| Receipt throughput, 15-minute window | 496/hour |
| Receipt throughput, 60-minute window | 468/hour |
| Receipt throughput, 120-minute window | 439.5/hour |
| Remaining loci | 327 |

The corresponding arithmetic finish projections are 0.66, 0.70, and 0.74 hours. These projections are descriptive and workload-sensitive, not completion guarantees. The frozen four-worker persistent process/batch setting remains unchanged.

## Reproduction

From the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```
