# Canonical LAVA v3 final runtime and throughput snapshot

Receipt snapshot: **2026-09-24 07:31 UTC**

Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`

Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. The frozen config and run-local lock both retain this setting. The coordinator completed under the same run ID and output root.

## Final receipt-bound family

The independent verifier accepted **2,495/2,495 locus receipts** and **17,465/17,465 canonical cells**, with zero invalid receipts. Final statuses: **13,745 TESTED**, **3,720 NOT_RUN**, **0 FAILED**. The frozen maximum is 873 untested cells (5%); the observed 3,720 is 21.30%, exceeding the maximum by 2,847. The immutable final decision is **`FAILED_QC_NOT_PROMOTED`**. Pairwise analysis and result promotion are not authorized.

| Trait | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|
| ADHD | 2,272 | 223 | 0 |
| Bipolar disorder | 2,327 | 168 | 0 |
| Insomnia | 1,824 | 671 | 0 |
| Long sleep | 1,204 | 1,291 | 0 |
| Major depressive disorder | 1,906 | 589 | 0 |
| Parkinson's disease | 1,859 | 636 | 0 |
| Schizophrenia | 2,353 | 142 | 0 |

## Runtime metrics

All 2,495 receipt-verified loci had exactly one matching compute metric; none were missing or duplicated. Metrics are joined through each accepted receipt to its worker configuration and sibling log.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 44.11 s |
| Median locus compute wall time | 35.54 s |
| 90th percentile locus compute wall time | 75.50 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.635 s |
| Receipt throughput, 15-minute window | 728/hour |
| Receipt throughput, 60-minute window | 571/hour |
| Receipt throughput, 120-minute window | 500/hour |

Throughput values describe the final rolling windows. No finish projection is applicable after completion.

## Immutable family decision

The fail-closed finalizer verified the complete 17,465-row aggregate against every receipt-bound cell and wrote `canonical_family_decision.json` in the run directory. Decision SHA256: `3569ce81d33479321f9b3adec7bd09c9b54a10a01f1fa4abd15888d9d2e66b22`. Aggregate SHA256: `ac450685970a35a4d0d06a9f02dd45ee5ac8f1af905ef411564b7e1f7f2f5352`. The decision records `qc_pass=false`, `promotion_permitted=false`, and `pairwise_stage_authorized=false`.

The coordinator's exit code 2 reflects the frozen family QC failure. It is not a license to alter the threshold or reinterpret underpowered cells as tested results.
