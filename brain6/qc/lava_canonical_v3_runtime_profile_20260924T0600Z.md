# Canonical LAVA v3 runtime and throughput snapshot

Snapshot time: **2026-09-24 06:00:30 UTC**  
Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`  
Execution policy: **4 workers**, **5 loci per persistent chromosome batch**. This is a read-only profile; no run configuration or output was modified.

## Receipt-bound progress

The full verifier accepted **1,695/2,495 locus receipts** (11,865/17,465 canonical cells). There were **9,514 TESTED**, **2,351 NOT_RUN**, **0 FAILED**, and **0 invalid receipts**. The untested lower bound exceeds the frozen ceiling of 873 cells by 1,478, so the family QC gate is guaranteed to fail. The run remains active for complete accounting; this profile is not a family decision and authorizes no pairwise promotion.

## Measured runtime

Locus times were joined only from each accepted receipt's recorded worker config and its matching log; 1,695 observed receipts had one matching timing metric each.

| Metric | Observation |
|---|---:|
| Mean locus compute wall time | 52.26 s |
| Median locus compute wall time | 41.79 s |
| 90th percentile locus compute wall time | 93.51 s |
| Median result write time | 0.002 s |
| Median reference load per worker startup | 0.687 s |
| Receipt throughput, 15-minute window | 420/hour |
| Receipt throughput, 60-minute window | 412/hour |
| Receipt throughput, 120-minute window | 407/hour |
| Remaining loci | 800 |

The corresponding arithmetic finish projections are 1.90, 1.95, and 1.97 hours for the 15-, 60-, and 120-minute windows. These are workload-sensitive projections, not completion guarantees. The main efficiency gains remain the frozen persistent R workers, cached chromosome reference per batch, and five-locus batching; the measured reference-load and write times are small relative to compute time. The 4-worker policy was not changed.

## Hash-bound context

| Object | SHA256 |
|---|---|
| Canonical family lock | `f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417` |
| Execution lock | `f37a212ff44472004549095a3b17d28da0d130e2ab2a9ea507feb87312b65d27` |
| Run identity document | `bf6727b6cdd9109fe99efe4135f6ba092ca9ee9e5bcdee86a6dcb306df5295e8` |
| Canonical worker script | `79b66315996e855be0e25f11bc2cb54c59532a1d833c0228262094779f764d1b` |
| Runtime profiler script | `185715fd8df95db2749b2943057d771198c9a5b00d99553b9f85c5ac1fee196d` |
| Materialized input provenance | `c50dae720237b548a282d931118f61348821a139f198c64d6d6876a015ce3d02` |
| Official reference provenance | `35ab371935b9c260dab21be248f4907b693565787ef9f62334e7d8030204270e` |

## Reproduction

From the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```
