# Canonical LAVA v3 production runtime profile

Snapshot: 2026-09-24 03:31 UTC
Run: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`
Frozen schedule: **4 workers**, **5 loci per batch**

## Verified production state

The read-only receipt auditor verified **842/842 observed locus receipts**, covering 5,887 of 17,465 intended canonical cells, with zero invalid receipts. Statuses were 4,754 `TESTED`, 1,133 `NOT_RUN` (all `LOW_LOCAL_H2_UNDERPOWERED`), and zero `FAILED`. The untested lower bound exceeds the frozen maximum of 873 by 260; canonical family QC cannot pass. The run remained active and no full-family audit or family decision existed at this snapshot.

## Measured compute and throughput

Per-locus timings were read from each verified receipt's `worker_config_path` and its matching batch log. All 842 loci had exactly one matching metric record; none was missing or duplicated.

| Measurement | Observed |
|---|---:|
| Per-locus compute wall time, mean | 64.90 s |
| Per-locus compute wall time, median | 51.45 s |
| Per-locus compute wall time, 90th percentile | 122.09 s |
| Sum of per-locus worker compute time | 15.18 worker-hours |
| Median result-write time | 0.002 s |
| Median reference-load time per worker startup | 0.751 s |

Receipt modification times give these recent rolling rates:

| Window | Verified receipts | Receipts/hour | Arithmetic projection for 1,653 remaining loci |
|---|---:|---:|---:|
| 15 minutes | 25 | 100 | 16.53 hours |
| 60 minutes | 103 | 103 | 16.05 hours |
| 120 minutes | 330 | 165 | 10.02 hours |

Observed per-locus costs vary substantially, and the recent 60-minute receipt rate is below the earlier 50-locus benchmark's projected whole-run rate. The rolling projections are planning arithmetic, not finish-time guarantees. These measurements do not identify the cause of the throughput variation; no CPU-utilization trace was available. Reference load and result writing are small relative to measured locus compute wall time. The currently locked production run was not changed.

## Reproduce

Run the read-only profiler from the repository root:

```sh
python brain6/scripts/profile_lava_canonical_v3_runtime.py \
  --run-dir work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b \
  --input-root "/Volumes/Extreme SSD/brain6-work/preparation-v1/lava-inputs-v1" \
  --reference-provenance "/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/reference.provenance.json"
```

The tool verifies the frozen run and receipts, joins only their recorded worker configs to logs, writes JSON to standard output, and does not modify the production directory. Window projections update over time; this file is a historical snapshot.
