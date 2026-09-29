# Brain6 canonical LAVA v3 scheduler benchmark

**Date:** 2026-09-23  
**Purpose:** validate the pair-independent 7 × 2,495 canonical-univariate worker architecture and select a resumable production schedule. These diagnostic subsets do not support a family decision or bivariate promotion.

## Locked analysis and runtime

- Family: `brain6-lava-canonical-v3`, exactly 17,465 planned trait × locus cells.
- Traits: ADHD, bipolar disorder, insomnia, long sleep, MDD, Parkinson disease, and schizophrenia.
- Gate: `p < 2.86286859433152e-06`; maximum untested fraction remains 5%.
- Runtime: R 4.3.3, LAVA 0.1.5, official UK Biobank LAVA v1.1 reference, roundoff-only block-matrix symmetrization.
- Both measured schedulers used the frozen **4-worker** setting, BLAS/OpenMP thread limits of one, atomic locus outputs, per-cell hashes, and input-bound receipts.
- The source-level insomnia/locus-2 check matched the public LAVA path: 1,260 SNPs, 179 components, observed h² `8.33079e-05`, latent h² `0.000147645`, and `p=0.0878159` in both paths.

## Measurements

| Schedule | Disjoint diagnostic sample | Wall time | Child CPU | Peak child RSS* | Receipt audit | Output footprint |
|---|---:|---:|---:|---:|---|---:|
| Four workers, 10 loci per chromosome batch | 50 loci / 350 cells, 5 chromosomes, 5 R batches | 1,245.4 s | 2,295.2 s | 764 MB | 50/50 loci valid; 279 `TESTED`, 71 `NOT_RUN`, 0 `FAILED`; 0 missing/invalid | 828 KB |
| Four workers, 5 loci per chromosome batch | 50 loci / 350 cells, 10 chromosomes, 10 R batches | 992.9 s | 2,615.3 s | 907 MB | 50/50 loci valid; 265 `TESTED`, 85 `NOT_RUN`, 0 `FAILED`; 0 missing/invalid | 1.7 MB |

The 10-locus run took 20.8 minutes and ended with one 10-locus batch after the first four batches finished. The 5-locus run took 16.5 minutes; work was split into ten chromosome batches and the final tail was two five-locus batches. The tested samples are disjoint and cover different chromosomes, so their wall-time difference is descriptive, not a controlled causal estimate of chunk-size speedup. Startup was about 1–3 seconds per R worker; reference load was under 1.8 seconds per batch. Input verification read and hashed 8.9 MB and 9.4 MB, respectively.

All 100 benchmarked loci have valid receipts. Every batch completed on its first attempt. The benchmark shows substantial per-locus cost variation, so the finer five-locus chunks provide more scheduling opportunities at a small measured worker-startup cost. The locked production batch size is therefore **5 loci**, with **4 workers**. The 5-locus observed rate projects to about **13.8 hours** for 17,465 cells; the 10-locus observed rate projects to about **17.3 hours**. Treat this 13.8–17.3 hour range as a planning estimate; the full family runtime may differ with locus composition, concurrent workloads, and checkpoint verification.

At the 5-locus sample’s measured output footprint, linear extrapolation is about **87 MB** for the full family, including per-batch logs and the full-length result table. The estimate is conservative only to the extent that production batch count and log volume scale similarly. `ru_maxrss` is the largest individual child-process peak, not the sum of simultaneous workers.

## Correctness and QC evidence

- 70 Brain6 script tests passed; 308 extension tests passed; 3 native tests were skipped because `susieR`, `TwoSampleMR`, and PLINK 1.9 are unavailable.
- Ruff, Python byte-compilation, R parsing, and `git diff --check` passed.
- Tests cover exact manifest uniqueness, the four-worker lock, deterministic chromosome batching, input/config/output tamper rejection, immutable run separation, resume verification, aggregate denominator, and process-failure retry limits.
- Both 50-locus subsets had `NOT_RUN` fractions above 5% (20.3% and 24.3%). They are diagnostic only; their expected subset QC failure blocks promotion. The full production audit will independently calculate all 17,465 statuses against the frozen 5% threshold.
- The first benchmark worker source snapshot is retained under `/private/tmp/lava-v3-benchmark-v4/7bb80075617d2ea8a3a7d03dddeca0ab81208e241ebceaca51bfc1563c37c9/source/`; the five-locus benchmark coordinator and worker source snapshots are under `/private/tmp/lava-v3-benchmark-v4-b5/06ca9dea959683dbca4f1c267527875f7bf9f51967d443680a8e68f4498471c7/source/`. Their immutable run identities and worker log hashes are recorded in each `latest_audit.json`.

## Production isolation

Production is running under the current five-locus/four-worker locks in the new output root `work/lava-canonical-v3-production`, run ID `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`. At 2026-09-24 01:19 UTC, a full receipt verifier validated all 453 observed locus receipts (3,171/17,465 cells) with zero errors, and the coordinator remained live. The earlier 00:28 UTC checkpoint had 179 verified loci, so the observed increase is 274 loci in 51 minutes (about 5.4 locus checkpoints/minute). If that short interval held for the remaining 2,042 loci, the arithmetic projection would be about 6.3 hours; this is a rough live-run projection, not a controlled benchmark, and locus cost varies. The original v2 baseline, its failed aggregation audit, and the roundoff run remain separate and untouched. No family-level decision or downstream pairwise output is authorized before the full 17,465-cell audit passes.

### Updated live production throughput (2026-09-24 02:26 UTC)

A later receipt audit validated 734/2,495 loci at 02:26 UTC (5,138/17,465 cells; zero invalid), compared with 453 at 01:19 UTC. This observed interval adds 281 loci over 67 minutes, or 4.19 receipt-verified loci/minute. At that rolling rate, the remaining 1,761 loci would take about 7.0 hours. This is a descriptive in-run projection from two verified checkpoints, not a benchmark or a guaranteed finish time; the previous 5.4 loci/minute interval and current rate show meaningful time variation. The production configuration remains four workers and five-locus batches. No scheduler, family threshold, or output path was changed.


### Historical production runtime profile (2026-09-24 03:31 UTC)

A receipt-bound, read-only log profile verified 842/842 current receipts (5,887/17,465 cells; zero invalid). It joined all 842 per-locus timing rows through each verified receipt's worker config to the corresponding five-locus batch log. Mean/median/p90 per-locus compute wall time was 64.90/51.45/122.09 seconds; median result writing was 0.002 seconds and median reference startup 0.751 seconds. Receipt rates were 100/hour over the last 15 minutes, 103/hour over 60 minutes, and 165/hour over 120 minutes. Arithmetic projections for the remaining 1,653 loci are 16.53, 16.05, and 10.02 hours for those respective windows; they are not guaranteed finish times. This rate variation has no established cause, and no CPU-utilization trace was available. The frozen four-worker/five-locus production schedule was left unchanged. Read-only verification at 04:26 UTC confirmed 1,027/1,027 receipts (7,189/17,465 cells; zero invalid), with 5,792 TESTED, 1,397 NOT_RUN (all `LOW_LOCAL_H2_UNDERPOWERED`), and zero FAILED. The 1,397-cell untested lower bound exceeds the immutable cap of 873 by 524, so v3 cannot pass family QC. A prior coordinator interruption from transient workspace disk exhaustion was resumed under the same run ID and output root after free space recovered; valid receipts were reused, no partial cell files remained, and the four-worker scheduler is active again. The historical profile is unchanged. Latest audit command and projections: `qc/lava_canonical_v3_runtime_profile_20260924.md`.
