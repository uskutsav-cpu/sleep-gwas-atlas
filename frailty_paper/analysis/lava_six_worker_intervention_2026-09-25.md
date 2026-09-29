# LAVA six-worker intervention — 2026-09-25 07:03 UTC

## Decision

The requested six-worker intervention was evaluated against the live run and host state. No coordinator or locus worker was started. The verified pause remains intact because the host is already under severe memory pressure; launching six jobs would violate the resource-safety condition in the intervention request. The previous measured six-worker interval also triggered the four-worker fallback and delivered lower throughput than the recent four-worker interval.

## Run and scheduler state

- Frozen family: 12 traits × 2,495 loci = 29,940 `(trait,locus)` jobs.
- Current receipts: 8,176; remaining jobs: 21,764.
- Runner state: `PAUSED_AFTER_WORKERS_EXITED`, last updated 2026-09-25 06:12:04 UTC.
- `runner.lock`: absent. Non-pause claims: 0. Exact pause holds: 21,764. All recorded worker process groups checked dead; no LAVA/R locus process was visible in the elevated process listing.
- Coordinator and workers running after intervention: 0. No trait/locus assignment changed.
- Existing scheduler already uses exclusive per-job claim creation, validates an existing receipt before scheduling, gives each worker its own temporary directory, and leaves global `runner_state.json` updates to the coordinator. No scheduler or scientific configuration edit was needed.

## Resource and throughput evidence

- Host: 8 logical processors, 8 GiB RAM.
- Sampled CPU: 23.7% idle; load average 4.28.
- Physical memory: 7,496 MiB used, 134 MiB unused; macOS reported 54% system-wide memory free by its pressure estimate.
- Swap: 7,443.56 / 8,192 MiB used (748.44 MiB free).
- Filesystem free space: internal 7.4 GiB; external analysis SSD 1.6 TiB.
- Disk I/O sample completed, but the sandboxed `iostat` could not read the kernel boot time; the elevated snapshot reported high system-disk activity and about 2.05 MiB/s on the external analysis disk during the second sample.
- Earlier six-worker throughput was about 392 receipts/hour and caused severe swapping. The recent four-worker interval completed 130 receipts in 12.8 minutes (~608/hour); that brief rate is not a new completion forecast. No post-intervention throughput or ETA exists because the run was not resumed.
- Current launch-failure count and duplicate-job count at the paused checkpoint: 0 and 0, respectively. Stale claims recovered: 0.

## Frozen inputs and QC

- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` (matches the prior verified value).
- Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` (matches the prior verified value).
- QC thresholds, loci, references, parameters, and receipt semantics were not changed. Chronotype and first-pair QC failures remain failures; no local-sharing claim is supported.
- Prior independent audit verified all 8,176 receipts, their identities/checksums, no duplicate receipt pairs, and exact equivalence between pause holds and unfinished jobs. No receipts or external runner state were modified during this intervention.

## Verification

- LAVA scheduler, pause, and receipt-audit tests: 16 passed.
- Available frailty suite excluding three modules that cannot import missing `requests` and `PyYAML` dependencies: 96 passed, 92 subtests passed.
- Full test discovery was attempted and stopped at collection because `requests` and `yaml` are missing from the active Python 3.13 environment. The repository's `frailty_paper/.venv` has no `pytest` module.
- `git diff --check` passed for this audit and progress entry.

The current shared-host swap pressure makes even the existing four-worker fallback unsuitable to resume now. Reassess after memory and swap recover; then use the existing four-worker safety cap unless a fresh measured six-worker interval can remain within the specified resource limits.
