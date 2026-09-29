# LAVA six-worker safety recheck — 2026-09-25 07:34 UTC

## Decision

The FI × sleep LAVA campaign remains paused. No process, claim, receipt, or configuration was changed and no worker was launched. The existing scheduler implementation already provides exclusive per-locus claims, receipt validation/skipping, per-worker scratch directories, and coordinator-only global state writes; its code and tests were previously committed in `87e6cee` (`Parallelize frailty LAVA runner with six safe workers`). The later explicit-pause/recovery safeguards are in commits `019fb9b` and `f8e420d`.

Resumption is unsafe in this snapshot. Four unrelated Brain6 R workers are actively using the shared host, swap is nearly full, physical RAM is nearly exhausted, and the system disk reported very high I/O. The user requested that unrelated Brain6 work remain untouched, so these processes were left running. No healthy process was terminated.

## Live campaign state

Read-only inspection of `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1` found:

- `runner_state.json`: `PAUSED_AFTER_WORKERS_EXITED`; 8,176 of 29,940 receipts; zero workers; last state timestamp 06:12:04 UTC.
- `concurrency.json`: `pause_requested=true`, pause ID `pause_20260925T060720Z_b646`, desired workers 4, pause reason is Brain6/FI scope separation.
- `runner.lock`: absent.
- Receipt files: 8,176.
- Claim files: 21,764; these are the verified pause holds from the prior integrity audit, not active jobs.
- No LAVA/R locus worker appeared in the live process listing. The four R workers visible were running `brain6/scripts/run_power_optimized_sleep_univariate_v1.R`.

The previously preserved independent receipt verification at `analysis/lava_pause_receipt_reverification_2026-09-25_0622.md` records 8,176 valid unique receipts, exactly 21,764 pause holds for unfinished jobs, no malformed/duplicate receipts, and unchanged frozen input hashes. This recheck did not rewrite any external state.

## Host snapshot

At approximately 07:34 UTC, the host reported 8 logical processors, load averages 13.57 / 10.12 / 7.92, and 73.4% CPU idle in the instantaneous sample. Physical memory was 7,414 MiB used, including 3,278 MiB compressed, with 215 MiB unused. Swap was 8,625 / 9,216 MiB used (591 MiB free). The internal volume had 5.9 GiB free; the external analysis SSD had 1.6 TiB free. Disk sampling showed high system-disk activity, including one 1-second sample at approximately 694 MiB/s reported aggregate throughput. These conditions do not support adding six memory-intensive LAVA workers alongside four active Brain6 R workers.

The prior measured six-worker trial was approximately 392 receipts/hour and triggered fallback; the latest four-worker LAVA interval was approximately 608/hour over 12.8 minutes, but was stopped by the scope pause. There is no current parallel throughput or completion estimate because LAVA is paused. Launch failures and duplicate jobs remain zero in the verified checkpoint; stale claims recovered remain zero.

## Scientific and recovery status

No loci, reference data, parameters, input hashes, QC rules, receipt semantics, or output settings were changed. Chronotype and first-pair QC failures remain failures; no local-sharing inference is supported. The pause holds remain in place. Reassess after the Brain6 workers finish and memory, swap, and I/O recover; resume only after the existing pause control is explicitly cleared, using no more than four workers unless a new measured six-worker interval satisfies the resource and throughput rules.

## Verification performed

- Read-only process, runner-state, pause-control, receipt-count, claim-count, memory/swap, disk-space, and disk-I/O checks.
- No files on the external analysis volume were modified.
- This report records the recheck only; scheduling code is unchanged.
