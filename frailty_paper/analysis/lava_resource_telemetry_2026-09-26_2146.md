# Frailty LAVA live resource telemetry — 2026-09-26 21:46 UTC

Read-only resource sample taken while the frozen FI × 12-sleep LAVA family continued. No process, claim, receipt, or scheduler-control file was changed by this observation.

## Host snapshot

- Host: 8 logical CPUs, 8 GiB physical memory.
- `top -l 1` at 21:46:23 UTC: load average 29.69 / 23.70 / 18.49; CPU 64.36% user, 35.63% system, 0% idle; 7,417 MiB physical memory used, 214 MiB unused, 3,439 MiB compressed.
- Swap: 5,452.12 / 6,144 MiB used (88.7%; 691.88 MiB free).
- `memory_pressure` reported 30% system-wide free memory.
- Two one-second `iostat` samples reported external `disk4` at 37.69 MB/s then 2.65 MB/s. These short samples do not establish sustained I/O contention.

## Runner snapshot

At 21:46:44 UTC, the external runner state reported 18,006 verified receipts of 29,940, two effective workers from six requested, zero launch failures and zero stale-claim recoveries. Workers 72535 and 72534 held separate claims for `sleep_efficiency` locus 1626 and `sleep_timing` locus 1626, respectively. The persisted concurrency control remains at two effective workers after prior swap-triggered downshifts. Two distinct LAVA R child processes were present in the 21:46:23 process listing; their active trait/locus claims had advanced by the subsequent runner-state read.

From runner count 17,842 at 21:20:14 to 18,006 at 21:46:44, the reported count increased by 164 in 26.5 minutes (~373/hour). This provisional rate implies about 32 hours for the 11,934 runner-count remainder; it is not a receipt-integrity scan and is not a completion guarantee. The latest full receipt integrity audit remains the 21:25 checkpoint, which independently validated 17,887 receipts with zero receipt or claim issues.

## Scheduling decision

Concurrency remains at two. The current sample shows no CPU idle capacity, a high host load, and swap use near 89%, following prior downshifts from six to four, then through three to two under the same guard. Increasing concurrency now would add pressure during a busy host period. Existing workers remain alive and progress between loci; the locked configuration, thresholds, input hashes, and receipt semantics are unchanged.
