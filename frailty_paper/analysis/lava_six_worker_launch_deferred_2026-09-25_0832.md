# Six-worker LAVA launch deferred — 2026-09-25 08:32 UTC

## Decision

The requested six-worker resume was not launched. The existing campaign remains paused under `pause_20260925T060720Z_b646`. The coordinator and worker scripts already support atomic per-locus claims and the six-worker mode; the last six-worker trial had already been reduced after severe swapping and lower measured throughput. A prior measured six-worker rate was about 392 receipts/hour; a later four-worker interval reached about 608/hour over a short sample. Neither rate justifies adding load to the current host.

## Live state and host sample

- `runner_state.json`: `PAUSED_AFTER_WORKERS_EXITED`, 8,176 / 29,940 receipts, 21,764 pause holds, zero LAVA workers, zero launch failures.
- No LAVA coordinator or LAVA R worker was present in process inspection; no `runner.lock` exists. The campaign control still requests four workers and remains paused.
- The most recent independent receipt/claim audit (06:23:41 UTC) validated the 8,176 saved receipts, exact coverage of all 21,764 unfinished jobs by pause holds, no duplicate claims/events, and unchanged locked input hashes.
- Host sample at 08:32: load average 9.87 on eight logical CPUs; 156 MiB physically unused; 3,397 MiB occupied by the memory compressor; swap 8,199.88 / 10,240 MiB used. Four unrelated Brain6 R workers and a separate SMFVI Python analysis were active. Disk sampling showed high system-wide traffic; the external SSD had about 1.6 TiB free.

These conditions are insufficient for a six-worker restart. Pause claims were not released. No receipts, claims, controls, logs or scientific inputs were changed by this check.

## Validation

- LAVA scheduler and pause unit tests: 7 / 7 passed under the pinned Python 3.11 environment.
- Frailty unittest discovery: 132 tests ran; 131 passed, with one import error from the untracked Brain6-only `test_audit_brain6_sleep_power_pilot.py`, which requires `pytest` absent from the pinned environment.
- The initial system-Python pytest attempt also failed collection because system Python 3.13 lacks `lxml`; the pinned Python environment does not include pytest.
- No full preflight or receipt recomputation was run. The prior independent receipt audit and external checkpoint were preserved.

## Resume gate

Re-sample host CPU, memory pressure, swap, competing processes and I/O before resuming. Start only at a completed-locus boundary; preserve all receipts and frozen analysis/input hashes. If six workers again cause material swapping, I/O contention, worker failures or lower throughput, reduce to four. If four workers are not safe either, remain paused.
