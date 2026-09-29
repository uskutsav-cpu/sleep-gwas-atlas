# FI×sleep LAVA intervention recheck — 2026-09-25 10:04 UTC

## Outcome

The coordinator accepted the request for six total workers, then its resource guard selected four at 09:51:56 UTC when swap reached 9,248/11,264 MiB (82%). Four workers launched at 09:52:33 UTC; no additional workers were started. The scheduler implementation already in commit `87e6cee77064d1b0637a42be591df906f4360580` uses deterministic enumeration, exclusive per-locus claims, valid-receipt skipping, per-worker logs/scratch, and coordinator-only global state writes.

At 09:57:38 UTC the shared control recorded a separate FI/Brain6 scope pause. The workers finished their in-flight loci and all 29,940 expected slots reconciled at 10:01 UTC to 8,206 receipts plus 21,734 `pause_hold` claims. The event log contained 3,333 receipt events with zero duplicate `(trait, locus)` events. Runner state recorded zero launch failures and zero stale-claim recoveries. All worker processes (PIDs 62401–62404) were confirmed defunct; no R locus process remained. The idle coordinator (PID 62290) still held `runner.lock` because its process-group probe treated unreaped zombie children as live. After verifying the exact receipt/hold set and dead workers, SIGINT was sent to that coordinator alone; its `finally` cleanup removed `runner.lock`. Receipts and pause holds were preserved. The pause request remains set as `pause_20260925T095738Z_d138`; no worker is active.

## Resource and throughput observations

The 8-core host was already saturated during the four-worker interval. At 09:56 UTC it reported 0.33% CPU idle, 84 MiB physically unused, and 11,427.75/12,288 MiB swap used (93%). Internal-disk samples later ranged about 445–768 MB/s aggregate; the external analysis SSD had 1.6 TiB free and samples up to 13.6 MB/s. A later sample after workers stopped reported 145 MiB unused and 10,766.62/11,264 MiB swap used (96%). These measurements rule out safely launching six workers and do not support resuming four until the host recovers.

Observed four-worker assignments at 09:56 UTC were worker/PID 1/62401 sleep_timing locus 493; 2/62402 sleepdur locus 911; 3/62403 sleep_efficiency locus 1331; 4/62404 longsleep locus 1754. The receipt count advanced 13 between 09:56:33 (8,193) and 09:58:25 (8,206), a short-interval rate of about 422 receipts/hour. A previous four-worker interval measured about 608/hour; the six-worker trial measured about 392/hour under severe pressure; the prior serial estimate was about 81/hour. At the latest short rate, 21,734 remaining slots imply roughly 52 hours of compute, but this is only a provisional rate-based estimate: the run is paused, and the observation interval is brief.

## Scientific integrity and validation

The frozen analysis lock remains SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; the prepared-input manifest file remains SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. No loci, reference data, parameters, input hashes, receipt semantics, or QC thresholds changed. Chronotype's locked QC failure remains unchanged and no local-sharing inference is promoted.

This recheck found a coordinator shutdown issue: `subprocess.Popen` children were not polled after launch, leaving zombie process groups visible to `killpg(..., 0)` and preventing normal lock cleanup. The coordinator now retains its launched `Popen` objects and polls them each loop so exited workers are reaped. A focused regression test covers that behavior. The frailty test suite passed 146 tests when excluding one untracked Brain6-only test module whose environment lacks `pytest`. An initial broad discovery under system Python 3.13 also failed to import a frailty test because `lxml` is absent there; the pinned Python 3.11 run completed the frailty tests successfully. The full `git diff --check` result is recorded below.

## Resume condition

Keep all 21,734 pause holds and the FI/Brain6 scope pause. Reassess only after the unrelated work and host memory/swap/disk pressure recover. Resume through the existing coordinator with the six-worker request; its recorded resource guard must continue to cap the run at four when swap is at or above 82%. Do not remove or rewrite any receipt or pause hold manually.
