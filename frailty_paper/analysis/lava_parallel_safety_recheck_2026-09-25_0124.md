# LAVA parallel safety recheck — 2026-09-25 01:24 UTC

## Decision

Keep the live 29,940-slot campaign at **four total workers**. Do not add workers 5–6 under the current host conditions. The scheduler supports six and its isolated six-worker trial is preserved in `lava_parallelization_intervention_2026-09-24.md`; that trial was reduced after memory pressure and throughput worsened. Current resource use already meets the request's stated fallback conditions at four workers.

## Live state

- Coordinator PID: `96341`; workers: `91492`, `91493`, `91494`, `91495` (all confirmed alive).
- At the final state read, the coordinator reported `6,872 / 29,940` verified receipts and four active, distinct claims. A receipt scan that overlapped worker writes validated 6,873 files; the one-receipt difference is an expected race with the coordinator's 10-second refresh, not a missing or duplicate receipt. A preceding snapshot reconciled exactly at 6,862 files/state receipts.
- Full receipt scan: all observed receipt files passed `verify_receipt`; zero duplicate `(trait, locus)` identities.
- Launch failures: `0`; stale claims recovered: `0`.
- Locked analysis SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Prepared input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The exact active claims changed as workers finished and picked up new loci; the coordinator state contains the latest assignments. No receipt, checkpoint, input, reference, threshold or scientific setting was modified by this recheck.

## Resource basis

At 01:20 UTC on the 8-logical-CPU host, system-wide `top` reported 134 MB unused physical RAM, 4,683 MB of 5,632 MB swap used, and 14.3% CPU idle. The SSD had approximately 1.7 TiB free; the short disk sample showed active throughput between about 5 and 11 MB/s. Four concurrent R locus jobs were individually using about 0.2–0.7 GB resident memory. Adding two such jobs at the present memory/swap level would exacerbate swapping. In the earlier measured six-worker trial, throughput was about 392 receipts/hour with 166 MB unused RAM and 4.77 GB swap used; four-worker observations later reached 546/hour and the latest interval measured about 622/hour. The current four-worker rate is therefore the safer and faster observed setting.

From 6,655 receipts at 01:03:09 UTC to 6,872 at 01:24 UTC, observed throughput was approximately 622 receipts/hour, implying about 37 hours for the remaining work if that interval holds. This is a short, host-load-dependent estimate.

## Validation

- The full frailty Python suite passed: `122/122` tests, including the atomic-claim race, deterministic job enumeration, and stale-claim recovery tests.
- `git diff --check` passed.
- No worker was interrupted; the existing four workers and coordinator remain active.
- Six-worker scheduling implementation remains in separate commit `87e6cee77064d1b0637a42be591df906f4360580` (`Parallelize frailty LAVA runner with six safe workers`). No duplicate implementation change was needed.

The current decision is operational rather than scientific. The frozen LAVA configuration, inputs, parameters, locus family and QC thresholds remain unchanged. Chronotype's locked QC failure remains in force and no local-sharing inference is promoted.
