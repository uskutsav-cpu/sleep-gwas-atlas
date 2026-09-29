# LAVA resource safety recheck — 2026-09-25 08:07 UTC

## Decision

The LAVA sensitivity campaign remains paused. No coordinator or LAVA locus worker was running at the recheck, `runner.lock` was absent, and no pause holds or receipts were changed. The user requested six total LAVA workers, with a reduction to four if resource pressure became substantial. The previous six-worker trial had already produced severe memory pressure and only about 392 receipts/hour. Four workers are the supported fallback, but current host conditions do not safely support resuming even that fallback.

## Current state and preservation

The external checkpoint `paused_checkpoint_20260925T0622Z_independent_reverification.json` reports 8,176 of 29,940 verified receipts and 21,764 remaining jobs held by pause claims. It verifies zero receipt hash or identity errors, zero duplicate receipt events or claims, zero active workers/locus processes, no runner lock, and an exact match between pause holds and expected unfinished jobs. The frozen analysis lock remains SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; the prepared-input manifest remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The checkpoint's last independent receipt verification was at 06:23:41 UTC. Its pause ID is `pause_20260925T060720Z_b646`. No runner state, receipt, claim, lock, or scientific input was modified during this recheck.

## Host resource snapshot

Read-only sampling at approximately 08:05 UTC found an eight-core host with four unrelated Brain6 R workers and an unrelated SMFVI Python analysis active, but no LAVA process. Memory pressure reported 30% system-wide free; the top snapshot showed 99 MiB physically unused and 3,432 MiB compressed. Swap was 9,227.94 MiB used of 10,240 MiB (about 90%). CPU was 21.32% idle with load average 8.30. The short disk sample showed approximately 251–524 MiB/s on the internal disk during subsequent intervals; the external LAVA SSD had 1.6 TiB free and its observed I/O subsided to zero after an initial 6.25 MiB/s sample. Internal filesystem free space was 5.9 GiB.

Given near-full swap, very low physically unused RAM, heavy internal disk traffic, and other active analyses, launching four additional LAVA workers would risk materially worsening resource pressure. Six workers remain ruled out by the prior measured six-worker period and the user's automatic-reduction criteria. No workers were launched. Do not release the pause holds until a fresh resource check shows adequate memory and I/O headroom; at that point resume with four workers, then reassess measured throughput and pressure before considering any increase.

## Provenance

The independent checkpoint records zero receipt validation errors and unchanged analysis/input hashes. The current scheduler and its six-worker atomic-claim implementation remain in the prior commits `87e6cee77064d1b0637a42be591df906f4360580` and `019fb9bcdb7b0c9bbb804f2eed09f80f1d60fb0f`. No scheduler code or scientific configuration changed in this recheck.
