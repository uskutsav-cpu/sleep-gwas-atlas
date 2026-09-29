# Independent LAVA paused-checkpoint receipt audit — 2026-09-25 06:22 UTC

## Scope and pause state

This is a read-only integrity re-audit of the FI × sleep LAVA sensitivity family. The active external pause control is `pause_20260925T060720Z_b646` (“Brain6/FI scope separation: stop new FI×sleep loci after existing claimed loci finish; preserve all receipts and checkpoints.”). No analysis job was started or changed. The coordinator and workers are stopped, `runner.lock` is absent, and `runner_state.json` reports `PAUSED_AFTER_WORKERS_EXITED` with zero workers.

## Receipt verification

The frozen family contains 12 traits × 2,495 loci = 29,940 expected jobs. The saved receipt manifest `receipt_manifest_20260925T0610Z.tsv` has 8,176 unique receipt paths. Each entry was checked against the SHA-256 of its current file, and each JSON receipt was passed through `frailty_paper/scripts/45_run_lava_sensitivity.py::verify_receipt` using the frozen trait order and locus identity.

- Receipt paths and manifest rows: 8,176 / 8,176
- Receipt checksum errors: 0
- Receipt identity/schema errors: 0
- Duplicate receipt pairs: 0
- Unfinished expected jobs: 21,764
- Pause holds: 21,764 unique holds, exactly equal to the unfinished `(trait,locus)` set
- Non-pause claims: 0; duplicate claim pairs: 0
- Duplicate receipt events: 0; malformed event rows: 0
- Worker failures and stale-claim recoveries: 0

The receipt-manifest file SHA-256 is `c8fb91f748a4c4b56edec61465e7aa5edcbfaedeae0a4c669956e2ebeb8da7db`. The previously saved 06:10 checkpoint remains unchanged at SHA-256 `664d79fe29d78506327172c8e0ac9fc1324be0dfd00f6c4f905aa41ed9624435`; its `analysis_lock_sha256` field was null. An additive checkpoint variant, `paused_checkpoint_20260925T0622Z_independent_reverification.json`, preserves the source checkpoint and adds the verified analysis-lock hash. Its SHA-256 is `1e3f299177829758e63107a351f546c540fffb2443f2a1756c672a2560737d4e`.

## Frozen inputs and interpretation

The analysis-lock file `frailty_paper/config/lava_frailty_sensitivity_v1.yaml` has SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. The prepared-input manifest has SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both match the recorded campaign values. Scientific settings and the frozen 1% / 24-locus completeness rule were unchanged. This verifies receipt identity and preservation only; it does not turn failed first-pair/chronotype QC gates into passing results, and it supports no local-sharing inference or PLACO analysis.
