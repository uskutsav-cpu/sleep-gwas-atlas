# Phase 1 integrity and live LAVA checkpoint — 2026-09-25
Compiled at 2026-09-25 00:47 UTC in the `frailty-paper-v1` worktree. This is an incremental checkpoint, not a full integrated preflight or a scientific eligibility decision.

## Acquired-resource integrity

- External-root manifest verification passed: 351/351 registered files matched recorded path, size and SHA-256.
- The 22-column metadata audit passed for all 351 rows.
- The manifest SHA-256 is `61cfa0052eb64009d1e34090e4f54d116457e1137af5e2203de2a9e6d0d04df1`. Captured command output and UTC interval are in `analysis/acquisition_manifest_reverification_2026-09-25_0031.log` (SHA-256 `f6f23ac9227180a962178569fc779183963a6bd7c1ef99ce1dedc9e6e93b0028`).
- The external workspace storage gate passed with 1,698 GiB free at 00:32 UTC. At 00:38:14 UTC, the internal filesystem had 2.1 GiB free and the external volume had 1.7 TiB free; keep large outputs and temporary data on the external workspace.

## Live sensitivity run

- At 2026-09-25T00:47:53.933828+00:00, the coordinator reported 6,493/29,940 valid receipts (23,447 remaining), four workers, zero launch failures, zero stale-claim recoveries, and lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Coordinator PID 96341 and workers 91492–91495 were confirmed live at 00:46 UTC. Active claims at the coordinator snapshot: longsleep locus 871 (PID 91493), napping locus 453 (PID 91492), shortsleep locus 1708 (PID 91495), sleep_apnea locus 1289 (PID 91494). The four `(trait, locus)` identities are distinct.
- From 6,143 receipts at 00:06:52 UTC to 6,493 at 2026-09-25T00:47:53.933828+00:00, the observed rate was about 512 receipts/hour. At that measured rate, the remaining work is about 45.8 hours; this is provisional.
- The current chronotype audit at 00:30:11 UTC found a contiguous prefix of 2,376/2,495, 1,326 flagged loci (53.15%), and 21 unexpected reference-SNP failures. The 24-locus frozen limit is exceeded. The diagnostic is `analysis/lava_interim_receipt_diagnostic_chronotype_2026-09-25_0027.json` (SHA-256 `d3497998eb108c484f56d530e76115983a74edd5390aaf2630ed610b8b332548`). It matches the previous prefix/count, so Figures 4–5 were not regenerated.
- The frozen configuration and thresholds are unchanged. No local-sharing inference or PLACO result is supported by the failed gate.

## Validation

- The pinned Python 3.11.11 frailty suite passed 119/119 tests. Captured output: `analysis/post_manifest_reverification_unittest_2026-09-25_0033.log` (SHA-256 `0e34fab36d312bd46e6fff0b019cf947a040ac9496394ee86289010b5734a9f2`).
- This snapshot does not resolve manual database exports/screening, Fried/HFRS full-statistics access, exact source cohort overlap, or physical-component lineage.
