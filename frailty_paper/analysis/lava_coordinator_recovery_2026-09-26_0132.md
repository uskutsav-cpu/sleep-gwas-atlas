# FI×sleep LAVA coordinator recovery — 2026-09-26 01:32 UTC

## Recovery evidence

The prior attached coordinator session handle was no longer available. At 01:27 UTC, external runner state was stale at 9,055/29,940 (last update 2026-09-25 22:52:36 UTC), while 9,365 receipt files and 4,491 unique receipt events existed. The latest event and receipt were timestamped 23:11:18 and 23:11:22 UTC. The stale runner lock named PID 52234; its process and registered worker groups 22177, 22178 and 52535 were dead. Both ordinary active claims belonged to dead process groups. Two additional ._ claim files were unreadable AppleDouble sidecars and were ignored by the scheduler; they were not treated as jobs.

The old worker logs end with PermissionError while appending the event log or removing a claim on the external run volume. Scoped write permission to the existing run root was restored. No old worker or locus process was terminated.

## Resume and verification

The adaptive coordinator was relaunched at 2026-09-26 01:31:13 UTC with four requested workers. Its startup validated the frozen inputs and all 9,365 existing receipts, then recovered exactly two stale normal claims. It launched four distinct workers: 74231–74234. Analysis-lock SHA-256 remains 74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70; input-manifest SHA-256 remains 8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2.

At 01:32:43 UTC, runner state and receipt files reconciled at 9,381/29,940. Four active atomic claims belonged to distinct (trait,locus) keys. Launch failures were zero and stale claims recovered remained two. A fresh resource check measured swap at 212.75/1,024 MiB and system-wide memory free at 47%; no downshift was triggered. External SSD free space was 1.6 TiB.

This recovery preserves all receipts and the frozen scientific configuration. The first-pair and chronotype interim gates still fail under locked thresholds; no local-sharing inference is supported. The overall paper remains NO-GO pending completed analyses, licensed search exports and independent dual screening.
