# FI×sleep LAVA live resume audit — 2026-09-25 22:35 UTC

## Live execution

Coordinator session handle 78739 was polled at 22:35 UTC and remained live. Runner state at 22:35:39 UTC reconciled with the receipt-file count at 8,993/29,940. Four workers are requested; two are effective under the adaptive guard after the deep-swap fallback (3→2 at 4,156/5,120 MiB at 22:05:47 UTC). Active claims were sleep_apnea locus 945 (worker PID 22178) and sleepdur locus 528 (worker PID 22177). There are zero launch failures and zero stale-claim recoveries.

## Receipt and event audit

All 8,993 receipt files passed the project's verify_receipt identity/status validator. The state and file counts matched. The append-only event log contained 4,168 parseable rows, including 4,121 receipt events with 4,121 unique (trait,locus) keys and zero malformed rows.

The frozen analysis-lock hash matches both runner state and input manifest: 74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70. The input-manifest SHA-256 remains 8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2. No scientific input, threshold, or output semantics changed. A fresh swap query was unavailable; retain the coordinator's current two-worker setting.

The first-pair and chronotype interim QC gates remain failed under locked criteria. No local-sharing inference or downstream gated analysis is supported. Licensed review exports and independent dual screening remain outstanding; overall readiness remains NO-GO.
