# FI×sleep LAVA live resume audit — 2026-09-25 22:32 UTC

## Live execution

The replacement coordinator session (handle 78739) was polled and remained live. At 22:32:23 UTC, runner state and the receipt-file count reconciled at 8,985/29,940. Four workers remain requested; the adaptive guard has reduced effective concurrency to two following the recorded deep-swap trigger (3→2 at 4,156/5,120 MiB at 22:05:47 UTC). Two live claims were recorded: shortsleep locus 528 (worker PID 22177) and insomnia locus 945 (worker PID 22178). Launch failures and stale-claim recoveries are zero.

## Receipt and event audit

All 8,985 receipt files passed the project's `verify_receipt` identity/status validator. The state count matched the file count. The append-only event log contained 4,159 parseable events, including 4,112 receipt events with 4,112 unique (trait,locus) keys and zero malformed rows. State was updated at 22:32:23 UTC.

The frozen analysis-lock SHA-256 matches both runner state and the input manifest: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. The input-manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. No scientific inputs, thresholds, or output semantics changed.

Host swap usage could not be freshly read from this shell; retain two workers under the coordinator's existing adaptive state. The external SSD reports 1.6 TiB free. The first-pair and chronotype interim QC gates remain failed under locked criteria; no local-sharing inference or downstream gated analysis is supported. The systematic review and full paper package remain incomplete.
