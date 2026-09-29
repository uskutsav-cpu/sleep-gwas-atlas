# FI×sleep LAVA live resume audit — 2026-09-25 22:27 UTC

## Live execution

The replacement coordinator PID 52234 remains active in session 78739. The run requested four workers; resource safeguards have left two effective workers (PIDs 22177 and 22178). At 22:26:57 UTC, runner state and receipt-path count reconciled at 8,970/29,940, leaving 20,970 jobs. Active atomic claims were napping locus 944 and sleep_timing locus 527. Launch failures and stale-claim recoveries were zero.

## Incremental receipt and event audit

From the reconciled checkpoint at 22:11:12 UTC (8,881 receipts) to 22:26:57 UTC (8,970), 89 receipts were added over 15.76 minutes, approximately 338 receipts/hour. The estimate for the remaining family is about 62 hours if that short-window rate persists; it is provisional. Eighty-three receipts created during the first incremental audit passed `45_run_lava_sensitivity.py::verify_receipt` (identity/status and nonempty-log checks). The remaining outputs in this interval were produced after that audit and passed the worker's same validator before their receipt events were appended. The complete event log then contained 4,144 rows, zero malformed rows, and zero duplicate receipt keys.

The frozen analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Scientific settings and QC thresholds were unchanged. At the 22:27 resource sample, swap was 6,774/9,216 MiB and system-wide free memory was 30%. Retain two workers; the absolute deep-swap trigger remains exceeded.

The first-pair and chronotype interim completeness/QC gates remain failed under the locked threshold. No family-level local-sharing inference or downstream PLACO, fine-mapping, colocalization, or molecular conclusion is supported. The systematic-review queue remains unscreened pending licensed exports and human dual screening. Overall readiness remains NO-GO.
