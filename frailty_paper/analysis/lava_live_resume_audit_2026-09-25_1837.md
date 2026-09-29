# Live FI×sleep LAVA checkpoint and full receipt audit — 2026-09-25 18:37 UTC

## Runner state

At 18:37:21 UTC, the live coordinator reported 8,776 of 29,940 expected receipts (21,164 remaining), four workers requested and two effective under the swap-pressure fallback. Two distinct jobs were active: sleepdur locus 517 (worker PID 22177) and sleepiness locus 934 (worker PID 22178). The coordinator session remained attached. Launch failures and stale-claim recoveries were zero.

## Receipt and event verification

A read-only audit using `frailty_paper/scripts/45_run_lava_sensitivity.py::verify_receipt` validated every receipt present during the audit: 8,776 identities and status fields passed, and every corresponding locus log was nonempty. The receipt-file count was 8,776 at the audit end and matched runner state. The full append-only event log had zero malformed rows and zero duplicate `(trait, locus_index)` receipt keys. During the audit, runner state advanced from 8,775 to 8,776; the audited receipt files and final runner count reconciled. This validates receipt structure and identity, not a fresh checksum rehash against a saved full receipt manifest.

## Resources and frozen configuration

The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Scientific inputs, parameters and QC thresholds were not changed. At 18:37 UTC, macOS `memory_pressure` reported 33% system-wide free memory and the external SSD had 1.6 TiB free. Fresh CPU, disk-I/O and swap-total readings were unavailable in this shell; the coordinator's last recorded trigger was 5,856/7,168 MiB swap used, above the 3,584 MiB deep-pressure threshold. Keep two effective workers under the existing adaptive guard.

## Throughput and scientific interpretation

The run advanced from 8,660 receipts at 17:36:32 to 8,776 at 18:37:21, about 114 receipts/hour. At that short-interval rate, the remaining 21,164 receipts would take about 7.7 days; this is provisional. First-pair and chronotype interim completeness gates still exceed the locked failure allowance. No local-sharing inference is supported. The systematic-review queue, source-access items and exact cohort-overlap questions remain open; project readiness is **NO-GO**.
