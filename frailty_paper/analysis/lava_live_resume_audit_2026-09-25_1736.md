# Live FI×sleep LAVA incremental audit — 2026-09-25 17:36 UTC

## Runner and receipts

At 17:36:32 UTC, the runner reported 8,660/29,940 verified receipts (21,280 remaining), two effective workers from a four-worker request, two active claims, zero launch failures and zero stale recoveries. The claims were napping locus 929 (worker PID 22178) and sleepdur locus 511 (worker PID 22177). The coordinator session remained attached, and the receipt-file count matched runner state.

Twenty-seven receipt events between the 17:23:47 and 17:36:32 snapshots passed the serial runner's `verify_receipt` identity/status checks; all corresponding locus logs were nonempty. The append-only event log contained 3,787 receipt keys, with no duplicate `(trait, locus_index)` keys or malformed rows. This is incremental validation, not a full receipt-manifest rehash.

## Resource state and provenance

At 17:36 UTC, `vm_stat` reported 1,413 free pages at 16 KiB/page (about 22 MiB); the external SSD had 1.6 TiB free. Fresh CPU, I/O and swap-total sampling was unavailable under host restrictions. The coordinator's last swap trigger was 5,856/7,168 MiB, above its 3,584 MiB deep-pressure threshold, so two workers remain appropriate. The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; scientific inputs, parameters and QC thresholds were unchanged.

## Throughput and interpretation

The runner advanced 27 receipts from 8,633 at 17:23:47 to 8,660 at 17:36:32, about 127 receipts/hour over this short interval. At that rate, 21,280 remaining receipts would take about seven days; the estimate is provisional. First-pair and chronotype interim gates still exceed the locked failure allowance. No local-sharing inference is supported by this incomplete family.
