# Live FI×sleep LAVA incremental audit — 2026-09-25 17:23 UTC

## Runner and receipts

At 17:23:47 UTC, the runner reported 8,633/29,940 verified receipts (21,307 remaining), two effective workers from a four-worker request, two active claims, zero launch failures and zero stale recoveries. The claims were shortsleep locus 510 (worker PID 22177) and snoring locus 927 (worker PID 22178). The coordinator session remained attached. A filesystem count immediately after the runner update also found 8,633 receipt files.

Twenty-two receipt events between the 17:12:24 and 17:23:47 snapshots passed the serial runner's `verify_receipt` identity/status checks, and all corresponding locus logs were nonempty. The append-only event log contained 3,761 receipt keys; they were unique, with no malformed rows, duplicate `(trait, locus_index)` keys, or incremental receipt/log errors. This is incremental validation, not a full receipt-manifest rehash.

## Resource state and provenance

At 17:22 UTC, `vm_stat` reported 3,862 free pages at 16 KiB/page (about 60 MiB); the external SSD had 1.6 TiB free. Fresh CPU, I/O and swap-total sampling was unavailable under host restrictions. The coordinator's last swap trigger was 5,856/7,168 MiB, above its 3,584 MiB deep-pressure threshold, so two workers remain appropriate. The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; scientific inputs, parameters and QC thresholds were unchanged.

## Throughput and interpretation

The runner advanced 22 receipts from 8,611 at 17:12:24 to 8,633 at 17:23:47, about 116 receipts/hour over this short interval. At that rate, 21,307 remaining receipts would take about 7.6 days; the estimate is provisional. First-pair and chronotype interim gates still exceed the locked failure allowance. No local-sharing inference is supported by this incomplete family.
