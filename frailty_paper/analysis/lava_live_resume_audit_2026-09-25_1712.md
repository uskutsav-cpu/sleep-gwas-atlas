# Live FI×sleep LAVA incremental audit — 2026-09-25 17:12 UTC

## Runner and receipts

At 17:12:24 UTC, the runner reported 8,611/29,940 verified receipts (21,329 remaining), two effective workers from a four-worker request, two active claims, zero launch failures and zero stale recoveries. The claims were longsleep locus 509 (worker PID 22177) and snoring locus 926 (worker PID 22178). The existing coordinator session was still attached at the snapshot. The scheduler remains at two after its deep-swap fallback.

Thirty-nine receipt events from the 16:49:32 runner snapshot window passed the serial runner's `verify_receipt` identity/status checks and each corresponding log was nonempty. The append-only event log contained 3,738 receipt keys at this snapshot; they were unique, with zero malformed rows, duplicate `(trait, locus_index)` keys, or incremental receipt/log errors. This is incremental validation, not a rehash of every receipt file.

## Resource state and provenance

At 17:12 UTC, `vm_stat` reported 3,964 free pages at 16 KiB/page (about 62 MiB); the external SSD had 1.6 TiB free. Fresh CPU, I/O and swap-total sampling was unavailable under host restrictions. The last coordinator-recorded swap sample was 5,856/7,168 MiB, above the 3,584 MiB deep-pressure threshold, so two workers remain appropriate. The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; scientific inputs, parameters and QC thresholds were not changed.

## Throughput and interpretation

The runner advanced 38 receipts from 8,573 at 16:49:32 to 8,611 at 17:12:24, about 100 receipts/hour over this short interval. At that rate, 21,329 remaining receipts would take about nine days; the estimate is provisional. First-pair and chronotype interim gates still exceed the locked failure allowance. No local-sharing inference is supported by this incomplete family.
