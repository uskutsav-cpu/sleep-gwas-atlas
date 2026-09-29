# Live FI×sleep LAVA incremental audit — 2026-09-25 17:08 UTC

## Runner and receipts

At 17:08:16 UTC, `runner_state.json` reported 8,603/29,940 verified receipts (21,337 remaining), two effective workers from a four-worker request, two active claims, zero launch failures and zero stale-claim recoveries. The live claims were sleep-efficiency locus 926 (worker PID 22178) and sleep-duration locus 508 (worker PID 22177). The coordinator session remained attached and active. The recorded resource fallback remains the deep-swap reduction from three workers to two.

Thirty receipt events after the prior 16:50 UTC checkpoint passed the serial runner's `verify_receipt` identity/status checks, and each corresponding locus log was nonempty. Across the append-only event log, all 3,730 receipt keys were unique; there were zero malformed rows, duplicate `(trait, locus_index)` keys, or incremental receipt/log errors. This incremental audit does not rehash all 8,603 receipt files.

## Resource state and provenance

At the 17:08 sample, `vm_stat` reported 3,798 free pages at 16 KiB/page (about 59 MiB). The external analysis SSD had 1.6 TiB free. Process inspection and `iostat`/`sysctl` sampling remain restricted by the host, so this note does not claim fresh CPU utilization, I/O rate, or swap totals. The last coordinator-recorded deep-pressure sample was 5,856/7,168 MiB swap, above its 3,584 MiB deep fallback threshold; keeping two workers is consistent with the configured 4→3→2 policy.

The live runner's analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. No scientific configuration or inputs were modified. Interim failures of the locked first-pair and chronotype QC gates remain in force; this incomplete sensitivity run supports no local-sharing inference.

## Throughput context

The runner advanced from 8,573 receipts at 16:49:32 UTC to 8,603 at 17:08:16 UTC: 30 receipts over about 18.7 minutes, or roughly 96 receipts/hour in this short interval. At that observed rate, the remaining 21,337 receipts would take about 9.3 days; this is a rough, variable-throughput estimate, not a completion commitment.
