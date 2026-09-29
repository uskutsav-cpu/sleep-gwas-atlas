# Live FI×sleep LAVA incremental audit — 2026-09-25 16:50 UTC

## Runner and receipts

At 16:49:32 UTC, `runner_state.json` reported 8,573/29,940 verified receipts, two effective workers, two active claims, zero launch failures and zero stale recoveries. The invocation requested four workers; the adaptive scheduler remains at two after the deep-swap fallback. At 16:50:39 UTC, process inspection confirmed coordinator PID 21938 and worker PIDs 22177 and 22178 alive; `runner.lock` was present.

Fourteen receipt events since the 16:42:35 checkpoint passed `verify_receipt` identity/status checks, and all corresponding locus logs were nonempty. A scan of the append-only event log found zero malformed rows and zero duplicate `(trait, locus_index)` keys. This is an incremental check, not a full receipt-manifest rehash.

## Resource state and provenance

Around 16:50 UTC, swap was 5,207.81/8,192 MiB; system-wide free memory was 33%; load averages were 20.25/19.64/20.89 on eight cores. Internal free space was 7.9 GiB and external SSD free space was 1.6 TiB. Although the swap fraction was 63.6%, use remained above the coordinator's absolute deep-pressure trigger; retain two workers.

The analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. No scientific input, parameter or QC threshold changed. The first-pair and chronotype interim QC gates continue to exceed the locked failure allowance, so no local-sharing inference is supported.

## Manuscript consistency

The manuscript and Reviewer 2 audit were updated to reflect that the incomplete family is actively running, rather than paused at the earlier 8,176-receipt checkpoint. Post-edit audits passed 28/28 selected quantitative claims and 92 reporting rows/93 locators with zero errors. These scoped audits do not establish full scientific validity or project readiness.
