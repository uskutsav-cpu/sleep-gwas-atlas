# Live FI×sleep LAVA resume audit — 2026-09-25 15:17 UTC

## Runner state

At 2026-09-25T15:17:26.013011Z, the external runner state reports 8,430 of 29,940 expected receipts, two effective workers, two active claims, zero launch failures and zero stale-claim recoveries. runner.lock is present and the state is advancing. The append-only concurrency record preserves the authorized 4→3 transition at 14:47:54 and 3→2 transition at 14:52:53 after swap thresholds were reached. The current pause request is false.

## Independent receipt-event audit

Audited each receipt event after the verified 15:00 checkpoint against its corresponding receipt file and locus log, using the same verify_receipt schema/identity validator as frailty_paper/scripts/45_run_lava_sensitivity.py. All 22 post-checkpoint events have unique (trait,locus) keys; every receipt passes identity/status checks and has a nonempty log. No malformed event lines, duplicate receipt keys or identity/log errors were found. The 15:00 receipt-manifest contains 8,408 receipt hashes; this audit intentionally validates only the 22 incremental receipts while the run continues.

Frozen analysis-lock SHA-256 remains 74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70; prepared-input manifest SHA-256 remains 8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2.

## Host sample and action

The system reports load averages 30.92/29.72/29.88 on eight cores and 39% system-wide free memory. The sysctl vm.swapusage query is denied in this shell, so no current numeric swap value is claimed. Continue with two workers under the existing adaptive fallback; do not increase concurrency based on this sample.
