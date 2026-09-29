# FI×sleep LAVA six-worker intervention — 2026-09-26 02:27 UTC

## Safe resume and measured outcome

The four-worker coordinator was paused at a locus boundary after its six active claims had finished. At 02:12 UTC, the paused state had no runner lock, no ordinary claims, and 9,486 receipt files. An independent pass through the runner's receipt validator confirmed all 9,486 identities and schemas. The receipt count and files agreed; no scientific inputs were changed. The analysis-lock SHA-256 remained `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`, and the prepared-input manifest SHA-256 remained `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The existing coordinator was resumed once with `--workers 6 --resume-paused` (six total workers, not six additional workers). It verified the prepared inputs and launched worker PIDs 78548–78553. The first live snapshot showed six distinct atomic claims across insomnia, longsleep, napping, shortsleep, and sleep timing. No locus was killed, no completed receipt was recomputed, and no launch failed.

From 9,486 receipts at 02:15:21 to 9,516 at 02:19:29, the six-worker interval produced 30 receipts in 4.13 minutes, about 436 receipts/hour. This short rate is provisional and is not a completion estimate. The scheduler then applied its existing swap guard at safe locus boundaries:

- 6 → 4 at 02:21:02, when swap reached 1,593/2,048 MiB and crossed the 1,536 MiB absolute trigger.
- 4 → 3 at 02:21:12, when swap reached 1,910/3,072 MiB.
- 3 → 2 at 02:22:26, when swap reached 2,772/3,072 MiB (90%).

At 02:27:19 UTC, the coordinator was still live with two effective workers (PIDs 78548 and 78549), 9,536/29,940 receipts, and two distinct claims: sleep-apnea locus 972 and snoring locus 546. The requested count remains recorded as six; workers 3–6 have exited. Launch failures and stale recoveries are zero. The campaign is 31.9% complete, with 20,404 receipts remaining. No stable two-worker throughput or ETA is available yet.

The latest host sample measured swap at 3,835/5,120 MiB, 34% system-wide free memory, 144 MiB physically unused, load averages 26.17/31.66/24.22, and 2.12% CPU idle. External SSD free space remained 1.6 TiB. The system remains under heavy host-wide pressure; the coordinator is at its two-worker floor. Historical fallback controls and scientific settings were preserved.

## Receipt, claim, and configuration audit

An independent full receipt validation after fallback passed all 9,535 receipt files with zero identity/schema errors. The runner's state and receipt-file count reconciled at the audit snapshot; each newly published receipt is also validated by its worker before its event is written. Across 4,662 receipt events, there were zero duplicate keys or malformed rows. The two ordinary claims had distinct `(trait, locus)` identities. There were zero launch failures, zero stale claims, and no evidence of concurrent duplicate jobs.

The analysis-lock SHA-256 and prepared-input manifest SHA-256 still match the frozen values above. The coordinator remains the only writer of `runner_state.json`; workers use atomic per-locus claims, immutable receipts, and the append-only event log. No loci, references, parameters, QC thresholds, or result semantics changed.

## Package validation during the run

- The full `make -C frailty_paper test` discovery ran 152 tests but reported one import error: the untracked Brain6-specific `test_audit_brain6_sleep_power_pilot` requires `pytest`, which is absent from the isolated Python 3.11.11 environment. No Brain6 files or dependencies were changed.
- The frailty-only suite, excluding that Brain6-only import failure, passed 151/151 tests under Python 3.11.11. The captured log is `analysis/test_suite_frailty_only_2026-09-26_0219.log` (SHA-256 `c1af2809c5a046be6f0c91caa4196ce460f5b4286baa18cf16232ff675f697f4`).
- `git diff --check` passed before these documentation updates.
- `make -C frailty_paper validate-plan` passed: locked v1, 12 sleep traits, 396-pair multiplicity family.
- The frozen FI audit passed (12 rows; nine all-396 FDR-significant frozen global-rg rows) and correctly reported that source reprocessing is incomplete because required h2/harmonized/munged inputs are unavailable.
- The broader `smoke` target stopped safely at the existing frozen-extension audit. It refused to overwrite the tracked audit because the current checkout has nine missing pinned core artifacts, including `results/tables/analysis_panel_provenance.tsv`, and one mismatch in the append-only `environment/tool_versions.tsv` rows. The diagnostic was written only under `/private/tmp`; no frozen extension or Brain6 artifact was modified.

## Remaining project gates

The interim first-pair and chronotype LAVA gates still fail the locked QC threshold, so there is no local-sharing inference. The complete 29,940-receipt campaign and full collator remain pending. The 56,117-record review queue has no human screening decisions; licensed Embase, Scopus, Web of Science, and PsycINFO exports and independent dual screening remain required. Exact cohort-overlap and source-eligibility questions also remain open. Overall readiness stays **NO-GO**.
