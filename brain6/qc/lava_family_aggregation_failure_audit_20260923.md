# LAVA v2 family aggregation failure audit

Audit time: 2026-09-23 20:58 UTC  
Baseline run: `13b80e6b64a5179fec70ba510a89c51241c9d58ec05833cc52161b435aee44fa`  
Separate roundoff run: `af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf`

## Finding

The baseline coordinator completed all 2,495 locus receipts, but exited with code 1 during the family aggregation stage. The error at `brain6/scripts/run_lava_family.py:422` is:

> Pair-specific inputs disagree on the shared univariate test: insomnia/2

The receipt audit validates all 2,495 baseline locus receipts with zero receipt errors. This is not evidence of a completed or accepted LAVA family: the runner did not create `family_decision.json`, and it stopped before writing the aggregate family artifacts.

## Why aggregation stopped

The aggregator keys univariate results by `(trait, locus)` and requires every occurrence to have an identical p-value. The worker runs `LAVA::process.input` and `LAVA::run.univ` separately for each selected trait pair. Consequently, a repeated trait's local test can use a different pair-specific shared-variant set and produce a different p-value. For insomnia at locus 2, the baseline receipts record `p=0.0742928` in `insomnia__adhd` and `p=0.00398049` in `insomnia__mdd`.

A read-only scan of all receipt-verified baseline worker outputs found:

- 18,407 univariate result rows across pairs.
- 14,632 distinct observed trait-by-locus keys, out of the frozen 17,465-test family.
- 2,804 repeated trait-by-locus keys; 2,802 have p-values differing by more than `1e-12` across pair contexts.
- 2,833 unobserved unique trait-by-locus tests: `2,833 / 17,465 = 16.22%`, above the frozen 5% maximum untested fraction.

Thus, changing the aggregator to accept a duplicate p-value would not make this family pass QC: the measured untested fraction independently exceeds its frozen cap. The observed bivariate slot counts are 12,356 `UNIVARIATE_UNDERPOWERED`, one `TESTED`, 111 `NO_OVERLAP`, and seven `FAILED`; the observed failure fraction is 7/12,475 (0.000561), but the family-level decision is unavailable and must not be inferred from that fraction alone.

## Run handling and consequence

The baseline receipts and worker outputs remain preserved. No family correction, detail materialization, or downstream local-rg result promotion was performed. The baseline run manifest remains at its runner-owned checkpoint value of 743; it was not edited after coordinator failure.

The independently identified roundoff run remains a separate active run. Its runner has the same unique-univariate aggregation assumption. Its receipts must remain separate, and its finalization must be reviewed for the same issue. Do not merge the baseline and roundoff run outputs.

## Required resolution

Before another primary LAVA family is frozen or run, define how the univariate gate is evaluated once per trait-by-locus when pair-specific harmonization produces different variant sets. The method must also meet the frozen univariate coverage/QC requirement, or be explicitly versioned as a new analysis family with a justified and stricter testing policy. Preserve this baseline family's aggregation/QC failure as observed; do not retroactively repair its p-values or relax the 5% cap.

## Demonstrated resolution path

A trait-only `LAVA::process.input` → `LAVA::process.locus` → `LAVA::run.univ` diagnostic was run for insomnia at locus 2 using only the insomnia summary-statistic input and no sample-overlap matrix. It returned one valid univariate result (`p=0.0878159`, 1,260 SNPs, 179 components), distinct from both pair-context p-values above. The reusable, standalone diagnostic is `brain6/scripts/diagnose_lava_single_trait_univariate.R`; it does not modify either frozen run or produce a family result. This established a feasible implementation path for the frozen `unique_trait_by_locus` test unit. At the time this audit was written, full-family coverage and bivariate re-execution remained to be implemented; the follow-up below records the distinct v3 implementation and its current QC status.

## Follow-up: pair-independent canonical v3 production (2026-09-24)

The resolution path above was implemented as a distinct immutable family; it does **not** repair or reinterpret the v2 receipts described in this audit. `brain6/scripts/run_lava_canonical_v3.py` runs each of the seven trait-by-locus univariate tests without pair context, in persistent four-worker chromosome batches of five loci. Family, execution, input, reference, runtime, script, per-cell and output identities are checksum-bound under run ID `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`. The architecture and five-locus batch choice are documented in `qc/lava_canonical_v3_benchmark_20260923.md`; the representative single-trait diagnostic matched the public LAVA path. The immutable v2 baseline and separate roundoff run remain in their original output roots.

At 2026-09-24 05:11 UTC, the read-only v3 verifier accepted 1,357/2,495 locus receipts, covering 9,499/17,465 canonical cells, with zero invalid receipts: 7,643 `TESTED`, 1,856 `NOT_RUN`, and zero `FAILED`. The frozen maximum is 873 untested cells, so the observed lower bound already guarantees that v3 will fail family QC. This is an underpowered local-h² status, not a reconciliation error; no threshold has been changed. The coordinator remains active to finish complete status accounting. A separate finalizer, `brain6/scripts/finalize_lava_canonical_v3.py`, refuses to write a decision until the exact frozen family and aggregate are complete and verified; a live incomplete-run check confirmed it creates no decision. Pairwise local-rg remains unauthorized unless the full family decision passes its frozen QC gate. This follow-up closes the implementation gap identified above while preserving v2's failed-analysis record; it does not establish a successful corrected family.
