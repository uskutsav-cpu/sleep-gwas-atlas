# Jansen UKB insomnia native-N pilot v2: terminal execution hold

**Status:** `FAILED_PREDECLARED_ZERO_FAILURE_CONDITION`. V2 was frozen and committed in `225e89ac` before materialization or outcomes. The v1 freeze was unused and has no output root. V2 materialized 88 loci and 207,643 rows by retaining every original receipt-verified SNP/A1/A2/Z row and replacing only the constant N with exact Jansen source per-SNP N. Four disjoint 22-locus workers launched from the dedicated SSD root. No canonical, v2, roundoff, protected candidate, or other production output was written by this pilot.

The predeclared rule required 88 receipt-bound outcomes and **zero execution failures** before any TESTED-fraction or low-local-h2 comparison. Worker 2 processed 22 loci but exited 1 at its final validation, `Error: Execution failure; do not promote`. Its R runner intentionally refused to write a TSV or summary because at least one per-locus result had status `FAILED`. Consequently, the aggregate and advancement decision files were **not produced**. The failed locus/reason cannot be recovered from the persisted worker-2 log; that is a diagnostic reporting limitation, not grounds to relax the criterion or infer its result.

| Worker | Exit | TSV/summary | TESTED | NOT_RUN | FAILED | Strict-gate pass |
|---|---:|---|---:|---:|---:|---:|
| 1 | 0 | yes | 12 | 10 | 0 | 0 |
| 2 | 1 | no | unknown | unknown | at least 1 | unknown |
| 3 | 0 | yes | 13 | 9 | 0 | 0 |
| 4 | 0 | yes | 12 | 10 | 0 | 0 |

The 66 retained worker-1/3/4 outcomes total 37 TESTED and 29 NOT_RUN. They are partial technical results and cannot be used to choose a model, claim a full-family improvement, or advance to a 2,495-locus screen. The source stays `ADMISSIBLE_SENSITIVITY_ONLY`; confirmatory/bivariate promotion remains withheld for the separate model and overlap reasons in the frozen protocol.

## Exact pilot receipts

SSD root: `/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2`

- Frozen config: SHA-256 `893b0b8251805a421355e3141b658c1997bed4783daf7990a378abeb32796401`.
- Frozen protocol: SHA-256 `4391e73face9c18279da7647e5106058ccbb778b6acf647b598a3bcb55d3e022`.
- Materialization receipt: SHA-256 `24ce34356859f1a3f58390c1dc9193c6849cabe2944b1ac947b08d23f53e4f83`.
- Launch receipt: SHA-256 `b5a3b7e0d60f43a38611150956e0fa6b302a40197885a88c16a7ea3d30039060`.
- Exit receipt: SHA-256 `4e3c578c526353d3a0ef8108b9e2915f8e2a47c599e0b30f0358b7a7bcfe1aa1`; exit codes `[0,1,0,0]`.
- Worker-2 log: SHA-256 `0bd8e5f5f28ca35e58d6fdb5764b0bee54e0d0f3d93d0985db3fa57ba8336830`.

A future *diagnostic* rerun, if separately frozen, would need to persist per-locus failed reasons before failing terminally. It must use a new version/output root and cannot retroactively change this pilot's rule or status. No such rerun is claimed here.

## Subsequent diagnostic addendum

A separately frozen v3 worker-2-only postmortem reproduced one failure at
locus 950; its original and native-N shards each contain zero rows. See
[worker-2 postmortem result](INSOMNIA_NATIVE_N_WORKER2_POSTMORTEM_V3_RESULT.md).
This explains the failure mechanism but **does not reopen or pass** v2's
predeclared zero-failure advancement criterion.
