# LAVA canonical v3 locus-worker diagnostic

Date: 2026-09-24. This is a one-locus implementation check, not a family result.

## Procedure

- Input: frozen locus 464 (chromosome 3), the seven checksum-bound materialized sum-statistic files, pair covariance matrices from the v3 family lock, and official UK Biobank LAVA v1.1 chromosome 3 reference.
- Runtime: pinned R 4.3.3 and LAVA 0.1.5; roundoff-only block-matrix symmetrization is enabled under `LAVA015_BLOCK_REDUCED_SYMMETRY_V1`.
- The v3 worker recomputed the seven single-trait LAVA tests, then applied the frozen strict gate `P < 0.05/17,465 = 2.86286859433152e-06` before preparing any pairwise inputs.
- The diagnostic configuration and outputs are isolated under `/private/tmp/lava-v3-locus464-config5.json` and `/private/tmp/lava-v3-locus464-diagnostic5/`; no production, v2, or roundoff output was written.

## Result

All seven worker single-trait P values exactly match the strings in the separately generated, receipt-verified canonical production row for locus 464. The production receipt passed the coordinator's full `verify_locus` check. At this locus, insomnia has `P=2.98078e-05`, above the frozen gate, so all five pair slots are retained as `UNIVARIATE_UNDERPOWERED` with reason `ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_FAILED_FROZEN_GATE`. No bivariate call ran for this locus.

The preserved v2 checkpoint had marked insomnia–ADHD at locus 464 `TESTED` using pair-specific univariate inputs. That result remains unchanged in v2; it is not promoted into v3 because the canonical insomnia test fails the v3 family gate. This illustrates the intended separation between pair-harmonized bivariate data and the unique trait-by-locus testing family.

The worker returned all 7 univariate rows and all 5 pair slots; zero worker failures occurred. This diagnostic does not exercise the `run.bivar` path because no pair passed the gate, and it does not establish whole-family eligibility or QC.

## Provenance

| Artifact | SHA256 |
|---|---|
| `brain6/scripts/run_lava_family_locus_v3.R` | `43341347166a252d7b87aebbe661fe236535a575f8f84f6f1fce5f6e849998fe` |
| `brain6/config/lava_family_canonical_v3.json` | `f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417` |
| `brain6/config/lava_execution_canonical_v3.json` | `f37a212ff44472004549095a3b17d28da0d130e2ab2a9ea507feb87312b65d27` |
| Diagnostic config | `1b78f1e7e891aa82f26ce1339bc2d2b10706b34e75c2d0cbcf1560dc395c4f57` |
| Diagnostic univariate output | `1c85c625caf5dc4c7d3cc911bcf7e96bf4b577238f4085e3cc7afb5e0c0a56d4` |
| Diagnostic pair output | `4b96d18cecf6dc5a29d98e6b49794e2acaf31d3352680a74123d60b0b56c01bf` |
| Canonical production row | `79f1be728b983c9fd6b9b7c9781a71c575728b1685fe34b6e864ec041a5f43c7` |
| Production locus receipt | `58c5fac14ba406d75b571920c4159df29ff61c33d451c4052f7fe148210a15ed` |

## Validation

The R worker parses successfully. `pytest brain6/scripts/tests/test_lava_canonical_v3.py -q` passes 15 tests, including a new source-order assertion that canonical univariate results and the strict gate precede pair input preparation and bivariate execution. The current full 17,465-cell run remains in progress; no family decision or downstream output is authorized until its frozen coverage/QC gate passes.

## Pairwise-stage optimization and fail-closed check

The pairwise worker was then changed to receive the seven canonical rows from the verified source receipt instead of recomputing the seven trait-only tests. It now loads pair reference and summary statistics lazily, only after a pair passes the gate. A second isolated locus-464 run supplied the production cell hash and receipt hash; its status record echoed both hashes and the canonical and bivariate run IDs, all seven P-value strings still matched production, and all five pair slots remained underpowered. No reference or summary-statistic load was needed because no pair passed. The outputs match the first diagnostic exactly.

A separate bivariate v3 lock now freezes the inherited LAVA 0.1.5 bivariate settings, the exact 12,475-slot family and BH denominator, the 1% execution-failure limit, the 5% canonical untested limit, four concurrent workers, and a new output root. It binds the canonical production run ID, both canonical lock hashes, pair-selection/overlap provenance, and current R worker source hash. The new coordinator re-verifies every canonical receipt and the full-family `qc_pass` audit before it creates an output directory. An actual invocation against the currently active canonical run stopped with “no completed full-family audit,” and `work/lava-local-rg-v3` remained absent.

Review of the 1% QC accounting found that the coordinator was incorrectly counting `NO_OVERLAP` as an execution failure. The lock calls this limit an execution-failure ceiling, and existing Brain6 audits record `NO_OVERLAP` separately from `FAILED`; `NO_OVERLAP` therefore remains a valid ineligible slot in the fixed BH family but does not consume the execution-failure allowance. Updated the coordinator so only `FAILED` and `NOT_RUN` count toward that ceiling, while any no-overlap slots remain reported and yield `INSUFFICIENT_EVIDENCE`. The frozen threshold was not relaxed, and no pairwise results existed when this correction was made.

The bivariate coordinator also reconstructs the exact ordered 17,465-row canonical aggregate from the verified per-locus cell outputs and compares every field before the pairwise gate. This closes a gap where the prior check compared the aggregate file hash with the audit record but did not independently compare its content to the receipts. The checksum-bound pair coverage parser now rejects negative counts, duplicate/missing keys, or statuses that disagree with the locked minimum-shared-variant threshold. A read-only validation of the actual table confirmed all 12,475 pair-locus rows, 0–3,855 shared variants, and consistent `READY`/`NO_OVERLAP_LT_MIN_K` labels. Regression tests cover valid and corrupted aggregates and malformed coverage. The focused canonical and bivariate tests pass 24 tests; the full Brain6 script suite passes 80. Python byte-compilation, R parsing, lock validation, and `git diff --check` pass. The bivariate output root remains absent and no pairwise production started while canonical-family QC is incomplete.

| Current pairwise-stage artifact | SHA256 |
|---|---|
| `brain6/scripts/run_lava_family_locus_v3.R` | `5049eb6491b55c89095d89135299d5065effcb07512d1a92afce46ceee3a4645` |
| `brain6/scripts/run_lava_bivariate_canonical_v3.py` | `1acbe70c4e6209141c76bded30aeb1cb5f648077be6f9b72c26a91420dd9c4e3` |
| `brain6/scripts/tests/test_lava_bivariate_canonical_v3.py` | `7cd32ba2d0018372ab9a0d2e74f8afed801827807f61d341ac228627fdc8d1cc` |
| `brain6/config/lava_bivariate_canonical_v3.json` | `91bc2f11abd4d594eb6132fb570ad008ddc01cdfe85b574e673fb76ed359dcc4` |
| Diagnostic config with source hashes | `960bf8b4b47077b0575cdd08ed4720cb05e45c94b7dbacfb465e8fac405de616` |
| Optimized diagnostic status | `6014d1e51eb18b6fee0aadda44fc3ad9a859dd1efba050e932000f20f543b8ee` |
