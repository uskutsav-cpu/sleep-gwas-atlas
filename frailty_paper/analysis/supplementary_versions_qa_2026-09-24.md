# Supplementary table and software-version QA

Review date: 2026-09-24.

## Table rebuild

Rebuilt the available supplementary tables with the pinned Python 3.11.11 environment. The replication audit contributes 15 candidate/source rows to S8. S16 retains the historical 91/91 and 94/94 validation records, the latest integrated preflight snapshot (97/97 tests; 351/351 registered resource files), and the separate 97/97 package-suite run after the S16 builder/regression update. It preserves the captured run identities rather than relabeling earlier records.

The S1 source is the 32 rows of `all_acquired_resources.tsv` whose `resource_type` is `GWAS summary statistics`. Provenance now hashes precisely that consumed subset, rather than unrelated resource-manifest rows. This permits reproducing S1 without tying its checksum to unrelated manifest metadata edits.

## Current test evidence

The standalone mixed-worktree suite snapshot ran in Python 3.11.11 and passed 95/95; its scope caveat is recorded in `current_worktree_test_validation_2026-09-24.json` and `current_worktree_unittest_2026-09-24.log`. A fresh integrated pinned preflight from source HEAD `6b906ba` at 07:01–07:03 UTC passed 97/97 tests, verified 351/351 registered resources, and passed the missing-abstract source audit and 28-claim manuscript audit; its captured output is `preflight_current_head_py311_2026-09-24.log`. The post-update standalone package suite passed 97/97; output is `current_head_py311_unittest_2026-09-24.log`. The earlier 06:21–06:26 capture remains at `preflight_post_abstract_audit_py311_2026-09-24.log`. Both preflight runs used a mixed worktree and are not clean-checkout attestations.

## Row counts and provenance

- S1 GWAS metadata: 32
- S4 SNP heritability: 20
- S5 available global rg: 168
- S7 latent-factor rg: 84
- S8 replication/source audit: 15
- S15 correction sensitivity: 96
- S15 sensitivity conclusion matrix: 119
- S16 software/resources/versions: 50
- Availability index: 16 table rows

All 24 recorded input checks passed, including the scoped manifest subset hash, latest integrated-preflight log, and post-update package-test log. All 9 generated-table/index output hashes match `supplementary_tables.provenance.json`. The table counts match the provenance record. S16 contains 50 rows and retains the earlier 94/94, 350-resource run as historical evidence while separately recording the 97/97, 351-resource preflight and post-update suite.

## Limits

This rebuild does not complete the supplementary package. Review-dependent tables S2–S3, source-provenance-blocked S6, and unsupported/not-justified downstream locus tables S9–S14 remain unavailable or gated. S16 remains partial because it does not represent every gated/future tool execution or the final end-to-end replay.
