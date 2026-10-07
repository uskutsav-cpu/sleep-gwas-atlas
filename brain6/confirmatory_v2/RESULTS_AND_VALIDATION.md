# Results and validation

Table audit identity: brain6-continuation-20261007-v1. Frozen protocol SHA256: d160e34f009543c650e894c035a8d30f0331572df92fc9e0231f1dd55643cf16. Starting public main commit: 303c905c229104b7bf623df6de7431abb87e0727; current remote was checked and remained this commit.

The original package has 77 available SHA256 matches, one missing original canonical source and zero mismatches. The original package auditor was run and correctly fails on that missing source; its failure is retained. Our audit reports PASS_AVAILABLE_TABLES_ORIGINAL_INPUT_REPLAY_INCOMPLETE. Neither duplicate exported tables nor a new validator is a replacement for the missing provenance source.

Actual independent numerical checks:

- Recounted 72 global rows/35 inherited significant, 25 candidates/20 regions and all local-method statuses.
- Replayed all 3,304 SUPERGNOVA covariance P values from rho/sqrt(var); max absolute difference 3.33e-16. Recomputed Bonferroni family values, max error 4.41e-13; recomputed same-pair protected interval overlaps: zero. Missing correlation values stay missing.
- Analytically reweighted all six H0–H4 vectors at lower p12; posterior agreement within 1e-10. This checks posterior arithmetic, not Bayes-factor derivation from raw inputs or the one-signal assumption.
- Recomputed BH in the original 54-tissue and 1,680-pathway families and empirical P=(extreme+1)/10001. Max BH difference <1.4e-12; no significant result. This checks archived permutation counts, not rerunning 10,000 null sets.
- Recovered source-first failure arithmetic: 2,847 recoveries needed; long-sleep alone leaves 1,291 NOT_RUN; unique smallest hypothetical perfect repair set is long sleep + insomnia + PD + MDD, leaving 533. No hypothetical recovery was admitted.
- Re-examined recorded LD factor/normalization audits. A diagonal reset without corresponding off-diagonal covariance normalization explains the mathematical scale failure in 19 diagnostic blocks. A normalized Gram matrix is PSD by construction; this does not validate genotype identity, alleles, cohort reference, or the originally exported LD. No empirical corrected LD or native SuSiE was run.

New diagnostic: 723 of 2,199 numeric local correlation estimates fall outside [-1,1]; 1,105 are absent. One significant chr11 insomnia–ADHD block has derived corr 1.103737. Its covariance test is retained with explicit inferential limits; no clipping or biological amplification claim.

Historical secondary results (GRCh37; same-source covariance, not replication):

| Pair | Interval | Covariance | Variance | P | Family-adjusted P |
|---|---|---:|---:|---:|---:|
| Insomnia–ADHD | chr11:112459489–114257728 | 0.0002691266 | 2.306485e-9 | 2.097239e-8 | 0.0001775313 |
| Insomnia–MDD | chr6:100630147–102636772 | 0.0001160579 | 5.911901e-10 | 1.813085e-6 | 0.0153477611 |
| Insomnia–MDD | chr11:112459489–114257728 | 0.0002255223 | 1.262670e-9 | 2.200532e-10 | 1.862750e-6 |

Software validation: 410 source-free contract tests pass with zero skips; 49 loaded test objects explicitly require archived sources/native runtimes (wildcard modules may contain more individual cases). The initial unpartitioned CI run failed 12 tests and errored on 39; this is an engineering availability problem, not evidence of scientific failure or success. Initial Brain6 pytest run: 10 failed, 534 passed, 3 skipped. Its ten real-data tests remain integration requirements; native tests remain separately labeled. The synthetic zero-family fixture now mocks available disk space; production disk checks and the insufficient-space test remain unchanged. No original assertion was removed. Python compilation, shell syntax, panel lock validation, phase0 strict declaration audit and synthetic end-to-end smoke pass. Phase0 reports 0/45 real harmonization/LDSC-ready traits in the new checkout; synthetic smoke has no biological meaning. Broad atlas report: 3/23 gates.

Final Brain6 source-free run: 544 passed after the claim-regression tests were added; 14 explicitly deselected integration/native cases. Eleven continuation tests pass. See qc logs and machine test partitions for numerical test counts, CI execution status and unresolved native checks. No failed native analysis was recast as a pass. The original full command remains `python -m unittest discover -s tests -v`; it must be run after exact archives/runtimes are restored. Source-free CI is a separate named gate.
