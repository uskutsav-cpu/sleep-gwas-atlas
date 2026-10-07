# Reproduce the statistics research package

All commands below run from the repository root. `PYTHON=/opt/miniconda3/bin/python` is the executable used in this environment; NumPy/SciPy/pandas/matplotlib/pytest are required. Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. The experiment uses one process and bounded batches128; output size is small.

1. Inspect `source_acquisition_manifest.json`, `source_acquisition_receipts.json` and `source/LICENSE`. Public sources are already retained; do not reacquire unnecessarily. `acquire_source.py` is the acquisition recipe and would regenerate acquisition timestamps/receipts if deliberately rerun.
2. Inspect `calibration_protocol.json` and `calibration_freeze.json`. The stored freeze must remain unchanged. `calibration.py freeze` deliberately refuses to overwrite an existing freeze.
3. Execute the locked simulation:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/calibration.py run
```

The code verifies the frozen hashes first, executes12 unchanged-upstream numerical fixtures, then24×10000 model-only replicates. It writes the same numerical tables deterministically, with a fresh completion timestamp in the summary. Captured historical execution output is `calibration_execution.log`.

4. Reproduce retrospective empirical and algebraic diagnostics:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/empirical_diagnostics.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/tail_diagnostics.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python brain6/translational_psychiatry_research_v1/research_v1/statistics/plot_diagnostics.py
```

The first command reads the archived block table and source protocol without editing them. It checks all3304 estimates, calculates fixed-threshold SE fragility, verifies scalar-N invariance, and reconstructs patched source hashes in a temporary directory that it removes. Tail diagnostics are explicitly post-result deterministic appendices, not edits to the frozen simulations. The plot contains only synthetic data and is a research diagnostic.

5. Verify meaningful regression tests:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/miniconda3/bin/python -m pytest -q brain6/translational_psychiatry_research_v1/research_v1/statistics/test_statistics.py
```

`regression_tests.log` records9 passed tests. They address source identity, frozen hashes, upstream concordance, N invariance, a deterministic estimator/variance counterexample, Gaussian-product cumulants, the unknown-covariance contrast bound, the unaltered family8465 and unclipped empirical determinant violations. Passing these tests does not count as native GWAS validation.

The optional `check_upstream_identity.py` records a fresh public HEAD observation without cloning. Network failure is recorded rather than interpreted as source identity.

No original empirical GWAS covariance calculation, genotype reconstruction, cohort replication or full native SUPERGNOVA pipeline ran in this package. Exact original inputs and pair receipts remain unavailable. See `SUPERGNOVA_REASSESSMENT.md` for the scientific limits and restart requirements.

Additional work completed after that primary calibration package is documented in `PARENT_ANALYSIS_ADVERSARIAL_REVIEW.md`: separate recomputation of72global contrasts and60molecular conditions,60cross-language retainedRreference-function calculations,twoGaussian quadratic-form checks,and60newfrozen all-palindromic-exclusion molecular conditions. The latter has its own `molecular_palindromic_diagnostic/protocol.json`; do not overwrite it. A storage-only amendment preserves the original code and redirects restricted ADHD-containing filtered slices to ignored local work. They must not enter Git or a public handoff. These use the same inherited sources and are not independent cohort replication.

`INDEPENDENT_NATIVE_FINE_MAPPING_REVIEW.md` additionally covers8newly saved native SuSiE fits,8credible sets and12coloc-SuSiE posterior conditions without refitting. The `test_parent_review.py` record now separates7source-free derived-evidence checks from2actual local-input identity integration checks, preserving the original nine-test record. Select `-m 'not integration'` versus `-m integration` explicitly. Missing actual local inputs fail the requested integration run, and must not be reported as completed science from a source-free checkout. See both review reports for commands and method limits.
