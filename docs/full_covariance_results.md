# Full 45-trait covariance result

The complete locked-panel GenomicSEM LDSC run finished on 2026-08-28 UTC in
38 minutes 47 seconds. It used the exact ordered `atlas-v1.0` panel, the pinned
European LD-score reference, R 4.3.3, GenomicSEM 0.0.5 at commit
`6b65ca5db39fdade08b0d811477be1cdd57b5039`, and 1,082 block-jackknife blocks.

The validated export contains 45-by-45 genetic-covariance, genetic-correlation,
and intercept matrices, the 1,035-by-1,035 sampling-covariance matrix, all 1,035
lower-triangle estimates, trait scales, pinned run metadata, the native RDS, and
matrix diagnostics. `scripts/26_validate_covariance.py` fails closed on the
dimensions, ordered IDs, completeness, finiteness, symmetry, metadata, and
explicit covariance-versus-correlation uncertainty fields.

## Diagnostics and interpretation constraints

- The genetic covariance matrix has three negative eigenvalues; its minimum is
  -0.00561345. The standardized genetic-correlation matrix likewise has three
  negative eigenvalues, with minimum -0.0952280. Sampling error can make an
  estimated LDSC matrix indefinite, so downstream SEM must report any smoothing
  or near-positive-definite correction and include the unsmoothed result as a
  sensitivity analysis.
- The sampling covariance matrix is finite and positive definite, but its
  absolute condition number is approximately 1.46 billion. GenomicSEM itself
  warns that the 1,082 blocks required for this 45-trait system can induce
  block dependencies and biased downstream test statistics. Q_SNP and model
  fit must therefore be treated as high-block-count sensitivity results.
- Relative to the earlier 396 sleep-by-non-sleep Python LDSC family, median
  absolute differences are 0.00157 for genetic correlation and 0.00137 for its
  SE. The largest differences occur for weak melanoma estimates.
- HDL is an implementation-sensitive QC case. Standalone Python LDSC invoked
  its two-step estimator and reported intercept 1.196, while GenomicSEM and the
  Python pairwise runs report intercepts above 1.20 (GenomicSEM 1.3435). A
  targeted GenomicSEM rerun with 200 rather than 1,082 blocks reproduced 1.3435,
  ruling out block count as the cause. The original Phase-1 verdict is retained
  for auditability, while HDL is excluded from confirmatory SEM under the same
  predefined intercept threshold.
- T2D remains excluded for intercept above 1.20 and melanoma for h2 Z below 4.
  Multiple sclerosis remains included with an explicit nonphysical
  liability-h2-above-one warning rather than being silently removed.

Run the structural validator and derive the SEM input ledger with:

```bash
python3 scripts/26_validate_covariance.py
python3 scripts/27_genomicsem_trait_qc.py
```
