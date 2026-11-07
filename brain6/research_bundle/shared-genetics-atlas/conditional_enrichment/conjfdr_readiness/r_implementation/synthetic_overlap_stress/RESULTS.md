# Synthetic one-trait-only and study-correlation stress test

In 10 paired simulated datasets, all 44,000 variants were null for the second trait; 5% had a fixed association with the first trait. Using the locally patched R implementation and a no-LD synthetic reference, **none of the 10 runs crossed conjunction FDR 0.05** when study-error correlation was zero, and **none of the 10 crossed** when it was set to 0.3. All returned values were finite and in [0,1].

| Simulated study-error correlation | Runs with any threshold crossing | Smallest conjunction value across runs | Median of run-wise minima |
|---|---:|---:|---:|
| 0 | 0 / 10 | 0.0793 | 0.3785 |
| 0.3 | 0 / 10 | 0.0652 | 0.1313 |

The paired runs reused the same underlying random errors and non-null locations; their minimum value was lower at correlation 0.3 in 8 of 10 pairs. These descriptive results do **not** prove that overlap causes a change in real data or that the method controls FDR. With zero events in only ten independent simulated datasets per scenario, an exact one-sided 95% upper bound for the probability of at least one crossing in this specific simulation is about 25.9%; the test cannot rule out a rate well above 5%. Correlation 0.3 is illustrative, not an estimate of actual overlap.

The [pre-run protocol](PROTOCOL.md), [simulation code](run.R), [per-run results](runs.tsv), and [log](run.log) preserve assumptions, seeds and settings. No real GWAS or AlphaGenome outputs were analyzed. The scenario lacks LD, realistic polygenicity, varying effects and sample sizes, source precision, actual phenotype overlap, and uncertainty in correlation. The reference has no edges, so one pruning iteration suffices for this artificial test but does not approximate an actual genome-wide run. Package warnings are retained in the log. Results should not be extrapolated to the real FI or FFS analysis. The package's sample-overlap correction remains unimplemented; a full realistic calibration and independent validation are still needed.
