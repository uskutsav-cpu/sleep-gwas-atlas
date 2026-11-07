# One simulated false finding: method-setting sensitivity

We independently reconstructed seed 69001 from the [linkage-aware null simulation](../RESULTS.md), where trait 2 is null everywhere and study errors are uncorrelated. The default calculation again returned one false conjunction-FDR crossing at **0.04247**. The simulated marker's marginal p-values were **0.00514** for trait 1 and **0.000343** for trait 2; it belongs to a truly non-null trait-1 block. The two directional conditional-FDR values were 0.04222 and 0.04247. Thus the crossing is not evidence of trait-2 causality: trait 2 was generated with zero mean at every marker. The conservative max-p multiplied by all 44,000 marker slots is 1 after capping, illustrating the different procedures; it is not a comparison of power or a claim that all source assumptions hold in real GWAS.

| Configuration on this same simulated dataset | Finite values | Crossings at 0.05 | Minimum conjunction value |
|---|---:|---:|---:|
| Default smoothing, monotonic adjustment and genomic correction | 44,000 | 1 | 0.04247 |
| No monotonic adjustment | 44,000 | 1 | 0.04247 |
| No smoothing | 44,000 | 0 | 0.38309 |
| No genomic correction | 44,000 | 1 | 0.03332 |

In this one selected dataset the false crossing depends on enabling smoothing, while disabling the monotonicity adjustment or genomic correction does not remove it. That observation does **not** prove smoothing is implemented incorrectly or that turning it off provides valid FDR control. The configurations are post-hoc, use only one seed, and could behave differently under other architectures. The source code's smoothed lookup and interpolation warrant further inspection and independent implementation comparison; the mechanism behind the wider simulation results remains unresolved.

The [pre-run protocol](PROTOCOL.md), [code](run.R), [setting results](sensitivity.tsv), [simulated crossing](default_hits.tsv), [log](run.log) and [verification](verification.json) retain the test. The log records warnings about discarded row names; their effect is unverified. This is a software diagnostic, not a real sleep–frailty association or a validated replacement for the failed LAVA analysis.

A [ten-dataset paired extension](../smoothing_replicates/RESULTS.md) found no unsmoothed crossings in this same selected simulation architecture, compared with six default runs with crossings. It does not establish unsmoothed FDR control or statistical power.

A [selected-marker tail audit](sparse_tail/RESULTS.md) finds only 0–2 pruned markers at least as significant on the conditioning trait in each iteration, and none meeting both exact target thresholds. These descriptive counts are not the package's actual grid-cell counts or a proof of the cause.
