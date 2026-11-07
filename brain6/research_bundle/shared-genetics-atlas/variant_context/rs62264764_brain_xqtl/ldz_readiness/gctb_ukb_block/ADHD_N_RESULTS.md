# ADHD sample-size and source-alignment follow-up

**The strongest residuals concentrate among variants with lower audited sample sizes, but their cause is not established.** The previous analysis used one scalar effective sample size for all 1,547 markers. Repeating the same diagnostic with the existing minimum `N_AUDIT` (64,699 instead of the median 128,214) barely changed its largest residual or the descriptive outlier counts:

| LD reference | Maximum residual, median N | Maximum residual, minimum N | Count above 5, either N |
|---|---:|---:|---:|
| GCTB, eigenvalue clip | 14.2498 | 14.2495 | 10 |
| GCTB, identity shrinkage | 14.2449 | 14.2446 | 10 |
| 1000 Genomes | 11.2733 | 11.2730 | 6 |

Replacing one scalar with another does not model per-variant samples. The minimum-N sensitivity, its three complete variant-level outputs and source-alignment check follow the [protocol](ADHD_N_PROTOCOL.md).

For the ten largest GCTB eigen-clipped residuals at median N, each original ADHD record uniquely matched the stored rsID, chromosome, GRCh37 position, allele pair and ALT-aligned effect and standard error; the input file matched its recorded SHA-256. The ten variants' audited effective sample sizes range from approximately **72,984 to 74,970**, compared with a median of **128,214** over the 1,547 markers. Absolute differences between source control ALT frequency and GCTB ALT frequency had a median of 0.0124 and a maximum of 0.0177. Source control and reference frequencies are not measurements from the same population, and close frequencies cannot establish correct LD or rule out all allele problems. The [source audit](adhd_source_audit.json) gives aggregate results without redistributing original ADHD source rows.

After observing this pattern, we explored a **post-hoc** threshold of N_AUDIT below 80,000. Under both GCTB adjustments, 10 of 58 lower-N variants and 0 of 1,489 higher-N variants had absolute standardized residual above 5. Under the 1000 Genomes reference, the corresponding counts were 5 of 58 and 1 of 1,489. These markers are correlated; the cutoffs and grouping were selected after inspecting results, so these are not independent tests, calibrated probabilities, or evidence that low N causes the residuals. [All counts](adhd_coverage_strata.json) are provided. Variable study coverage and participant overlap are plausible hypotheses for further source-level investigation alongside LD mismatch, association-model differences and other possibilities; per-marker N does not reveal pairwise sample overlap.

No variants were removed or flipped, and no new fine-mapping, colocalization or LAVA conclusion is admitted. The [lower-N results](adhd_min_n_summary.tsv), [variant-level values](adhd_min_n_details.tsv), [code](adhd_n_sensitivity.R), [source-audit code](adhd_source_audit.py), [stratification code](adhd_coverage_strata.py) and [arithmetic verification](adhd_verification.json) document the analysis. The source GWAS is limited to noncommercial research; see the original project source documentation for applicable reuse conditions.

A subsequent [study-direction audit](ADHD_DIRECTION_RESULTS.md) found that all ten GCTB outliers share a question mark at the same source direction-string position, as do 38 non-outliers. It is a possible study-coverage signature, not a verified cohort identity or explanation of the LD–Z mismatch.
