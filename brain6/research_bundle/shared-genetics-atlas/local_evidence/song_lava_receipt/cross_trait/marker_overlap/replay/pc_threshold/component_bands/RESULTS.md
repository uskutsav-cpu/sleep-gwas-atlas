# Where the LAVA variance estimate changes

In the 99.5%-retention replay, cumulative squared FI projections at 73, 107, 142, 216 and 306 components agree with the five separate threshold runs to within 10⁻¹². This supports comparing the added component bands in this particular replay.

The first 142 components, retained at 98%, have a combined squared projection of 0.0010156. The 74 components added to reach 99% contribute 0.0002162, or **2.92×10⁻⁶ per component**. At the 99% setting the estimated sampling term is **5.71×10⁻⁶ per component**. Adding that band therefore lowers the average projected signal below the sampling subtraction, making the FI variance slightly negative. The next 90 components, added at 99.5%, have a mean squared projection of 4.65×10⁻⁶, again below that setting's sampling term. The [complete band table](bands.tsv) includes all five bands.

This describes the numerical mechanism of the sign change, not its statistical cause. Components are combinations of markers, not individual variants. The sampling term is estimated separately at each threshold, and these post-selected bands are not independent experiments. Sampling variation, GWAS–reference LD differences and other input/model issues remain possible; these values are not calibrated tests or evidence for dropping particular components. The original LAVA family gates remain failed, and no bivariate correlation was calculated.

An [idealized noise comparison](idealized_null/RESULTS.md) shows that the total 216-component variance estimate is close to its numerical boundary and describes the selected band under strong hypothetical Gaussian assumptions. Those probabilities are not calibrated tests of the actual data and do not establish the cause of the failure.

[Protocol](PROTOCOL.md), [replay code](run.R), [squared projections](projections.json), and [independent arithmetic verification](verify.py) document the calculation. The full reference and source limitations from the [parent replay](../../RESULTS.md) apply; reproduction does not establish historical input identity or valid local inference.
