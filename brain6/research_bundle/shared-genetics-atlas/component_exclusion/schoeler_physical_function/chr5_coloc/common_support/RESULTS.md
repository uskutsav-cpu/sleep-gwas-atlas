# Chromosome 5 colocalization on common markers

Restricting both analyses to the same 2,359 previously matched markers changes the post-selected single-signal results very little. The insomnia input lost 92 matched markers and the ADHD input lost one. The outcome values and alleles were identical on the intersection.

| Comparison | H4, p12=10⁻⁵ | H4, p12=10⁻⁶ | H2, p12=10⁻⁵ / 10⁻⁶ |
|---|---:|---:|---:|
| Physical function–insomnia | 0.497 | 0.090 | 0.443 / 0.801 |
| Physical function–ADHD | 0.533 | 0.102 | 0.411 / 0.790 |

H4 is the single-signal model's posterior for one shared causal variant; H2 is association with the insomnia or ADHD trait alone. This same-marker sensitivity does not make the datasets independent or eliminate uncertainty due to phenotype differences, source processing, selection, multiple signals, physical-effect prior and overlap. The earlier conclusion remains: these prior-sensitive results do not establish shared causality or a sleep-excluded frailty mechanism. H4 values for the two traits are not formally compared.

The [protocol](PROTOCOL.md), [code](run.R), [marker positions](positions.tsv), and [posterior table](posteriors.tsv) preserve the analysis. The [independent check](verification.json) confirmed the previously recorded input hashes, exact position intersection and direct-ABF arithmetic (maximum difference 1.78×10⁻¹⁵); it is not a validation of model assumptions. See the [parent report](../RESULTS.md) for original source, matching, prior sensitivities and limitations.
