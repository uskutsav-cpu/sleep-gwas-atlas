# Pair-independent LAVA gate probe (diagnostic only)

This two-cell diagnostic probes the ambiguity between the frozen unique trait-by-locus univariate gate and pair-context univariate values. It is a new, explicitly non-promoted diagnostic; it does not replace the failed canonical v3 result, repair the v2 or roundoff family, or establish a complete local-rg family.

The frozen v3 trait-only gate uses 17,465 tests and requires both traits' locus-specific p-values to be strictly below `2.86286859433152e-6`. Applying that rule across the five preselected local pairs admits exactly two pair-locus cells, both at locus 2207: long sleep–SCZ and long sleep–Parkinson. The probe then uses the existing checksum-bound pair inputs, overlap matrices, LAVA 0.1.5, roundoff-only numerical patch, and locked bivariate parameters to run those two bivariate calculations even though the pair-context univariate p-value for long sleep exceeds the gate.

| Pair, locus 2207 | Pair-context SNPs | Canonical trait-only p-values | Pair-context univariate p-values | Local rg (95% CI) | p | Full 12,475-slot BH q |
|---|---:|---|---|---|---:|---:|
| Long sleep–SCZ | 2,397 | long sleep 8.08649e-7; SCZ 7.98000e-7 | long sleep 2.38316e-3; SCZ 8.30126e-7 | 0.387708 (−0.11166, 1.00000) | 0.120592 | 1.000000 |
| Long sleep–Parkinson | 69 | long sleep 8.08649e-7; Parkinson 1.02223e-18 | long sleep 3.83645e-6; Parkinson 1.02223e-18 | 0.884357 (0.62499, 1.00000) | 2.21035e-5 | 0.275741 |

The family correction uses the fixed 12,475 pair-by-locus denominator, assigning p=1 to the other slots that did not pass the pair-independent univariate gate. Neither probe result passes that correction. The Parkinson result is a nominal signal only and remains exploratory; the wide interval and tiny shared-reference set (K=15) limit interpretation. The SCZ estimate is not nominally significant. These two cells do not rescue the roundoff family's aggregation failure or the canonical v3 family QC failure.

The contrast itself establishes why a cosmetic aggregation fix is invalid: at the same locus the pair-context long-sleep p-value differs substantially by partner, as the pair-specific shared-SNP sets differ (2,397 vs 69). The probe is a targeted diagnostic of the unique-gate interpretation, not a complete or promoted sensitivity family. No downstream locus analyses are eligible from these results.

Reproduction uses `brain6/scripts/probe_lava_unique_gate_v1.R` with the receipt-bound roundoff locus-2207 worker config, canonical v3 family lock and aggregate shown in `lava_unique_gate_probe_v1.provenance.json`. The script refuses to overwrite its output and fails if the frozen unique-gate candidate set changes.
