# Hospital-frailty follow-up of selected marginal SNPs

All 368 SNPs selected by the earlier insomnia/FI/FFS marginal screen had one chromosome- and allele-compatible HFRS lookup in the existing filtered HapMap3 source. Selection and limitations are frozen in [PROTOCOL.md](PROTOCOL.md). The earlier screen used a 3-trait cutoff; here, the four marginal values are compared descriptively with the stricter 4-trait scale 0.05/(4×1,138,262) = 1.098×10⁻⁸. The first three values are conservative bounds from rounded Z-scores; HFRS uses its displayed source p-value.

| Crossing pattern among the 368 selected SNPs | SNPs |
|---|---:|
| FFS alone | 288 |
| FI alone | 19 |
| Registry insomnia alone | 23 |
| FI and FFS | 3 |
| None at the stricter four-trait threshold | 35 |
| Any pattern including HFRS | 0 |

The smallest HFRS p-value among these selected SNPs is 1.47×10⁻⁷ for rs11764337, which is below 0.05 but above the four-trait threshold. These are SNP counts, not LD-independent loci; 35 previously selected SNPs no longer cross because the cutoff changed. This is a selected lookup: it does **not** evaluate HFRS-only associations outside the original 368 markers or exclude shared causal signals that fail any marginal threshold. Unmatched or differently represented variants could affect conclusions; matching rsID, chromosome and alleles across builds is not a full reference audit. HFRS and R13 insomnia share FinnGen participants, FI and FFS overlap UK Biobank, and different ascertainment and power can change significance patterns. “Alone” describes a threshold pattern, not proven trait specificity.

`selected_with_hfrs.tsv` retains the four p-values, all allele strings, HFRS beta before and after alignment, and matching statuses. `summary.json` contains input hashes, counts and the threshold; `lookup.py` reproduces the lookup from the existing external source. This is not an independent two-trait replication, a PLACO test, enrichment-based conjFDR, or a component-excluded frailty analysis. The separate [insomnia–HFRS conjunction screen](../../hfrs_conjunction/RESULTS.md) evaluates its larger common HapMap3 pair and uses a different correction family.
