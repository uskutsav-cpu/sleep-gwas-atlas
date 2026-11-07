# Registry-insomnia source-shift stress test

The multi-start model was trained on other cells in the 12×8 UKB-oriented sleep/brain/frailty matrix, with each tested UKB-insomnia–FI or UKB-insomnia–FFS cell held out. We compared its ten existing masked-cell predictions for each outcome with the separately estimated FinnGen R13 registry-insomnia correlation. **The model has not learned an R13 insomnia embedding**; this comparison asks what happens if the UKB-insomnia representation is used as a proxy. It is not a fully independent validation because the FI and FFS outcome GWAS are reused and the two sleep definitions, populations, variant supports and ascertainment differ.

| Outcome | UKB-insomnia correlation used as model target | R13 registry-insomnia correlation (SE) | Mean multi-start prediction; range across 10 dependent masks |
|---|---:|---:|---:|
| FI | 0.6396 | 0.4263 (0.0389) | 0.4916; 0.4496–0.5350 |
| FFS | 0.4833 | 0.2445 (0.0371) | 0.6299; 0.5735–0.6858 |

The FFS model predictions are noticeably higher than the R13 estimate. This shows a limitation of treating one sleep source as a substitute for another; it does not isolate the causes of the difference, quantify population transfer, or prove that any estimate is wrong. In particular, training did not include the held-out UKB-insomnia–FFS value but did retain other correlations involving the same sleep trait and frailty outcome. The ten masks reuse the same source estimates, and are not independent replications. No uncertainty interval for the predictions or statistical significance of their difference from R13 is asserted.

The existing column-mean, additive and original single-start predictions are retained for context in `registry_transfer_summary.tsv` and `registry_transfer_predictions.tsv`; they were developed or reviewed before this source-shift analysis and are not independent model-selection evidence. The R13 point estimates use FI, FFS and insomnia on 1,134,953 common SNPs; see [the paired sensitivity audit](../frailty_sensitivity/FRAILTY_COMPARISON.md). FFS and FI include sleep-related items or exhaustion, cohort overlap is not fully known and the Finnish reference fit requires scrutiny. A meaningful generalization study would first harmonize definitions and alleles, verify overlap and population-matched LD, freeze the model and metrics, and evaluate genuinely independent sleep **and** frailty sources. No AlphaGenome outputs were used to fit this model.

The [pre-result protocol](REGISTRY_TRANSFER_PROTOCOL.md), `registry_transfer.py` and `registry_transfer_receipt.json` specify the saved predictions, exact LDSC reference files and checksums. The original failed LAVA gates remain unresolved.
