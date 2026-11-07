# UKB hospital frailty–FinnGen insomnia: chromosome 7 coloc pilot

A subsequent [common-marker insomnia-source sensitivity](insomnia_source_sensitivity/RESULTS.md) obtains a different, prior-sensitive questionnaire-insomnia posterior. Questionnaire insomnia and HFRS share the UK Biobank cohort, so the model's independence assumption is a serious unresolved limitation.

In the previously selected chromosome 7 region, a single-signal colocalization model favors **distinct dominant HFRS and insomnia signals** over one shared signal under the tested settings. This is exploratory evidence, not proof that the traits have different causal variants: the model can miss secondary shared signals, and the HFRS effect-size scale and appropriate prior are not independently verified.

The model used **1,212 matched variants** in the GRCh37 chr7:113571035–114571035 window. For an illustrative HFRS raw-effect prior SD of 0.15, insomnia log-odds prior SD of 0.2, and p1=p2=10⁻⁴:

| Shared-variant prior p12 | H1: HFRS only | H3: distinct signals | H4: shared signal |
| --- | ---: | ---: | ---: |
| 10⁻⁵ | 0.3157 | 0.6335 | 0.0508 |
| 10⁻⁶ | 0.3309 | 0.6638 | 0.00532 |

H0 and H2 are small under those specific settings. H1 retains substantial probability, consistent with uncertainty in the insomnia signal. Varying the **illustrative** HFRS effect-prior SD from 0.05 to 1.50 gives H4≈0.0505–0.0508 at p12=10⁻⁵ and ≈0.00529–0.00532 at 10⁻⁶. With HFRS prior SD 0.15, varying insomnia's prior SD to 0.1 or 0.3 gives H4≈0.0610 or 0.0438 at p12=10⁻⁵. These choices are sensitivity analyses, not proof that any prior is correctly calibrated for the raw HFRS scale. The minimum matched p-values were 2.99×10⁻¹⁰ for HFRS and 2.33×10⁻⁶ for registry insomnia.

Conditional **only if H4 is true**, rs7783012 has the largest SNP contribution, 0.131, under the 0.15/0.2 effect-prior setting; the conditional 95% set contains 22 variants. The previously highlighted rs1476535 has conditional contribution 0.0126. These values are not unconditional probabilities of causality. The analysis does not assign FOXP2 or any other gene as the effector, validate AlphaGenome, establish a sleep–ADHD–frailty signal, or resolve the failed LAVA gate. The earlier [FI–registry-insomnia coloc](../RESULTS.md) uses a different frailty definition and source; comparing its H4 directly would conflate phenotype, source and model choices.

## Input checks and limitations

The complete 280,558,723-byte UKB HFRS gzip passed its publisher MD5 and was streamed without retaining the full archive. Of 2,080 regional rows, 1,374 passed the specified UKB filters and **1,212 (88.2%)** had a unique matching, eligible R13 row after liftover and frequency alignment. The exclusions were 706 UKB ineligible, 37 without a unique R13 match and 125 R13 ineligible; no unmatched variant is a negative result. The maximum HFRS source-p versus BETA/SE normal-p discrepancy in the matched set was 0.000099 on the log10 scale. The underlying source allele/build conventions are supported by the match, not fully authenticated.

The authors' UKB file reports N=407,643, while its study description reports 407,463; the reason is unresolved. Association statistics may share participants or relatives, and cross-cohort independence is not verified. The British and Finnish sources have different ascertainment and LD; one causal variant per trait is a strong assumption for this region. The posterior is conditional on the chosen variants, priors and model, and is not adjusted for the many prior project analyses. It cannot rule out secondary sharing or a sleep-independent frailty effect.

The [pre-result protocol](PROTOCOL.md), [matched rows](../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-chr7-matched-20261028.tsv`), [all 14 posterior fits](../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-chr7-posteriors-20261028.tsv`), [conditional SNP table](../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-chr7-conditional-20261028.tsv`), and [source and matching receipts](../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-chr7-match-receipt-20261028.json`) preserve inputs and results. The SSD also contains the extraction, matching, coloc and verification scripts and their receipts. A separate Python implementation reproduced all 14 posterior vectors to within 2.4×10⁻¹⁵ and checked the conditional probabilities sum to one; this verifies arithmetic, not the model's biological validity.
