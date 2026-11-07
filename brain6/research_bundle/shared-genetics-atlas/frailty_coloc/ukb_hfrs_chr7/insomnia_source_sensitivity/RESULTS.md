# Chromosome 7 HFRS colocalization: insomnia-source sensitivity

A subsequent [correlated-statistic stress test](overlap_stress/RESULTS.md) shows that the questionnaire HFRS–insomnia cluster can disappear under hypothetical positive correlation settings in a different method. The real residual correlation remains unknown.

On **1,125 identical variants**, the single-signal model gives substantially different results for questionnaire and registry insomnia when paired with the same UK Biobank hospital-frailty statistics. These are exploratory posterior probabilities conditional on the model and priors, not a formal test of the difference or proof of shared or distinct causal variants.

| Insomnia source | p12 | H1: HFRS only | H3: distinct signals | H4: shared signal |
| --- | ---: | ---: | ---: | ---: |
| Brain6 questionnaire | 10⁻⁵ | 0.0080 | 0.1204 | 0.8715 |
| Brain6 questionnaire | 10⁻⁶ | 0.0372 | 0.5586 | 0.4042 |
| FinnGen registry | 10⁻⁵ | 0.3172 | 0.6317 | 0.0510 |
| FinnGen registry | 10⁻⁶ | 0.3325 | 0.6621 | 0.00535 |

These rows use an illustrative HFRS raw-effect prior SD of 0.15, an insomnia log-odds prior SD of 0.2, and p1=p2=10⁻⁴. With p12=10⁻⁵, the questionnaire fit favors H4; with p12=10⁻⁶ it favors H3. The registry fit favors H3 under both. Varying HFRS prior SD from 0.05 to 1.50 barely changes the displayed H4 values, but the appropriate HFRS scale and prior remain unverified. Same marker support removes one difference in coverage; phenotype definition, ancestry, ascertainment, GWAS model, sample size and cohort dependence still differ together.

**The questionnaire insomnia and HFRS studies both use UK Biobank, yet these coloc calculations assume independent association statistics.** Exact participant overlap and the effect of that assumption are unknown; a shared-signal interpretation may be biased. The Finnish registry study uses another population and a different insomnia definition. Consequently the comparison cannot show whether the posterior difference is caused by cohort overlap, phenotype, ancestry or another source difference. It is not independent replication, a three-trait model, a gene assignment, a sleep-item-excluded frailty analysis, or a rescue of the failed LAVA gate. The one-signal-per-trait assumption can also miss secondary shared signals.

Conditional on H4 under the 0.15/0.2 setting, the questionnaire model's highest SNP contribution is rs2894699 (0.322), with rs12536335 (0.147) and rs1476535 (0.105) next; its conditional 95% set contains 18 SNPs. The registry model's top candidate is rs7783012 (0.131), and its conditional set has 22 SNPs. These are conditional contributions, not unconditional causal probabilities, and regional LD can make multiple markers tag a signal. FOXP2 remains a positional candidate rather than an established effector.

## Source checks

Of 1,212 previously matched HFRS–registry rows, 1,125 matched the frozen questionnaire source by unique rsID, position and allele pair; 87 lacked a unique questionnaire match. Its source lacks allele frequencies, so frequency compatibility with the questionnaire cohort could not be checked. The HFRS source's N discrepancy and other source limitations are retained in the [parent report](../RESULTS.md). An independent Python calculation reproduced all 20 posterior vectors to within 2.3×10⁻¹⁵ and checked both conditional SNP distributions; this verifies arithmetic, not statistical calibration or model assumptions.

The [pre-result protocol](PROTOCOL.md), [matched table](../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-source-chr7-matched-20261028.tsv`), [all posterior fits](../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-source-chr7-posteriors-20261028.tsv`), [conditional SNP table](../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-source-chr7-conditional-20261028.tsv`), [matching receipt](../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-source-chr7-match-receipt-20261028.json`) and [independent verification](../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-source-chr7-verify-20261028.json`) document the analysis. Source-matching, R model and verification scripts are on the SSD with the same filename prefix. Posterior probabilities are not adjusted for selection across the atlas.
