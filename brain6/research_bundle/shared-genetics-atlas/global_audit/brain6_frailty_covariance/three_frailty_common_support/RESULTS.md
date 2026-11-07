# Original Brain6 insomnia compared with FI, FFS and HFRS

On 1,099,360 matched SNPs, all three frailty measures show positive genetic correlations with the original Brain6 insomnia GWAS. FI has the largest estimate in this dataset; the FFS–HFRS difference is not statistically distinguishable from zero in the post-hoc three-contrast family.

| Frailty phenotype | r_g (SE) | p | Bonferroni p, 3 correlations |
|---|---:|---:|---:|
| FI | 0.6404 (0.0236) | 1.65×10⁻¹⁶² | 4.95×10⁻¹⁶² |
| FFS | 0.4840 (0.0269) | 4.07×10⁻⁷² | 1.22×10⁻⁷¹ |
| HFRS | 0.4355 (0.0255) | 1.75×10⁻⁶⁵ | 5.26×10⁻⁶⁵ |

| Post-hoc contrast | Difference | Paired SE | Unadjusted 95% CI | Bonferroni p, 3 contrasts |
|---|---:|---:|---:|---:|
| HFRS − FI | −0.2049 | 0.0328 | −0.2692 to −0.1406 | 1.23×10⁻⁹ |
| HFRS − FFS | −0.0485 | 0.0335 | −0.1141 to 0.0171 | 0.442 |
| FI − FFS | 0.1564 | 0.0281 | 0.1013 to 0.2115 | 7.81×10⁻⁸ |

The contrasts use paired 200-block jackknife SEs and the four-decimal LDSC point estimates, so their statistics are approximate. The p-value corrections cover only the three post-hoc correlations and, separately, the three post-hoc contrasts; they do not cover all analyses in the atlas. Intervals are unadjusted for multiple comparisons. A nonsignificant FFS–HFRS difference does not establish equivalence.

These are same-data analyses of different phenotypes and source cohorts. FI includes insomnia and tiredness; FFS includes exhaustion; HFRS derives from weighted hospital diagnoses. Differences cannot be assigned to physical frailty, causal pathways, or one particular cohort feature. FI and FFS have UK Biobank participation; HFRS is from FinnGen, and exact overlaps and relatedness are not established. The HFRS univariate LDSC intercept is **1.1875 (SE 0.0139)**; ancestry/reference mismatch, ascertainment and the use of study-level total N rather than confirmed per-variant N remain concerns. The original insomnia–HFRS cross-trait intercept is approximately zero (0.000004, SE 0.0062); this does not prove no participant overlap or model bias.

## Checks and context

The HFRS input was taken from the previously audited HFRS support and matched to the frozen FI/FFS/insomnia subset. All four ordered SNP/allele keys and file hashes were checked. Each LDSC fit retained all 1,099,360 markers and the shared insomnia h² delete vectors were identical. Reconstructed per-trait SEs agree with LDSC's displayed SEs, and [an independent R check](check_paired.R) reproduced all three paired SEs and corrected p-values. These computational checks do not validate study ancestry, sample overlap, per-SNP N or causal interpretation.

The inputs are a filtered common support selected for previous analyses, not all genomic variants; the HFRS coordinate audit identified 537 unresolved markers on the larger prior support and its restricted-support sensitivity is documented [here](../../insomnia_frailty_same_support/coordinate_sensitivity/RESULTS.md). A separate [FinnGen registry-insomnia comparison](../../insomnia_frailty_same_support/RESULTS.md) used different insomnia data and obtained a different ranking; the differences in phenotype, support and cohorts prevent attributing that change to a single cause. No sleep-item-excluded GWAS, independent replication, local genetic correlation, or rescue of failed LAVA QC follows from this comparison.

The [protocol](PROTOCOL.md), [HFRS matching code](prepare_hfrs.py), [receipt](matching_receipt.json), [run script](run.sh), [HFRS log](HFRS.log), [correlation table](correlations.tsv), [contrast table](contrasts.tsv), [delete values](delete_rg.tsv), [paired code](paired.py), and [paired receipt](paired_receipt.json) document the analysis. The earlier FI/FFS inputs and logs are preserved in the [two-frailty comparison](../paired_common_support/RESULTS.md); source-derived matched inputs remain on the SSD. Earlier receipts capture the code, reference and source hashes.

A [same-marker insomnia N sensitivity](insomnia_n_sensitivity/RESULTS.md) leaves the three r_g estimates and paired contrasts nearly unchanged at displayed precision, while neither N convention is validated by that agreement.
