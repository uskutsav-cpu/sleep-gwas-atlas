# PIP-weighted AlphaGenome accessibility: exploratory sensitivity

This post-hoc analysis combines previously saved AlphaGenome/JBrowse adult frontal-cortex DNase predictions with conditional single-signal fine-mapping PIPs, for the six previously selected insomnia–ADHD regions. It is descriptive, restricted to the selected variant union and eligible trait-specific PIPs, and subject to the approved non-commercial terms. No model was trained on AlphaGenome outputs. See [protocol](PROTOCOL.md), [code](analyze.py), [full region summary](region_summary.tsv), [variant contributions](variant_contributions.tsv), and [input hashes](receipt.json). Input hashes match the earlier ranking analysis.

In the chromosome 3 region, rs62264764 contributes **40.6%** of the insomnia-weighted and **38.2%** of the ADHD-weighted sums. Its conditional PIPs are **0.0402** and **0.0312**, respectively, and its predicted ALT–REF DNase summary is **−0.2734 log₂**. Its influence comes from the combination of these PIPs and a relatively large predicted absolute change; the contribution shares are not the probability that this variant causes either phenotype. The chromosome 3 weighted mean is 0.0273 for insomnia and 0.0235 for ADHD, versus 0.0202 unweighted on the eligible selected variants.

| Region (GRCh37 bounds) | Trait | PIP available / selected | Selected PIP mass | Unweighted mean absolute DNase | PIP-weighted mean absolute DNase | Largest contributor | Share |
|---|---|---:|---:|---:|---:|---|---:|
| chr11:28076261–29076261 | adhd | 88/88 | 0.954 | 0.0170 | 0.0177 | rs11030386 | 8.4% |
| chr14:98137479–99137479 | adhd | 69/69 | 0.951 | 0.0140 | 0.0193 | rs74469887 | 22.5% |
| chr3:117142005–118142005 | insomnia | 128/128 | 0.991 | 0.0202 | 0.0273 | rs62264764 | 40.6% |
| chr3:117142005–118142005 | adhd | 128/128 | 0.950 | 0.0202 | 0.0235 | rs62264764 | 38.2% |
| chr5:103447968–104447968 | insomnia | 71/71 | 0.951 | 0.0119 | 0.0119 | rs161645 | 19.7% |
| chr5:103447968–104447968 | adhd | 70/71 | 0.990 | 0.0120 | 0.0079 | rs77960 | 41.6% |
| chr5:87185500–88185500 | adhd | 3/3 | 0.973 | 0.0090 | 0.0120 | rs4916723 | 94.1% |
| chr7:113571035–114571035 | insomnia | 66/66 | 0.965 | 0.0117 | 0.0116 | rs1476535 | 13.5% |
| chr7:113571035–114571035 | adhd | 66/66 | 0.968 | 0.0117 | 0.0181 | rs9969232 | 43.5% |

One ADHD PIP is unavailable in the chromosome 5:103447968–104447968 selected union; it is excluded rather than set to zero, and the unweighted comparison uses the same 70 eligible variants. Three region–trait pairs did not pass their association gate and are omitted. Selected mass is the sum of available PIPs in the selected union, not an unconditional probability of causality. Contribution shares were checked to sum to one.

These values combine separate models and cannot establish a functional effect, a disease-risk direction, a causal gene, tissue enrichment, or independent validation. PIPs are conditional on the assumed single-signal model and its priors; LD, selected variants, the post-hoc choice of track and settings, missingness and source overlap constrain interpretation. Scores are predictions rather than observed molecular effects. The [setting sensitivity](../alphagenome_window_sensitivity/RESULTS.md) and [rank comparison](../prediction_rank_comparison/REPORT.md) help place this result in context.

## Sensitivity to removing the largest contributor

A post-hoc leave-one-out calculation removes the largest PIP × absolute-prediction contributor in each eligible region–trait pair, then renormalizes the remaining PIPs. For chromosome 3 this removes rs62264764 for both traits: the weighted mean falls from 0.0273 to **0.0169** for insomnia and from 0.0235 to **0.0150** for ADHD. The unweighted mean of those same 127 remaining variants is **0.0182** for both traits. Therefore the original comparison above the unweighted mean depends heavily on this one prediction; it should not be described as evidence of a general regional enrichment. This is a sensitivity calculation, not evidence that the variant is wrong or right, and the removal was chosen after observing the data.

The full nine-pair [table](leave_one_out.tsv) and [reproducible script](leave_one_out.py) are provided. For the chromosome 5:87 Mb ADHD set, only three variants were selected and removal of its top contributor leaves PIP mass 0.0897, so that particular remaining mean is especially sensitive to selection.
