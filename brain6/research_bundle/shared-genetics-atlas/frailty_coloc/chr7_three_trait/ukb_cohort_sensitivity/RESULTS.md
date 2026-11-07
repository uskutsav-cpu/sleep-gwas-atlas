# Chromosome 7: UK Biobank-only FI source sensitivity

On the same 1,055 variants as the [parent pilot](../RESULTS.md), replacing the FI meta-analysis with its UK Biobank cohort component left the default HyPrColoc cluster as **FI–insomnia**, with posterior 0.8794 versus 0.8822 and candidate rs12671330 in both. No three-trait cluster was returned. At stricter thresholds neither source returned a cluster. The posterior is conditional on the single-signal model, priors, selected region and an **unverified independence assumption**; this is not an independently replicated signal.

| Model and setting | Meta-analysis FI | UKB-only FI |
|---|---:|---:|
| HyPrColoc FI–insomnia posterior, conditional prior 0.02 | 0.8822 | 0.8794 |
| Same, conditional prior 0.01 | 0.7870 | 0.7824 |
| Same, conditional prior 0.005 | 0.6428 | 0.6363 |
| Pairwise FI–insomnia H4, p12=10⁻⁵ / 10⁻⁶ | 0.9528 / 0.6688 | 0.9516 / 0.6628 |
| Pairwise FI–ADHD H4, p12=10⁻⁵ / 10⁻⁶ | 0.0441 / 0.00460 | 0.3689 / 0.0552 |
| Pairwise FI–ADHD H3, p12=10⁻⁵ | 0.9557 | 0.6310 |

The FI–ADHD pairwise posterior is sensitive to the FI source and sharing prior. Under UKB-only FI, neither H3 nor H4 dominates decisively at p12=10⁻⁵; H4 drops at p12=10⁻⁶. The three-trait and pairwise models answer different conditional questions, so their posteriors cannot be multiplied or directly equated. The source change combines sample composition, precision and potentially other source-processing differences; it cannot establish the cause of the shift. UKB is part of the meta-analysis, and both analyses reuse the insomnia and ADHD sources.

The UKB-only source avoids interpreting the TwinGene beta for this particular fit; it does **not** resolve the [source orientation ambiguity](../../../placo_diagnostic/completed_audit_fi_only/prior_literature/cohort_orientation/RESULTS.md) or prove the meta-analysis is biased. UKB FI includes insomnia and tiredness, and UKB FI and the insomnia source share UK Biobank, so the assumed independence is particularly uncertain. The previously examined [study-correlation stress test](../overlap_sensitivity/RESULTS.md) shows that plausible illustrative dependence can change the cluster. These fits do not rescue failed LAVA QC, estimate local genetic correlation, test component-excluded FI, or demonstrate a causal gene or absence of secondary signals.

The [protocol](PROTOCOL.md), [matching receipt](matching_receipt.json), [cluster results](clusters.tsv), [pairwise posteriors](pairwise.tsv), [code](run.R), and [independent arithmetic verification](verification.json) document this exploratory comparison. The complete 16,557,800-row UKB/TwinGene cohort source was scanned with gzip integrity checking; all 1,055 previously matched SNPs had unique, position- and allele-consistent valid UKB estimates. The matched file stays on the SSD because it contains restricted ADHD-derived inputs. The parent meta clusters reproduce exactly, and an independent calculation matches all eight pairwise posterior vectors within 2.7×10⁻¹⁵. These checks do not validate the statistical assumptions or biological interpretation. [Atkins cohort data](https://figshare.com/articles/dataset/Genome-Wide_Association_Study_of_the_Frailty_Index_-_Atkins_et_al_2021/9204998).

A [variant-weight decomposition](variant_diagnostic/RESULTS.md) shows that the FI–ADHD posterior shift is driven by small model-conditional overlap outside the individual 95% candidate sets. This does not determine why the source inputs differ or rule out secondary shared signals.
