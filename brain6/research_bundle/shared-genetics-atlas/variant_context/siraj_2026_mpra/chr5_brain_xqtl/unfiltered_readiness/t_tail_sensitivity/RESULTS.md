# Sensitivity to the methylation p-value tail

The [preceding pilot](../RESULTS.md) treated the source's beta/SE values as Gaussian Z statistics, although the reported methylation p-values closely match a t distribution with an empirically inferred 370 degrees of freedom. To assess sensitivity to that difference, we constructed a **heuristic** comparison: each source two-sided p-value was converted to an equal-tail standard-normal score, preserving the sign of beta, while the original standard error and effect-size prior were kept unchanged. This substitution is **not a coherent t-based likelihood or a calibrated Bayes factor**; it is only a check of how the earlier score responds to weaker evidence from heavier tails.

| Pair | Original H4 range over tested priors | Heuristic p-matched H4 range |
|---|---:|---:|
| CpG–insomnia | 0.000075–0.00160 | 0.000085–0.00252 |
| CpG–ADHD | 0.0000020–0.000418 | 0.0000055–0.000911 |

Under both calculations, H3 remains the largest component across the tested settings. The result is not proof that causal variants differ, that secondary sharing is absent, or that these posterior values are calibrated. It also does not resolve source overlap, causal-variant multiplicity, missing GWAS frequencies, region selection, phenotype scale or the failed LAVA quality gates. A properly specified likelihood, model-matched data and independent validation are needed for stronger inference.

[All comparisons](comparison.tsv), [code](compare.py), and [receipt](receipt.json) preserve the inputs, priors and output hashes. Independently pairing the 16 original rows with the earlier R output reproduced all posterior components within 4.5×10⁻¹⁶; this is an arithmetic check only.
