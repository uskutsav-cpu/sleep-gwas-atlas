# Insomnia source sensitivity across frailty and physical function

On the same 1,020,637 markers, the original Brain6 questionnaire-insomnia and FinnGen registry-insomnia GWAS show different patterns of genome-wide genetic correlation with three frailty measures and two physical-function traits. This is a post-hoc analysis of existing data, not independent replication or an item-deletion experiment. Negative physical-function correlations assume the original source coding; the raw phenotype transformations were not independently reconstructed.

| Outcome | Questionnaire insomnia r_g (SE) | Registry insomnia r_g (SE) | Questionnaire minus registry (paired SE) | Contrast p, Bonferroni across five |
| --- | ---: | ---: | ---: | ---: |
| FI | 0.6443 (0.0244) | 0.4301 (0.0372) | 0.2142 (0.0414) | 1.18×10⁻⁶ |
| FFS | 0.4869 (0.0284) | 0.2435 (0.0358) | 0.2434 (0.0397) | 4.48×10⁻⁹ |
| HFRS | 0.4363 (0.0253) | 0.6131 (0.0326) | −0.1768 (0.0397) | 4.15×10⁻⁵ |
| Walking pace | −0.2900 (0.0264) | −0.0608 (0.0319) | −0.2292 (0.0374) | 4.62×10⁻⁹ |
| Left grip | −0.1680 (0.0251) | −0.0947 (0.0345) | −0.0733 (0.0328) | 0.1270 |

FI, FFS, HFRS and walking-pace contrasts pass this bounded five-contrast correction; the grip contrast does not. Point differences use the rounded LDSC estimates and the p-values use paired 200-block jackknife SEs. These contrasts cannot tell us why the sources differ: phenotype definitions, cohort makeup, participant overlap, ancestry, ascertainment and source models change together. FI contains insomnia and tiredness, and FFS includes exhaustion; HFRS is based on hospital diagnoses. Questionnaire insomnia shares UK Biobank data with FI, FFS and physical measures; registry insomnia and HFRS both use FinnGen, with exact overlap unknown. The HFRS intercept is elevated in previous audits, and Finnish reference fit and effective N remain concerns. None of these r_g values is an individual-level phenotypic correlation or a causal effect.

For a separate post-hoc family of **ten individual correlations**, all FI, FFS and HFRS values and the two questionnaire physical-function correlations meet Bonferroni 0.05. Registry walking and registry grip do not: their corrected p-values are 0.566 and 0.060, respectively. The same registry-grip result met the earlier [two-test pilot](../../RESULTS.md); the difference in status reflects the larger correction family here, as well as slightly different common support. Neither family is a correction for all atlas-wide analyses.

## Verification and source data

All three added frailty files matched previously recorded hashes and had ordered SNP and allele keys identical to the archived questionnaire-insomnia file. Their exact subset was matched to the earlier four-source common support. Every fit retained the same 1,020,637 valid-allele SNPs. For every outcome the two source fits have identical outcome-heritability delete vectors; for each insomnia source the five fits have identical insomnia-heritability delete vectors. The reconstructed jackknife SEs match the logs to printed precision. These checks support the paired arithmetic, not the original GWAS models or ancestry assumptions. Per-marker N semantics, rounded source Z values, and the prior source warnings remain as documented in the [physical-function sensitivity](../RESULTS.md) and [frailty source audits](../../../../global_audit/insomnia_frailty_same_support/RESULTS.md).

The [pre-run protocol](PROTOCOL.md), [matching code](prepare.py), [matching receipt](prepare_receipt.json), [LDSC script](run.sh), [calculation](summarize.py), [full numerical output](results.json), and per-pair logs and delete vectors document this analysis. Matched input files remain on the SSD. The earlier three-frailty comparisons used different marker support and correction families, so their exact estimates are not expected to match this table. The matched sleep/exhaustion-excluded GWAS and a valid genome-wide local-correlation family remain outstanding.
