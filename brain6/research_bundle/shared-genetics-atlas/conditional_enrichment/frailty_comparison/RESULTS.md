# Conditional insomnia patterns across FI, FFS and HFRS

On the same 1,137,927 allele- and chromosome-matched HapMap3 markers, 2.33% have approximate insomnia p≤0.01 overall. Among markers selected separately by each frailty phenotype at p≤0.001, the fraction is:

| Conditioning phenotype | Selected SNPs | Also insomnia p≤0.01 | Fraction |
|---|---:|---:|---:|
| FI | 9,024 | 553 | 6.13% |
| FFS | 12,132 | 684 | 5.64% |
| HFRS | 16,963 | 1,459 | 8.60% |

Another 335 of the original 1,138,262 common-support SNPs were held by the joint join and are not treated as null. The groups of selected SNPs differ; these percentages are **not genetic correlations or estimates of causal effects**. Registry insomnia and HFRS use overlapping FinnGen releases, while FI and FFS overlap UK Biobank. Power, ascertainment, phenotype content, ancestry and source methods also differ. HFRS uses source-reported p-values, whereas the other traits' p-values are normal approximations from rounded Z. These factors can change the apparent ordering.

Paired physical-block jackknife comparisons gave FI−FFS = +0.49 percentage points, FI−HFRS = −2.47 points, and FFS−HFRS = −2.96 points. Their descriptive 1.96-SE intervals using 5 Mb blocks were respectively [−1.86, +2.84], [−4.96, +0.01] and [−5.42, −0.51] points; the 10 Mb results were similar, except the FI−HFRS interval's upper end was just below zero. These are **not calibrated significance tests**: the thresholds were chosen after inspecting related analyses, physical blocks may not capture long-range LD, and no adjustment for the several comparisons or source differences was made. A nominal interval excluding zero cannot be promoted as a finding of biological difference or a causal frailty-specific pathway.

The [protocol](PROTOCOL.md), `run.R`, `summary.tsv` and `block_counts.tsv` specify the shared support, test thresholds and paired jackknife. `verify.py` independently recomputes all six differences and SEs from block counts, with source and output hashes in `verification.json`. This analysis does not perform enrichment-based conjFDR, identify independent loci, remove sleep or exhaustion items, or replicate any two-trait association independently.

A [regional sensitivity](regional_sensitivity/RESULTS.md) shows the HFRS fraction changes little when the broad chr19 APOE-containing window is excluded, whereas excluding all chromosome 5 lowers the FI- and FFS-conditioned fractions. These post-hoc SNP counts are not independent loci or causal evidence.
