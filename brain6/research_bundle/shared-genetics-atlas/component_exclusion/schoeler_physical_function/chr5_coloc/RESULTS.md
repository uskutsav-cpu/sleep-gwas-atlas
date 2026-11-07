# Chromosome 5: physical function versus insomnia and ADHD

A post-selected, single-signal colocalization sensitivity does **not establish a shared causal variant** between the public physical-function composite and either Brain6 trait. The physical-function association in the matched region is comparatively weak, and the results depend strongly on the prior probability of sharing.

| Comparison | Matched SNPs | H4 at p12=10⁻⁵ | H4 at p12=10⁻⁶ | H2 at p12=10⁻⁵ / 10⁻⁶ |
|---|---:|---:|---:|---:|
| Physical function–insomnia | 2,451 | 0.494 | 0.089 | 0.442 / 0.796 |
| Physical function–ADHD | 2,360 | 0.533 | 0.102 | 0.411 / 0.790 |

H4 is the model's posterior probability of one shared causal variant; H2 is the probability of association with the sleep/ADHD trait only, in the specified region and model. These probabilities are conditional on the priors, data and single-signal assumptions, not empirical replication rates. At the default physical-effect prior SD 0.15 and case-control SD 0.20, H3 (distinct signals) was 0.063 / 0.114 for insomnia and 0.056 / 0.107 for ADHD at the respective p12 values. The minimum matched physical-function p was 1.86×10⁻⁵; minimum sleep/ADHD p values reconstructed from beta/SE were 1.97×10⁻⁹ and 2.47×10⁻¹³. These are selected regional minima, not independent tests. The two match sets differ, so H4 values are not a formal cross-trait comparison.

The physical outcome averages standardized component scores; the composite's resulting standard deviation is not established as one by that statement. Consequently the default 0.15 effect-size prior is an assumption on the reported beta scale. A post-hoc [sensitivity](prior_sensitivity.tsv) using fixed physical prior SD 0.05–0.30 gives H4 0.343–0.701 for insomnia and 0.378–0.734 for ADHD at p12=10⁻⁵, and 0.050–0.190 and 0.057–0.216 respectively at p12=10⁻⁶. This further demonstrates dependence on model assumptions, not a means of selecting a favorable result.

The physical source is a lung-function/grip/height composite, not an exhaustion-removed frailty instrument. The original sleep and ADHD sources may overlap with UK Biobank, and the region was chosen after previous analyses. The frozen GWAS inputs are sparse relative to the full region and include earlier source assumptions, including case-control processing and rounded effects; multi-signal structure, cohort overlap and unresolved LAVA local QC can invalidate simple single-signal interpretations. No individual-variant, gene, clinical or AlphaGenome validation claim follows. The conditional SNP lists are conditional on H4 and **must not** be read as unconditional causal probabilities.

Of 3,615 physical rows in the indexed region, 2,796 passed allele/numeric filters. The source inputs contained 2,455 eligible insomnia and 2,360 eligible ADHD rows, with 2,451 and 2,360 uniquely matched by lifted coordinate and non-palindromic allele pair. See [protocol](PROTOCOL.md), [matching receipt](match_receipt.json), [matched data](insomnia_matched.tsv), [ADHD data](adhd_matched.tsv), [posterior table](posteriors.tsv), [prior sensitivity code](prior_sensitivity.R) and [verification](verification.json). An independent direct-ABF implementation reproduced H0–H4 with maximum absolute difference 1.67×10⁻¹⁵; this verifies arithmetic, not model validity. The GRCh38 source was accessed by partial indexed range query; its whole-file MD5 was not checked here. Source input and chain hashes were checked against existing receipts.

A subsequent [regional source-orientation audit](original_audit/RESULTS.md) found agreement between original and harmonized physical-function values after effect-allele alignment for every retained row in both matched sets. This does not validate the colocalization model or establish independent replication.

A post-hoc [common-marker sensitivity](common_support/RESULTS.md) restricted both comparisons to 2,359 identical positions. Results remained similarly prior-sensitive; this does not establish shared causality or an independent trait comparison.

A post-hoc [scale and numeric diagnostic](scaling_audit/RESULTS.md) found reported physical beta/SE and p-values internally consistent on the common markers. A simplified per-variant residual-scale proxy had median 0.475, but its assumptions prevent treating it as a measured outcome SD or using it to tune the prior.
