# Comparing two FI–physical-frailty GWAS pairs across insomnia sources

A second frailty study provides a useful **construct sensitivity**. [Mak et al. 2023](https://doi.org/10.1111/acel.13868) released UK Biobank GWAS for a 49-item FI and a modified Fried physical phenotype (FP), each described as 403,041 participants. Equal reported sample sizes do not establish identical participant sets. The FP contains exhaustion, and the FI retains sleep-related items; neither is an item-exclusion analysis. The original project pair is the Atkins FI meta-analysis and Ye FFS.

For each insomnia source, four prior native LDSC fits used a common marker set. This independent arithmetic audit reproduced all 12 stored jackknife SEs and verified identical ordered SNP hashes, shared-insomnia h² deletion vectors and block boundaries within each source. Values below are r_g (SE):

| Insomnia source | Atkins FI | Ye FFS | Mak FI | Mak FP | Common SNPs |
|---|---:|---:|---:|---:|---:|
| FinnGen R13 registry | 0.4292 (0.0393) | 0.2438 (0.0360) | 0.4009 (0.0284) | 0.2370 (0.0353) | 1,129,146 |
| UK Biobank symptoms | 0.6409 (0.0236) | 0.4841 (0.0269) | 0.6412 (0.0193) | 0.4933 (0.0261) | 1,098,887 |
| MVP clinical insomnia | 0.5516 (0.0393) | 0.3965 (0.0359) | 0.5356 (0.0288) | 0.4043 (0.0335) | 1,124,391 |

Both physical-phenotype estimates are lower than their corresponding FI estimates in these comparisons. The paired differences, approximate unadjusted 95% CIs and conservative Bonferroni p-values for the **six** within-pair comparisons are:

| Insomnia source | Mak FP − FI (paired SE) | Bonferroni p, 6 tests | Ye FFS − Atkins FI (paired SE) | Bonferroni p, 6 tests |
|---|---:|---:|---:|---:|
| R13 | −0.1640 (0.0264) | 3.00×10⁻⁹ | −0.1854 (0.0399) | 2.04×10⁻⁵ |
| UKB | −0.1479 (0.0206) | 3.94×10⁻¹² | −0.1569 (0.0281) | 1.45×10⁻⁷ |
| MVP | −0.1313 (0.0262) | 3.23×10⁻⁶ | −0.1551 (0.0355) | 7.69×10⁻⁵ |

These corrections are for a **post-hoc selected family**, not all analyses across the atlas. The three source-specific tests of the difference between the two FI-versus-physical contrasts do not reject zero after a separate three-test Bonferroni correction; that does not prove the contrasts are equal. Source-specific supports differ and a direct comparison across sleep sources was not performed. The paired jackknife accounts for covariance in these recorded block estimates under the model; it does not remove systematic bias.

The frailty studies reuse or may reuse UK Biobank participants, and the UKB insomnia source also overlaps that population. FinnGen and MVP source intersections and relatedness have not been independently established. FinnGen uses a general-European rather than Finnish-specific LD reference and a total N; the prior MVP run had an elevated univariate intercept (~1.114), while its source uses p-derived Z and rounded reported odds ratios. Mak GWAS coding, the assumed GRCh37 position mapping and actual variant-level N require their source documentation. These findings do not show sleep-free physical frailty, causal effects, clinical-score correlations or independent two-trait replication.

The [protocol](PROTOCOL.md), `audited_contrasts.tsv`, `audit_receipt.json`, 12 copied source exact records and support receipts preserve the comparisons. Full source reports and diagnostics were previously generated in the user's separate `/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/outputs` project. See [source review](PRIOR_SOURCE_REVIEW.md), noting that its paths are from that source project, and the [item-exclusion protocol](../../component_exclusion/CONSTRUCTION_PROTOCOL.md) for what is still missing.
