# Candidate-variant FI and physical-frailty results side by side

The [22-marker comparison](shared22_mak_fi_fp_effects.csv) presents aligned, per-allele effects, standard errors, approximate 95% confidence intervals and p-values from Mak et al.'s public UK Biobank FI and modified Fried physical-phenotype (FP) GWAS. Both outcomes were reported for 403,041 people, but equal sample counts do not prove identical participants for each trait or SNP. The FP is a 0–5 Fried-style count and may not be identical to the Ye FFS. The FI and FP effect coefficients use different outcome scales; **do not compare their raw magnitudes or subtract them as a test of stronger association**.

For illustration, effects below are aligned to the score's counted allele and kept on their respective source scales. Intervals are approximate normal intervals, not absolute disease-risk estimates.

| Candidate; counted allele | FI beta (95% CI); p | Physical-phenotype beta (95% CI); p |
|---|---|---|
| rs12536335 A, FOXP2-labeled | 0.1292 (0.0964 to 0.1620); 1.21 × 10⁻¹⁴ | 0.00561 (0.00193 to 0.00928); 0.00279 |
| rs28528997 A, MAD1L1-labeled | 0.1896 (0.1483 to 0.2308); 2.09 × 10⁻¹⁹ | 0.01730 (0.01268 to 0.02192); 2.12 × 10⁻¹³ |
| rs2155292 G, NCAM1-labeled | 0.1163 (0.0834 to 0.1492); 4.56 × 10⁻¹² | 0.00157 (−0.00212 to 0.00526); 0.403 |
| rs6452792 A, CTC-498M16.4-labeled | 0.03532 (0.00260 to 0.06804); 0.0344 | −0.00207 (−0.00573 to 0.00160); 0.269 |

For rs2155292, the FP interval includes zero and small effects in either direction. For rs6452792 the negative FP point estimate is also imprecise; its sign is not evidence of a protective physical-frailty effect. A result crossing a threshold for one outcome and not the other is **not a formal test that the effects differ**. Such a comparison needs a meaningful common outcome scale and the covariance of estimated effects, which may be influenced by sample overlap. These intervals are individual unadjusted intervals and do not account for having selected or examined many markers.

Across the 22 selected candidates, all 22 Mak FI estimates are positive and nominal p < 0.05; 11 have nominal p < 0.05 for both Mak outcomes. This descriptive pattern is strongly affected by selecting candidates from earlier insomnia–FI results, reuse or potential reuse of UK Biobank, correlated markers, and a nominal threshold: it is not 22 independent confirmations or evidence that half the variants are FI-specific. Several candidate rsIDs have multiallelic representations and require source-level verification.

## Reproducibility

The comparison is extracted from [the eight-source aligned matrix](shared22_eight_source_evidence_matrix.csv), itself based on the [Mak public GWAS release](https://doi.org/10.6084/m9.figshare.24591573.v1). The extraction verified all 88 reported confidence interval bounds against the exact normal 97.5th-percentile multiplier and the published beta and SE values. Alleles are aligned to the unweighted shared-direction count; nearby-gene names are positional labels, not proof of causality. Candidate selection, LD and genotype-encoding limits are detailed in the [mentor table](shared22_mentor_risk_table.md). These selected-variant results complement, but do not replace, the [genome-wide FI–FP genetic-correlation analysis](mak2023_same_cohort_fi_fp_ldsc_review.md). Neither validates a genetic risk score.
