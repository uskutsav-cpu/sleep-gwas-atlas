# UK Biobank insomnia and the two FI–physical-frailty pairs

Using the project's UK Biobank insomnia-symptom GWAS and a common set of 1,098,887 aligned SNPs, the genetic correlation is higher with FI than with physical frailty for both GWAS pairs:

| Frailty GWAS | UKB insomnia genetic correlation r_g (SE) |
|---|---:|
| Mak FI | 0.641 (0.019) |
| Mak modified Fried physical phenotype (FP) | 0.493 (0.026) |
| Atkins FI meta-analysis | 0.641 (0.024) |
| Ye Fried frailty score (FFS) | 0.484 (0.027) |

For Mak FP minus FI, the paired difference is **−0.148** (SE 0.021; approximate 95% CI **−0.188 to −0.108**). For original FFS minus FI, the difference is **−0.157** (SE 0.028; CI −0.212 to −0.102). The difference between these two contrasts is +0.009 (SE 0.018; CI −0.027 to +0.045; unadjusted p = 0.624), which does not establish equality. These comparisons were made after multiple earlier analyses and p-values are unadjusted.

The same pattern is observed using the project's [FinnGen](mak2023_insomnia_fi_fp_comparison.md) and [MVP](mak2023_mvp_insomnia_fi_fp_comparison.md) insomnia sources. Those runs used different SNP supports and participants; no formal between-source comparison was made. The original and Mak frailty sources reuse or may reuse UK Biobank participants, as does the UKB insomnia study; they are not independent replications. FI includes sleeplessness and tiredness, while the physical phenotype includes exhaustion. The shared source population and overlapping measures make it particularly important not to interpret these correlations as causal, sleep-free frailty associations, individual risk probabilities or clinical score correlations.

## Allele exclusions and diagnostics

Of the **1,129,146** variants in the prior four-frailty common set, **30,259** have blank A1, A2, Z and N fields in the existing munged UKB insomnia file. The local munge log and full-file scan show that these are consistent with reference-list rows without association statistics, rather than evidence of missing raw-source alleles. They were excluded; the remaining 1,098,887 matched on SNP identifier and directly, swapped or strand-complemented nonpalindromic allele pairs. [Padding audit](ukb_insomnia_padding_audit.md). The munged insomnia file does not contain positions, so this is not an independent insomnia-source position check. The individual reasons for unavailable association results have not been established; this check does not demonstrate a problem in the original raw GWAS. The original previous project LDSC automatically required valid alleles, but this analysis does not reconstruct or repair that source file. [Allele-QC counts](mak_ukb_insomnia_allele_qc.json) and [matched-file hashes](mak_ukb_matched_support.json) record the exclusions.

The four LDSC fits use the same ordered SNP hash, LD reference and 200 matching block boundaries; insomnia h² and its delete-block estimates agree across fits. The individual SEs and paired contrasts were recomputed from the block estimates. UKB insomnia h²/SE is ~22; this is a precision check, not proof that the model is unbiased. Its cross-trait intercept with Mak FI is 0.280 and with Mak FP 0.145; shared participants and correlated errors may contribute, but these values are not clinical correlations or counts of overlapping people. The intercepts were freely estimated. The common-support comparison was additionally restricted by prior HFRS and FinnGen insomnia availability; its observations should not be treated as independent GWAS or a valid individual score.

[Exact numerical results and contrasts](mak_ukb_matched_contrasts.json) and the `mak_ukb_*.exact.json` and `.log` files preserve the results. The [Mak FI–FP analysis](mak2023_same_cohort_fi_fp_ldsc_review.md) documents frailty source acquisition and allele checks. Scripts: `work/mak2023_ukb_common.py`, `work/ukb_insomnia_allele_qc.py`, and `work/mak_ukb_contrasts.py`. No participant-level outcome or prediction test was conducted.
