# MVP insomnia and two pairs of frailty measures

Using the European-ancestry Million Veteran Program (MVP) clinical-insomnia GWAS, the genetic correlation is again higher for FI than for the physical-frailty measure in both study pairs. All four LDSC fits below use the same **1,124,391** allele-aligned SNPs:

| Frailty GWAS | MVP insomnia genetic correlation r_g (SE) |
|---|---:|
| Mak FI | 0.536 (0.029) |
| Mak modified Fried physical phenotype (FP) | 0.404 (0.034) |
| Atkins FI meta-analysis | 0.552 (0.039) |
| Ye Fried frailty score (FFS) | 0.396 (0.036) |

For the Mak pair, FP minus FI is **−0.131** (paired SE 0.026; approximate 95% CI **−0.183 to −0.080**). For the original pair, FFS minus FI is **−0.155** (SE 0.036; CI −0.225 to −0.085). The difference between those two within-pair contrasts is +0.024 (SE 0.024; CI −0.024 to +0.072; unadjusted p = 0.331). This does not establish equality. The calculations are exploratory and unadjusted for the broader sequence of analyses.

The earlier [matched FinnGen-insomnia comparison](mak2023_insomnia_fi_fp_comparison.md) also found higher correlations with FI for both frailty pairs. The insomnia sources differ in diagnosis, samples and SNP support, and no paired test was performed across insomnia sources; these results do not establish that the MVP and FinnGen effects differ or agree exactly. The two frailty pairs reuse or may reuse UK Biobank participants, while exact overlap within the Mak pair is unverified. These are not independent replications. FI includes sleep-related items and the physical phenotype includes exhaustion. Item removal has not been tested, and genetics, phenotype definitions, selection, ascertainment and statistical assumptions may all matter. Genetic correlation is not individual risk or causality.

## Diagnostics and traceability

All four native fits used the same SNP and allele order, same European LD reference, and identical 200 block boundaries; the MVP heritability and all deletion estimates matched across fits. Each individual correlation SE was independently reproduced from its delete-block estimates; paired contrasts use the same deleted blocks. MVP's observed-scale h² was 0.0291 (SE 0.0021), with univariate LDSC intercept **1.114**. Mak FI and FP intercepts were 1.064 and 1.041. Heritability precision exceeds the project's ≥4 screening gate, but that does not establish absence of bias; the elevated MVP intercept warrants model and source-QC review. The MVP summary data use per-variant sample sizes where available; Z was constructed in the existing project workflow using p-values and effect direction because standard errors were not supplied. Rounded odds ratios can make directions ambiguous for effects near one, and the existing preprocessing and filters matter. Sample overlap and relatedness across cohorts are unknown; European-ancestry group labels do not guarantee perfect LD-reference fit.

The common set intersects the four frailty sources with the project's filtered MVP insomnia file and was already restricted by earlier HFRS and FinnGen-insomnia support. This is a sensitivity on selected SNP support, not a new GWAS. Source and allele details for the Mak frailty files are in [the Mak comparison](mak2023_same_cohort_fi_fp_ldsc_review.md); prior MVP preparation is described in the project’s [source comparison](insomnia_frailty_risk_briefing.md). [Input hashes](mak_mvp_matched_support.json), [calculation receipt](mak_mvp_matched_contrasts.json) and the four `mak_mvp_*.exact.json` files and logs preserve results. Preparation and summary scripts: `work/mak2023_mvp_common.py` and `work/mak_mvp_contrasts.py`. No individual-level outcome data were used, and the sparse score remains unvalidated.
