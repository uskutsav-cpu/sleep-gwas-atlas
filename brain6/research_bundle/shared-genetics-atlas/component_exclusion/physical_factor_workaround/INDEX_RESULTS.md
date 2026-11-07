# Explicit genetic physical-function index

This post-hoc workaround defines an equal-weight index of the **standardized additive genetic components** of walking pace, left grip strength and moderate-activity days. Each source measure has no explicit sleep or exhaustion item. For physical genetic-correlation matrix R, the index is `(Gwalking + Ggrip + Gactivity) / sqrt(1'R1)`. Its correlation with insomnia is the sum of the three insomnia–indicator genetic correlations divided by the same normalization. Here the pre-normalization genetic variance is 4.017.

The definition fixes the weights and does not impose the latent factor's zero indicator-residual covariance assumption. It does not make the index biologically independent of sleep, a measured clinical frailty score or an individual's polygenic score. Equal weighting applies to standardized genetic components, not raw walking categories, kilograms or activity days.

| Insomnia source | Index genetic correlation | Paired block-jackknife SE | Unadjusted normal 95% interval | p adjusted for these two source estimates |
|---|---:|---:|---:|---:|
| Registry diagnosis | −0.0736 | 0.0350 | −0.1423 to −0.0050 | 0.0711 |
| Questionnaire | −0.2309 | 0.0260 | −0.2818 to −0.1799 | 1.34×10⁻¹⁸ |

The registry estimate is nominally negative (unadjusted p≈0.0355) but does not meet the two-source Bonferroni benchmark. The questionnaire estimate is more negative; questionnaire and physical data share UK Biobank, and the source comparison is not independent replication. The paired questionnaire-minus-registry difference is −0.1572 (SE 0.0415; unadjusted normal p≈0.000152). Phenotype, population, ascertainment, sample overlap and reference fit differ together, so this does not identify the reason for that difference. Neither correction covers the broader exploratory project or the post-hoc choice to analyze this index.

The calculations use the existing 1,053,614-marker common support and 200 aligned block deletes. Numerator and normalization are recomputed in every delete. All inspected point and delete correlation matrices were positive semidefinite, and the index correlations were within bounds. An independent R calculation reproduced both estimates and jackknife SEs to within 10⁻¹². The [protocol](INDEX_PROTOCOL.md), [Python calculation](explicit_index.py), [numerical results and input hashes](joint_fit/explicit_index.json), and earlier [source verification](VERIFICATION.md) document the analysis. Original GWAS files were not rehashed in this run because the SSD was disconnected; the local-summary reconstruction verified recorded log hashes and common per-trait delete vectors.

This gives a reproducible, explicit genetic physical-function sensitivity without a latent zero-residual assumption. It remains a distinct construct from Fried frailty without exhaustion or FI without sleep/fatigue. The original LD-score assumptions, uncertain source overlap, Finnish–general-European LD fit, rounded summary values, phenotype-selection effects and lack of independent validation remain. No new genome-wide index association file, individual-risk predictor, causal estimate or matched item-deletion result was produced, and these inspected outcomes were not added to ML training.

A [weight-definition sensitivity and scientific figure](WEIGHT_RESULTS.md) show that emphasizing activity shifts estimates toward zero, while the equal-weight index remains unchanged. The optimized point-estimate ranges are descriptive and should not be used to select a stronger association.

The [local reproduction workflow](REPRODUCE.md) successfully regenerated the saved-summary analyses and figures without the external SSD. This reproduces the calculations, not the original GWAS/LDSC pipeline or independent validation.

After the SSD reconnected on 6 November, the [retained-input LDSC replay](VERIFICATION.md) also reproduced all nine genetic correlations and SEs and all 27 block-delete vectors exactly. All ten retained input hashes matched their archived receipts. This verifies the LDSC stage with the same inputs and reference; it does not reproduce raw GWAS generation or provide independent validation.
