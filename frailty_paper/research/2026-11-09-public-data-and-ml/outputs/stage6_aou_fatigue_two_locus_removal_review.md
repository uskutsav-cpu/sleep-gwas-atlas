# Fatigue weighted-Z sensitivity after removing insomnia-associated loci

29 October 2026. The [public All of Us EUR malaise/fatigue weighted-Z check](stage6_aou_eur_fatigue_weighted_z_review.md) gave positive conditional summary statistics for all seven exploratory score constructions. We repeated the calculation after removing rs113851554 and rs12536335, which were selected in the prior insomnia follow-up. The 15-marker prototype contains both; the 22-marker sets contain only rs12536335. Removal changes the score and its denominator. It is not a comparison of predictive performance, and it is not evidence that deleting a variant improves a score.

| Construction | Full-set conditional Z | Z after removing the available selected loci |
|---|---:|---:|
| 15-marker FI | 3.217 | 2.835 (13 markers) |
| 15-marker FFS | 3.182 | 2.645 (13 markers) |
| 15-marker insomnia | 1.963 | 2.376 (13 markers) |
| 22-marker equal-weight | 4.659 | 4.441 (21 markers) |
| 22-marker FI | 4.713 | 4.539 (21 markers) |
| 22-marker FFS | 3.412 | 3.110 (21 markers) |
| 22-marker insomnia | 4.573 | 4.359 (21 markers) |

All seven shortened constructions retain positive conditional Z values. This shows that the conditional signal does not disappear when those terms are removed. Other retained variants may tag the same underlying signals through linkage disequilibrium; this does not establish independent effects, a genuine score–fatigue association in people, or biological specificity. The two loci, candidate lists and outcomes were chosen post hoc. Their omission changes both numerator and reference-LD variance, so an increase in Z cannot establish an improved model or that the omitted variant has a protective effect. The full insomnia-weighted 15-marker Z was sensitive to deletion of rs113851554, highlighting that the particular weighting and marker set matter.

Coded fatigue can overlap directly with insomnia and the exhaustion component of frailty. It is not a frailty index or FFS. The source ALT-effect assumption, source-sample independence, exact multiallelic contrasts, imputed candidate quality, HWE variance and representative reference LD remain unverified as described in the original review. FFS-weight direction is source-supported but inferred. No multiple-testing-adjusted or prespecified replication conclusion follows.

The [full sensitivity table](stage6_aou_fatigue_two_locus_removal.csv) includes each single-locus removal, and the [receipt](stage6_aou_fatigue_two_locus_removal_receipt.json) records source and output hashes. Run `work/stage4_env/bin/python work/aou_complete/fatigue_locus_removal_sensitivity.py` to reproduce it. It reconstructs each full-set statistic from the published component terms and checks it against the prior file before calculating removals. No individual-level data were used.
