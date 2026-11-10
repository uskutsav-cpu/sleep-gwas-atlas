# FOXP2-region candidate across four All of Us EUR coded phenotypes

The [public All by All](https://allbyall.researchallofus.org/) EUR-group genome results for rs12536335 use GRCh38 `chr7:114403104 A>G` in each of four previously retrieved candidate lookups. The following values are **marginal variant associations**. Under the [source-code-supported but not independently authenticated ALT-effect convention](stage6_hunt_aou_allele_direction_review.md), negative beta means the G allele is associated with the lower-coded binary outcome relative to A. The 22-marker equal-weight prototype counts A; the separate 15- and 22-marker weighted prototypes use their own allele and weight conventions.

| Coded outcome | API beta for G, conditional | Nominal p |
|---|---:|---:|
| Insomnia (`NS_333.2`) | −0.0513 | 0.00000553 |
| Difficulty walking (`NS_325.2`) | −0.0667 | 0.00361 |
| Malaise/fatigue (`SS_807`) | −0.0170 | 0.0434 |
| Generalized weakness (`SS_807.21`) | +0.00562 | 0.777 |

The negative signs for the first three are consistent on the assumed G-effect allele, but this is a **postselected descriptive observation**. In the [33-variant insomnia–walking max-p screen](stage6_aou_insomnia_walking_eur_conjunction_review.md), this variant's BH q was 0.119, not significant at 0.05. Its fatigue max-p screen also did not pass the within-33 correction. The weakness association is imprecise and nonsignificant; this does not prove no weakness effect, nor does a sign difference prove different mechanisms. The [public PhecodeX vocabulary](stage6_aou_phecodex_mapping_overlap_review.md) maps weakness codes within the broad fatigue category, although actual All of Us case sets and local mapping are not established. The outcomes can share participants and diagnoses; they are not independent cohorts. Walking and weakness codes are not measured FFS, and the broad fatigue phenotype overlaps with sleep and exhaustion.

The [All by All walkthrough](https://allbyall.researchallofus.org/walkthrough) cautions that displayed beta magnitude is not readily interpretable given association-testing weights. Do not compare these raw beta magnitudes as causal effect sizes or calculate clinical risks from them. The allele contrast, source coding, participant and relative overlap, different phenotype definitions and potential LD tagging require further study. Postselected candidate results cannot identify a causal variant or assign causality to nearby FOXP2. Separate [regional colocalization](stage6_foxp2_independent_cohort_coloc_review.md) remains reference- and prior-sensitive and does not establish a causal gene. None of these findings validates an individual genetic score.

The [four-row extract](stage6_foxp2_aou_eur_direction_context.csv) and [receipt](stage6_foxp2_aou_eur_direction_context_receipt.json) reproduce with `python3 work/aou_complete/foxp2_eur_direction_context.py` using the existing query files; exact coordinate, allele pair, ancestry and data type are asserted. No participant data were used.
