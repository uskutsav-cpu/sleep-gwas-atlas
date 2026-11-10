# All of Us AFR-group fatigue: exploratory score and variant sensitivity

29 October 2026. The public All by All AFR genetic-similarity group has 8,680 coded malaise/fatigue cases and 45,674 controls in its analysis metadata; the corresponding EUR group reports 39,924 and 115,544. We queried all 33 distinct candidates for `SS_807` using the browser's whole-genome slow Hail route. Each exact GRCh38 position/REF/ALT was returned. These are summary associations, not individual-level data; the group label is not equivalent to self-identified race and does not imply a homogeneous population.

Using the same conditional weighted marginal-Z method as the earlier [EUR fatigue sensitivity](stage6_aou_eur_fatigue_weighted_z_review.md), with 661 1000 Genomes AFR-labeled reference individuals, results were:

| Construction | Conditional Z | Nominal p |
|---|---:|---:|
| 15-marker FI | 0.707 | 0.479 |
| 15-marker FFS | 0.637 | 0.524 |
| 15-marker insomnia | 0.635 | 0.526 |
| 22-marker count | 2.101 | 0.0357 |
| 22-marker FI | 2.118 | 0.0342 |
| 22-marker FFS | 2.004 | 0.0450 |
| 22-marker insomnia | 1.925 | 0.0543 |

The three nominal p-values below 0.05 are unadjusted; none passes BH q < 0.05 when these seven correlated exploratory constructions alone are treated as a family. That limited correction cannot account for the wider project search. These scores, outcome and follow-up were selected post hoc; the conditional statistics are not observed patient score effects, predictive performance, or clinical validation.

We also joined these 33 associations to the previously retrieved AFR insomnia results using exact coordinate and alleles. No candidate has p < 0.05 for both; none passes BH on the 33 max-p values, and the smallest max-p is 0.102. This max-p intersection–union screen is a distinct question from the weighted summary; neither proves that effects are absent or that causal variants differ by trait. As in EUR, fatigue is a broad clinical symptom and can overlap with insomnia or exhaustion. No direction or colocalization test was performed in the max-p screen.

A [public vocabulary audit](stage6_aou_phecodex_mapping_overlap_review.md) shows that standard PhecodeX v1.0 and v1.1 fatigue mappings include the two codes for generalized weakness. The exact deployed mapping and patient overlap remain unverified, so fatigue and weakness are not independent replication outcomes.

Do not compare significance thresholds across EUR and AFR as evidence of differing effects or score transferability: the groups have different sizes, ascertainment, allele frequencies, LD and sampling uncertainty, and no formal between-group test was done. The 1KG AFR reference is not the All of Us target sample. HWE genotype variance, ALT effect coding, exact biallelic contrasts, common effective N, sample/relative independence and source allele quality remain assumptions; multiallelic rsIDs and the FFS-weight convention require further checks. `SE × HWE dosage SD` has CV 1.42% (15) and 1.20% (22), a limited diagnostic rather than proof of calibration.

Reproduce weighted results with `work/stage4_env/bin/python work/aou_complete/afr_fatigue_weighted_z_sensitivity.py`, and the conjunction with `python3 work/aou_complete/insomnia_fatigue_afr_conjunction.py`. Outputs: [seven score sensitivities](stage6_aou_afr_fatigue_weighted_z_sensitivity.csv), [per-marker terms](stage6_aou_afr_fatigue_weighted_z_components.csv), [variant conjunctions](stage6_aou_insomnia_fatigue_afr_conjunction.csv), and their [score](stage6_aou_afr_fatigue_weighted_z_sensitivity_receipt.json), [conjunction](stage6_aou_insomnia_fatigue_afr_conjunction_receipt.json), and [query](stage6_aou_afr_fatigue_candidate_association_receipt.json) receipts. No individual score or risk was computed.
