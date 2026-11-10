# Follow-up across HUNT, Partners and All of Us: six variants

**Follow-up:** A [source-code and allele-direction audit](stage6_hunt_aou_allele_direction_review.md) supports the working ALT-effect assumption for the EUR results; the comparison remains conditional and postselected.

29 October 2026. The [published HUNT and Partners insomnia tables](stage6_lane2019_candidate_replication_review.md) cover six project variants. The same six have public All of Us All by All EUR-group associations for coded insomnia and generalized weakness. The [joined results](stage6_lane_aou_six_variant_comparison.csv) retain each source p-value and allele representation, without pooling effect sizes.

| Variant | HUNT insomnia p | Partners insomnia p | All of Us EUR insomnia p | All of Us EUR weakness p |
|---|---:|---:|---:|---:|
| rs113851554 | 0.000182 | 0.0002891 | 4.21 × 10⁻⁷ | 0.674 |
| rs4886140 | 0.0119 | 0.508 | 0.000408 | 0.912 |
| rs1942262 | 0.00414 | 0.081 | 0.0761 | 0.783 |
| rs9845387 | 0.0108 | 0.468 | 0.280 | 0.287 |
| rs35881094 | 0.177 | 0.108 | 0.664 | 0.470 |
| rs10865954 | 0.462 | 0.532 | 0.125 | 0.605 |

The two variants with unadjusted p < 0.05 in both HUNT and All of Us insomnia results are rs113851554 and rs4886140. rs113851554 also has a nominal Partners association, but the All of Us variant lookup returned two different allele records at that position; the extracted record is G/T, matching the Lane table allele pair, while the other is G/A. See the [allele check](stage6_meis1_aou_multiallelic_lookup_review.md). Exact representation and effect coding still need further verification. The HUNT/Partners odds ratios for the listed source allele appear in the CSV, but effect direction was not harmonized or meta-analyzed against All of Us. The source positions are treated as GRCh37 for Lane’s published tables and GRCh38 for All of Us. An [Ensembl assembly-mapping audit](stage6_lane_aou_six_liftover_review.md) found that all six coordinates map uniquely on the forward strand to the All of Us positions, with matching allele pairs. The Lane source-build assumption, exact multiallelic contrasts, sample overlaps and effect orientation still need checking before a directional replication claim.

None of these six has unadjusted p < 0.05 for the All of Us weakness code. This is consistent with the previously reported [33-variant shared-association screen](stage6_aou_insomnia_weakness_conjunction_review.md); it does not demonstrate the absence of a weakness or frailty effect. Weakness is an EHR-coded phenotype, not a Fried frailty score or deficit-accumulation FI; the samples and ascertainment differ between sources. These postselected findings are not a prespecified replication test, proof of shared causality, or validation of an individual-level score. The source studies may overlap with discovery or candidate selection, and independence has not been established.

Sources: [Lane et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6415688/) and [All of Us All by All documentation](https://allbyall.researchallofus.org/about); the local [retrieval review](stage6_aou_full_candidate_association_review.md) explains the queried codes. Run `python3 work/cross_cohort/compare_hunt_aou.py` to regenerate the CSV and its [checksum receipt](stage6_lane_aou_six_variant_receipt.json).
