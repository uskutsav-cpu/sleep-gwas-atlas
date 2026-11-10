# MEIS1-region rs113851554: checking the two All of Us records

On 29 October 2026, the public All by All EUR-group insomnia interval query at GRCh38 chromosome 2:66,523,432 returned two rows, not two identical copies of the selected variant. The rows were:

| All by All variant ID | Allele frequency reported | Insomnia p-value |
|---|---:|---:|
| chr2-66523432-G-T | 0.054609 | 4.209443 × 10⁻⁷ |
| chr2-66523432-G-A | 0.00000904072 | 0.5949036 |

The existing project query selected the **G/T row** by exact variant ID; its `returned_count=2` reports the number of rows in the queried interval, not duplicate G/T matches. The Lane HUNT table reports the same base pair as T/G, with T as the published effect allele and T frequency 0.05612 in HUNT. This frequency is broadly similar to the All by All G/T row's reported frequency; it is only a consistency check, not proof of harmonization or sample independence. The 15-variant prototype counts T. Matching alleles do not resolve genome-build conversion, strand, multiallelic dosage treatment, sample overlap, or the All by All beta's effect-allele semantics. G/A must not be substituted for G/T.

The [All by All documentation](https://allbyall.researchallofus.org/walkthrough) says the single-variant beta indicates direction but that its magnitude is not readily interpretable given the association-testing weights. The walkthrough does not, on the cited page, define the allele to which that beta refers. Do not interpret the HUNT OR and All of Us beta as independently harmonized effects, or as predictive-score weights, on the basis of this lookup. The risk of restless-legs confounding and the evidence limitations in the [six-variant comparison](stage6_lane_aou_six_variant_comparison.md) remain.

The public query is [the All by All interval API](https://allbyall.researchallofus.org/api/variants/associations/interval/chr2:66523431-66523433?analysis_id=NS_333.2&ancestry_group=eur&sequencing_type=genome&query_mode=slow). Query settings were phenotype NS_333.2, EUR, genome and slow interval. The first comparison's exact-variant filtering is implemented in `work/aou_complete/query_candidates_eur.py`; HUNT values come from Table S7 of the [Lane et al. supplement](https://pmc-oa-opendata.s3.amazonaws.com/PMC6415688.1/NIHMS1519714-supplement-2.xlsx). Neither these results nor this allele check validate the full score.
