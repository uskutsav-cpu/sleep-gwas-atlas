# All of Us insomnia and malaise/fatigue: exploratory shared-association screen

29 October 2026. We queried public [All by All](https://allbyall.researchallofus.org/about) EUR genetic-similarity group whole-genome results for `SS_807`, described as **malaise and fatigue**, for all 33 distinct project candidate variants. The public analysis metadata reported 39,924 cases and 115,544 controls in this group. We joined the exact GRCh38 position/REF/ALT pairs to the previously retrieved EUR coded-insomnia (`NS_333.2`) associations. All 33 rows were returned and matched the insomnia records. This is an observational association follow-up, not individual data or score validation.

For each variant, we computed `max(p_insomnia, p_fatigue)`, a two-sided intersection–union test of the null that at least one of the two trait associations is null. Under valid source p-values this maximum is valid without assuming independent phenotypes. We applied Benjamini–Hochberg (BH) across these 33 candidate max-p values; its standard FDR guarantee depends on the joint dependence structure. Five variants have both source p-values < 0.05, but **none has BH q < 0.05 within this screen**. The smallest max-p is 0.0114 (rs11790060), with q = 0.188. These are postselected variants and outcomes within a larger series of project searches, so the 33-test adjustment does not account for the broader search.

| Variant | Insomnia p | Malaise/fatigue p | Max p | BH q within 33 |
|---|---:|---:|---:|---:|
| rs11790060 | 0.0114 | 0.000653 | 0.0114 | 0.188 |
| rs4886140 | 0.000408 | 0.0114 | 0.0114 | 0.188 |
| rs28528997 | 0.0213 | 0.0205 | 0.0213 | 0.234 |
| rs1032524 | 0.00438 | 0.0310 | 0.0310 | 0.256 |
| rs12536335 | 0.00000553 | 0.0434 | 0.0434 | 0.287 |

A [public vocabulary audit](stage6_aou_phecodex_mapping_overlap_review.md) shows that standard PhecodeX v1.0 and v1.1 fatigue mappings include the two codes for generalized weakness. The exact deployed mapping and patient overlap remain unverified, so fatigue and weakness are not independent replication outcomes.

Malaise/fatigue is a broad clinical symptom code. It is relevant to the exhaustion component of some physical-frailty definitions, but is **not** the FFS, an FI, or a joint sleep–frailty outcome; fatigue can also directly accompany insomnia, so phenotype overlap is a major alternative explanation. A variant with associations to both coded outcomes does not establish a shared causal variant, equal direction, a causal gene or a mechanism. Here, direction was not tested. Association tests can include overlapping samples and conditions; their true joint case counts and cohort/relative overlap with discovery or selection data are unknown. The EUR group is a genetic-similarity category and does not represent all ancestry groups. No lack-of-association conclusion follows from a nonsignificant p-value.

The full [33-variant table](stage6_aou_insomnia_fatigue_eur_conjunction.csv), [query results](stage6_aou_eur_fatigue_candidate_association.csv), and [checksummed calculation receipt](stage6_aou_insomnia_fatigue_eur_conjunction_receipt.json) are available. Reproduce the join and BH calculation with `python3 work/aou_complete/insomnia_fatigue_conjunction.py`; query via `python3 work/aou_complete/query_candidates_eur_fatigue.py` (cached rows are reused). Source details are in the [query receipt](stage6_aou_eur_fatigue_candidate_association_receipt.json). The prior [insomnia–weakness screen](stage6_aou_insomnia_weakness_conjunction_review.md) found no nominal joint associations among these candidates in either EUR or AFR. Differences in results do not demonstrate outcome-specific biological effects, because fatigue and weakness differ in prevalence, definition and statistical power.
