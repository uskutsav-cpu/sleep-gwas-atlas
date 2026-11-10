# How fatigue and exhaustion can overlap in UK Biobank frailty measures

The physical-frailty phenotype in Mak et al. scores exhaustion using a question about tiredness or low energy during the past two weeks. That question matches the wording of UK Biobank [field 2080](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2080), which the published [Williams FI derivation script](https://raw.githubusercontent.com/dylwil/ukbfrailty2/main/Stata%20script%20for%20FI%20derivation%20in%20UKB%20data%20v3.do) uses for the FI fatigue item. This gives a concrete mechanism for overlap between the two measures, beyond their being generally related constructs.

**Limit of verification:** Mak's supplementary FI item list calls the item fatigue but does not identify its field or provide its recoding. The evidence supports a strong inference that the same question contributes to Mak's FI and physical phenotype, but it does not independently establish that the exact Williams script was used to construct Mak's FI. Verify Mak's original construction code before treating the example below as its actual coding. Ye's FI and FFS implementations likewise require source-specific verification.

The following table compares the Williams script with the physical-frailty exhaustion rule reported in Mak Supplementary Table 5, assuming that both use the same field 2080 response for a participant. Response labels are from [UK Biobank coding 100484](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100484).

| Response over the past two weeks | UKB code | FI fatigue item in Williams script | Mak physical-frailty exhaustion point |
|---|---:|---:|---:|
| Not at all | 1 | 0 | 0 |
| Several days | 2 | 0.25 | 0 |
| More than half the days | 3 | 0.75 | 1 |
| Nearly every day | 4 | 1 | 1 |

The measures use different transformations even in this matching-question example. The FI is a proportion of available deficits; the physical phenotype counts positive physical domains. Missing or special responses cannot be treated as zero merely because they do not match a positive category. A published FI table elsewhere reports a different intermediate value for response 3, so the values above are specifically those in the linked script, not a verified universal FI rule.

Separately, Mak's FI item list includes insomnia; UK Biobank [field 1200](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1200) asks about difficulty falling asleep or waking during the night. These direct or potential shared items complicate interpretation of an insomnia–FI genetic correlation and of FI–physical-frailty overlap. They do not tell us how much of either observed genetic correlation is attributable to the items. That requires matched analyses of reconstructed phenotypes or suitable summary statistics with the relevant items removed.

Mak's published modified FI removes insomnia, fatigue and six other mental-wellbeing items together; its full summary statistics were not identified in the checked release. It could support a broader sensitivity analysis if obtained, but would not isolate either item's contribution. See [the modified-FI lead](mak2023_sleep_excluded_fi_lead.md), [the FI–physical-frailty comparison](mak2023_same_cohort_fi_fp_ldsc_review.md), and [the physical-frailty crosswalk](mak_hanlon_physical_frailty_crosswalk.md).

Sources: Mak et al. [Supplementary Tables 4 and 5](https://doi.org/10.1111/acel.13868), workbook `ST4!A12:C20` and `ST5!A3:C5`; UK Biobank field and coding pages; Williams published FI derivation script. No participant-level records were analyzed for this comparison.
