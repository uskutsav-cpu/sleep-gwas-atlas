# Regional context for absent Foote factor records

For each selected marker absent from a factor's harmonized indexed release, we queried a ±250 kb GRCh38 interval around the lifted coordinate. The table counts markers whose interval returned zero source rows; multiple markers can occupy the same region and must not be interpreted as independent gaps.

| Factor | Selected markers absent | Those with zero returned rows within ±250 kb |
|---|---:|---:|
| General frailty | 55 | 50 |
| Limited social support | 11 | 0 |
| Unhealthy lifestyle | 11 | 0 |
| Multimorbidity | 18 | 8 |
| Metabolic problems | 48 | 38 |
| Poorer cognition | 11 | 0 |
| Disability | 23 | 4 |

In the general-factor file, all selected missing markers on chr7 (8), chr11 (9), chr16 (4), chr17 (14), and chr22 (15) had empty ±250 kb windows. By contrast, the four missing chr5 markers and one missing chr18 marker had 937–1,249 and 976 returned neighboring records, respectively. The missing chr5 group includes rs60271, rs10073929, rs40465 and rs4703307; an isolated missing target can arise for several reasons, and the indexed query does not establish which one applies.

The [published Foote methods](https://www.nature.com/articles/s41588-025-02269-0) state that, for downstream analyses, Q-heterogeneity-significant variants and ±1 Mb neighborhoods were pruned. Such pruning could produce regional gaps, but input coverage, harmonization, an indexed-release issue or other steps are also possible; the actual cause of each gap is **not established**. The publisher metadata for the submitted file reports removal of significant Q values and does not by itself prove the exact set of flanking windows in this released file. Neither an empty interval nor a missing marker means no association, trait specificity, or a null genetic effect.

Individual query counts and nearest returned distances are in the per-factor `*_absence_context.tsv` files; receipts record queried intervals and input hashes. The query is descriptive, uses remote indexed data without full-file checksum verification, and concerns highly selected markers. See [protocol](ABSENCE_PROTOCOL.md) and the [original crosswalk](EXTENDED_RESULTS.md).

A [partial original-file check](ORIGINAL_MISSING_RESULTS.md) did not observe the five isolated targets in sampled ranges, but the original is not guaranteed sorted; their presence or absence in the complete source is still unresolved.
