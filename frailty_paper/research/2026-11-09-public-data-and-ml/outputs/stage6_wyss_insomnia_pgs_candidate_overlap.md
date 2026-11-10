# Candidate overlap with published 2026 insomnia scores

Checked 29 October 2026. We streamed the complete PGS Catalog scoring files for [PGS019928](https://www.pgscatalog.org/score/PGS019928/) (UK Biobank-derived, 1,069,748 rows) and [PGS019929](https://www.pgscatalog.org/score/PGS019929/) (MVP European-ancestry-derived, 1,085,967 rows). Both are genome-wide insomnia scores evaluated for sleep phenotypes by [Wyss et al.](https://doi.org/10.1093/sleep/zsag200), not the project's postselected scores. The project has used the MVP insomnia GWAS in follow-up, and UK Biobank contributes to discovery; these score files are not independent confirmations of our candidate selection.

| Published score | Candidate file | Rows identified | Derived GRCh38 coordinate and allele pair match |
|---|---|---:|---:|
| PGS019928 | 15-marker | 5 of 15 | 5 |
| PGS019928 | 22-marker | 7 of 22 | 7 |
| PGS019929 | 15-marker | 5 of 15 | 5 |
| PGS019929 | 22-marker | 7 of 22 | 7 |

These are candidate-level file matches, not performance results. PGS019928 provides positions without rsIDs, so matching used project-derived GRCh38 coordinates and then checked the two-allele pair. PGS019929 matching checked both its supplied rsID and the project-derived coordinate, followed by allele-pair checks. The project coordinate file for the 22 candidates contains positions, but Ensembl verification was incomplete for some entries; a coordinate/allele-pair match is not complete verification of a multiallelic source contrast. A second complete scan on 29 October also queried position regardless of rsID and found no additional matching rows; a differently named same-position row would now be flagged for allele review. This still does not establish that candidates are biologically absent because build mapping, multiallelic normalization and proxy variants need separate consideration. It also does not examine proxy markers. Directly verified matching rows are in [the extract](stage6_wyss_pgs_candidate_overlap.csv); [the streaming script](../work/compare_public_insomnia_pgs.py) records matching and counts.

The two published scores contain roughly one million variants, use their own LD-aware weights and were evaluated as complete scores; association or performance of those scores **cannot be assigned to these few overlapping variants** or transferred to our 15-/22-marker sums. Different signs or magnitudes of conditional and marginal weights would not by themselves establish agreement, disagreement, causality or validation. This is an explanatory benchmark and provenance check, not an individual-level score test or an independent evidence of shared insomnia–frailty risk.
