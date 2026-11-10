# Conditional direction comparison with the clinical-insomnia meta-analysis

For the previously extracted clinical-insomnia meta-analysis, we compared each project weight's sign with the association sign after expressing the source beta relative to the project's counted allele. This is a **postselected consistency check in overlapping cohorts**, not replication, score testing, genetic correlation or prediction. It assumes LocusZoom's beta refers to ALT, and that switching from ALT to REF gives a valid reversal for the same biallelic diploid contrast. Several project variants have multiallelic reference representations; source-specific coding and target representation have not been established. The FFS weights' original effect-allele convention also remains inferred.

| Research quantity | Conditional sign agreement | Markers |
|---|---:|---:|
| 15-marker FI weights | 14 | 15 |
| 15-marker FFS weights | 11 | 15 |
| 15-marker insomnia weights | 14 | 15 |
| 22-marker FI weights | 22 | 22 |
| 22-marker FFS weights | 21 | 22 |
| 22-marker insomnia weights | 22 | 22 |
| 22-marker equal-weight count | 22 | 22 |

These numbers are not independent confirmations or comparisons of predictive performance. The candidates were selected using related sleep/frailty results, and UK Biobank, FinnGen and MVP contribute to the source meta-analysis. High agreement is therefore unsurprising and could reflect selection and sample overlap. Different weight sets, scales, allele conventions, effect precision and LD make simple sign counts inadequate to rank scores. The insomnia outcome is based on diagnosis codes and differs from the symptom definitions of other studies. No uncertainty test or score effect is calculated; weak and strong effects count equally, and no adjustment for correlated variants is implied.

[Candidate-level conditional values](stage6_kunorozva_meta_direction_comparison.csv), [receipt](stage6_kunorozva_meta_direction_comparison_receipt.json) and [script](../work/kunorozva2025/compare_candidate_weight_directions.py) document the exact rows, weight files and allele transformations. Run `python3 work/kunorozva2025/compare_candidate_weight_directions.py`; it verifies all 33 unique source entries and allele-pair compatibility. The [source extraction](stage6_kunorozva_meta_candidate_extract_review.md) documents phenotype, source sample and remaining data quality limitations. None of the resulting signs validates either project score.
