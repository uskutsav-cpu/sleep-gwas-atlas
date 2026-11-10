# Public nonagenarian blood-expression follow-up

The public [Vitality 90+ expression series GSE40366](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40366) includes a 151-person nonagenarian subset with reported frailty categories and PBMC expression. Its metadata identify 29 frail, 38 non-frail and 84 pre-frail samples; we compared frail with non-frail and excluded pre-frail. This is observational bulk blood expression, not a genetic association, clinical insomnia result, cell-type-specific result or score validation. Vitality 90+ also appears as a cohort in the project's earlier longevity GWAS follow-up, so these analyses should not be counted as wholly independent studies without participant and family overlap checks.

Four candidate genes were chosen because of prior project hypotheses, not a prespecified cohort analysis. For genes with multiple platform probes we selected the highest mean-expression probe across all 151 category-labeled subjects, without using group labels to select probes. On the 67 frail/non-frail subjects we regressed the series' submitted log₂ quantile-normalized values on group and sex, using HC3 standard errors. Exact individual ages are not available in these metadata beyond nonagenarian status; no adjustment for age, cell composition, smoking or other health factors was possible. BH adjustment covers only these four postselected genes.

| Gene | Sex-adjusted frail minus non-frail difference | HC3 SE | Nominal p | BH q within four genes |
|---|---:|---:|---:|---:|
| FOXP2 | −0.059 | 0.051 | 0.249 | 0.714 |
| SMAD5 | −0.052 | 0.084 | 0.540 | 0.714 |
| RBM6 | −0.074 | 0.091 | 0.419 | 0.714 |
| MEIS1 | +0.074 | 0.200 | 0.714 | 0.714 |

We did not detect differences in this small, selected cross-sectional sample. These estimates do not exclude effects in other populations, cell types or tissues and do not resolve the baseline-study confounding in the earlier muscle dataset. Differences in phenotype definitions, age, cohort selection, platform and tissue preclude treating the blood and muscle results as direct replications or pooling their coefficients. Selection to age 90 introduces survival selection, and genetic confounding and causal direction were not assessed.

The original GEO SuperSeries includes 327 expression samples from two subseries processed using different normalization; the 151 frailty-labeled samples are all from the log₂/quantile-normalized subset, rather than the separate variance-stabilized subset. The count of 151 is not 151 frail cases. Reproduce with `work/geo_vitality/fetch_candidate_platform.py` and `work/geo_vitality/analyze_candidates.py`; probe mapping and aggregate results are in `outputs/stage6_gse40366_candidate_probe_map.json` and `outputs/stage6_gse40366_frailty_expression.json`. The downloaded matrix SHA-256 is `4e8af4f5e45cd5b29046edaec72eda15150dc3b15b54b1ef61d1795e0441b68d`.
