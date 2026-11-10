# Older-adult insomnia PBMC expression: four-candidate exploratory comparison

The public [GSE208668 series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208668) provides submitted log₂-transformed, quantile-normalized, gene-symbol-level PBMC expression for 17 older adults with insomnia disorder and 25 without. It does **not** provide genotype data or a frailty outcome and cannot validate the project scores or establish variant-to-gene mechanisms. The [original study](https://www.sciencedirect.com/science/article/am/pii/S0889159122003440) says all samples were assayed in a single batch with the same protocol; nonetheless, the insomnia and control samples were recruited from different studies, so insomnia status and recruitment source are perfectly confounded. Differences in enrollment and unmeasured study or participant characteristics remain inseparable from the insomnia contrast. GEO also says raw data were lost; this limits independent reassessment of normalization and probe-level QC.

Using four previously selected project genes and **no outcome-driven filtering**, we fitted each submitted value on insomnia group, age and sex, with HC3 robust standard errors and residual-degrees-of-freedom t intervals. Multiple-testing adjustment is Benjamini–Hochberg across these four genes only, not across the many project analyses. The following are insomnia-group minus control-group differences on the submitted log₂ expression scale; these are exploratory and cannot resolve study-source confounding.

| Gene | Difference | HC3 SE | Nominal p | BH q among four |
|---|---:|---:|---:|---:|
| FOXP2 | −0.156 | 0.077 | 0.0488 | 0.0488 |
| SMAD5 | +0.475 | 0.070 | 4.36×10⁻⁸ | 1.74×10⁻⁷ |
| RBM6 | +0.457 | 0.138 | 0.00207 | 0.00276 |
| MEIS1 | +0.640 | 0.166 | 0.000444 | 0.000888 |

Small p-values do not distinguish insomnia-associated expression from recruitment-source effects, establish robustness across analyses, or support causal interpretation. No detection p-values, probe-level provenance, cell composition, genetic covariates or independent replication were checked. FOXP2's modest difference is near the stated four-gene threshold and its observed values cluster near 7 on this normalized scale; without probe-level data and background detection we cannot establish reliable expression of this gene in these PBMCs. The selection is post hoc in the overall project.

The previous [GSE40366 frailty PBMC comparison](stage6_gse40366_frailty_expression_review.md) found no BH-significant differences for these four genes in 29 frail versus 38 non-frail nonagenarians. That result neither replicates nor contradicts this insomnia analysis: cohorts, age, disease, assay and modeling differ; absence of significance is not evidence of no effect, and neither study has both phenotypes and genetic scores. Participant overlap across sources has not been checked. A single-cell or QTL study would need to connect measured variant, cell type, gene and phenotype before a mechanism is inferred.

Reproduce using `work/stage4_env/bin/python work/geo_insomnia/analyze_candidates.py`. [Aggregate estimates](stage6_gse208668_candidate_expression.csv) and the [source/checksum receipt](stage6_gse208668_candidate_expression.json) are provided; no individual expression rows are exported. Source matrix SHA-256 is recorded in the receipt. The analysis uses the submitted symbol-level matrix, not individual probe intensities; the separate platform file is not required.
