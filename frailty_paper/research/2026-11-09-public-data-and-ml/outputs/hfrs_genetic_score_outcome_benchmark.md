# Published HFRS genetic score: outcome benchmark

[Mak et al. (2025)](https://www.nature.com/articles/s43587-025-00925-y) trained a genome-wide genetic score using the FinnGen hospital frailty risk score (HFRS) GWAS and tested it in UK Biobank. Their score uses approximately 1.1 million variants and continuous-shrinkage weights. It is not the project's 15- or 22-variant prototype, and its results cannot be attributed to those candidates.

For each standard-deviation increase in the full HFRS genetic score, their Supplementary Table 15 reports:

| UK Biobank outcome | Reported association | Approximate 95% CI |
|---|---:|---:|
| Standardized HFRS | β = +0.074 | +0.071 to +0.077 |
| HFRS >5 before age 65 | OR = 1.250 | 1.236 to 1.265 |
| All-cause mortality | HR = 1.083 | 1.072 to 1.094 |

The intervals are Wald calculations from the paper's beta and SE, not independently refitted estimates. HFRS is a weighted hospital-diagnosis measure, and the before-65 outcome is defined by HFRS >5, not by the project's FI or FFS. The OR and HR are relative associations per score SD in the reported models, not absolute risk or the effect of a treatment. The study also reports associations with hospitalizations. For the HFRS score without dementia, the mortality HR is 1.076 (approximate 95% CI 1.065–1.087); removing dementia does not remove all potential sleep or exhaustion overlap.

The paper also compares models in Supplementary Table 16. For all-cause mortality, adding the HFRS score to birth year, sex and ten ancestry PCs changed the reported **C-index from 0.648 to 0.650**, a difference of 0.002. For hospitalization count, the reported R² changed from **0.03659 to 0.04219**, a difference of 0.00560 (0.56 percentage points of outcome variance). These are model-level, within-study comparisons, not proof of a clinically meaningful improvement, calibration, prospective benefit, or performance in a new population. We have not established whether the C-index and R² were evaluated in a separately held-out subset of UK Biobank. A very small association p-value does not by itself establish useful prediction.

## Source caveats

Supplementary Table 15 says models include birth year, sex, smoking and ten PCs, whereas the article's PRS methods paragraph lists birth year, sex and PCs without smoking. Table 16 labels its baseline and score models by birth year, sex and PCs. These are the source descriptions; exact covariates should be checked before reusing results. Table 15 stores several very small p-values as numerical zero while separately reporting −log10(p), so zero must not be interpreted as a literal p-value of zero. We avoided copying the HFRS p-value from a PDF text extraction that dropped an exponent digit; the supplementary −log10(p) is about 521.7.

FinnGen was the HFRS-score discovery cohort and UK Biobank the reported evaluation cohort, but UK Biobank contributes to the project's FI, FFS and UKB insomnia discovery results. Thus this study cannot validate the project's selected candidate sets, and UK Biobank is not automatically a clean test for them. Finnish and UK Biobank recruitment, HFRS coding and ascertainment differ; source sample overlap and per-analysis sample sizes should be checked. This genetic HFRS score also differs from the **clinical HFRS**: associations of measured HFRS with mortality must not be substituted for associations of the polygenic score with mortality.

## Traceability

The full-sample rows in `hfrs_score_published_outcomes.csv` are transcribed from the [published supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43587-025-00925-y/MediaObjects/43587_2025_925_MOESM1_ESM.xlsx), ST15; OR/HR intervals were calculated using exp(beta ± 1.95996×SE). The incremental model statistics are from ST16 rows 5–6 and 11–12. Workbook SHA-256: `123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625`. Output CSV SHA-256: `a2d741be20d16cbc343f86beb5791fa1e5914dfa6acd298ebd91b8dd322aca8d`. No individual data were accessed, no published score was applied, and no prediction model was refitted.
