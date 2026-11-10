# Open muscle-expression follow-up for four candidate genes

The open [GSE144304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144304) dataset reports vastus-lateralis RNA-seq and Fried-related frailty measurements. Its submitted metadata contain 24 older adults labeled frail (including pre-frail), 30 fit older adults and 26 young adults. We compared only the first two groups; young samples were excluded. This is gene-expression context in bulk muscle, not a genetic score, a test of insomnia, a cell-type-specific result or proof of causal mediation. Cohort independence from every project source has not been individually established.

We selected FOXP2, SMAD5, RBM6 and MEIS1 because of earlier project candidate hypotheses; selection was not prespecified. We joined 80 unique sample prefixes between the released count matrix and GEO metadata and used only 54 older participants. Each gene's log₂ counts per million was regressed on group, centered age and sex using HC3 robust standard errors; BH adjustment was applied across these four postselected genes. This is a simple exploratory calculation, not a substitute for a fully specified RNA-seq differential-expression pipeline.

| Gene | Adjusted frail/pre-frail minus fit difference in log₂ CPM | HC3 SE | Nominal p | BH q within four genes |
|---|---:|---:|---:|---:|
| FOXP2 | −0.035 | 0.142 | 0.806 | 0.900 |
| SMAD5 | +0.013 | 0.066 | 0.841 | 0.900 |
| RBM6 | +0.032 | 0.117 | 0.786 | 0.900 |
| MEIS1 | +0.016 | 0.130 | 0.900 | 0.900 |

No association was detected in this small exploratory comparison; this does **not** show that these genes are unrelated to frailty or that effects are absent in other cells or tissues. Library-size scaling and the log transformation do not model count dispersion as DESeq2 or edgeR would; group labels combine pre-frailty and frailty, frailty measures overlap muscle function, batch and tissue cell composition are unadjusted, and the study is small. BH here corrects only the four genes in this post hoc screen and not the wider project search. Estimates do not test a variant-to-gene link or individual score performance.

The separate open [GSE117525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117525) series contains 259 expression samples but 187 distinct subject IDs, including paired training follow-ups. Its metadata report 61 (pre-)frail and 73 fit older baseline samples and 53 young baseline samples. Do not treat all 259 samples as independent or pool follow-ups with baseline. These data may support a separately designed replication after participant overlap and probe mapping are verified; no such replication is claimed here.

Reproduce the GSE144304 check with `work/geo_muscle/analyze_gse144304_candidates.py`; outputs are `outputs/stage6_gse144304_candidate_expression.json`. Compressed count matrix SHA-256 `a0867af116e23ea5d3aa5b84f20e6d0770d832f754c1fa6839753c823b2601ac`; series-metadata SHA-256 `4d79175238dfec042c5435f953689e7611a88d8ae280784d5b0b0e99ae1c0a0c`.

## Count-model follow-up for comparison with GSE167186

We reanalyzed the same 24 (pre-)frail and 30 fit older participants with edgeR 4.10.7, using `filterByExpr`, TMM normalization and robust quasi-likelihood modeling adjusted for age and sex. Of 58,283 count rows, 14,966 passed filtering; all four target Ensembl IDs were retained. The estimated frail-minus-fit log2 fold changes were FOXP2 −0.00565 (p=0.959), SMAD5 +0.00300 (p=0.946), RBM6 +0.0187 (p=0.833) and MEIS1 +0.0113 (p=0.915); each BH-adjusted p across these four genes was 0.959. [Results](stage6_gse144304_edger_candidate_results.csv); run `Rscript work/geo_muscle/gse144304_edger_comparison.R`. Original input hashes are above.

This provides **no observed FOXP2 replication** of the positive GSE167186 result, but a nonsignificant p-value does not prove no effect or a difference between studies. A direct heterogeneity test was not performed. Sarcopenia versus Fried-related pre-frailty/frailty, sex composition, covariates, populations and technical conditions differ, and participant/relative overlap is unverified. Both analyses are postselected, adjust only four genes, and have unresolved batch and cell-composition confounding. This is not a test of insomnia, genotype regulation or score prediction.

### Male-only sensitivity

Because the GSE167186 comparison was reported in men, we repeated the GSE144304 edgeR fit in its 13 men labeled frail/pre-frail and 15 men labeled fit, adjusting for age. This post hoc subset is small and was not used to select the original finding. FOXP2's frail-minus-fit log2 fold change was **−0.226** (nominal p=0.0772; BH p=0.309 across the four genes), and no other target passed that four-test threshold. [All four results](stage6_gse144304_edger_male_sensitivity.csv); run `Rscript work/geo_muscle/gse144304_edger_male_sensitivity.R`. The result differs in direction from the +0.402 estimate in GSE167186, but it is imprecise; we did not calculate a between-study heterogeneity test or a sex interaction. It cannot establish an opposite effect, a sex difference or a failed clinical-score prediction. Group definitions, age distributions, sequencing and unmeasured covariates remain different.

The published [GSE144304 cohort description](https://pmc.ncbi.nlm.nih.gov/articles/PMC10332038/) reports a Netherlands-based study in 54 older adults and gives the same 28-male count. GSE167186's source article describes a different skeletal-muscle study. This supports their being separately described studies, but we have not established that no participants or close relatives overlap, and phenotype definitions still differ.
