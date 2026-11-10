# FI and Fried frailty: a same-cohort comparison

## Why this matters

The project's LDSC analysis estimates genetic correlation, not how similarly the two measures classify actual people. To assess the latter, I reconstructed a two-by-two table from a published study that measured both instruments in the same participants. This is external clinical evidence, not an analysis of the project's UK Biobank participants.

## Classification in the Toledo Study for Healthy Aging

The study included 1,563 community-dwelling Spanish adults aged 65 or older (mean 74.5). The authors defined frailty using a 40-item FI ≥0.25 and a Fried physical phenotype of at least 3 of 5 criteria. These specific thresholds and instrument versions matter; the result does not describe the correlation of the continuous scores, and cannot be transferred directly to the project's GWAS phenotypes. Source: [Phenotypes of frailty and their mortality and disability implications](https://pmc.ncbi.nlm.nih.gov/articles/PMC13264346/), methods and Table 1.

| | Fried frail | Fried not frail | Total |
|---|---:|---:|---:|
| FI frail | 103 | 125 | 228 |
| FI not frail | 18 | 1,317 | 1,335 |
| Total | 121 | 1,442 | 1,563 |

Reconstruction: the source also used a third measure, FTS-5. Both FI and Fried: 92 frail on all three + 11 frail on FI and Fried alone = 103. FI alone relative to Fried: 70 FI + FTS-5 + 55 FI alone = 125. Fried alone relative to FI: 16 Fried + FTS-5 + 2 Fried alone = 18. Neither FI nor Fried: 223 FTS-5 alone + 1,094 non-frail on all three = 1,317. 'Not frail' includes people classified as prefrail; it does not mean healthy.

Derived from those counts:

- 103/121 = **85.1%** of those classified frail by Fried were also classified frail by FI.
- 103/228 = **45.2%** of those classified frail by FI were also classified frail by Fried.
- 103/246 = **41.9%** were classified as frail by both among those classified as frail by either (Jaccard overlap).
- Overall binary agreement was 1,420/1,563 = **90.9%**, largely because 1,317 people were classified as not frail by both. The chance-corrected Cohen's kappa calculated from this table is **0.544**; kappa depends on prevalence and thresholds. Neither kappa nor Jaccard is the continuous phenotypic or genome-wide genetic correlation.

The asymmetry is meaningful: many people met the FI cutoff without meeting the physical phenotype cutoff. Neither measure is a gold standard here, so these numbers do not establish sensitivity or specificity or show that one instrument is correct for an individual.

## Subsequent adverse outcomes

The study's separate marginal models reported all-cause mortality hazard ratios for frail versus non-frail, adjusted for age, sex, Charlson comorbidity index, and polypharmacy: **2.80 (95% CI 1.93–4.06) for FI** and **1.73 (1.14–2.64) for Fried**. The groups being compared differ, the models are not necessarily mutually adjusted for both instruments, and no direct paired difference test is given here. Therefore these estimates do not prove FI has greater prognostic accuracy or that the hazard ratios differ statistically. Adding disability to the models changes estimates; disability can be part of the frailty pathway and conditioning on it changes interpretation. The published 5-year mortality model discrimination for the combined category model versus covariates was 0.7983 versus 0.77, with p=0.3 for that comparison; this is not a standalone FI-vs-Fried AUC comparison. Some combination groups were tiny, including two participants frail by Fried alone across all three instruments, limiting conclusions about rare categories.

For **worsening disability**, the same paper reports separately modeled adjusted odds ratios of **5.01 (95% CI 2.94–8.52) for FI frailty** and **3.05 (1.61–5.77) for Fried frailty**, relative to the respective non-frail groups. Its Table 3 states that the models include baseline disability, age, sex, Charlson index and polypharmacy. These are associations for different frailty classifications, not an adjusted head-to-head comparison or absolute probabilities; overlapping FI items and adjustment for comorbidity or disability also affect interpretation. The reported combined-category model AUC was 0.837 versus 0.7972 for the base model, with p=0.058, which does not establish a significant improvement at the conventional 0.05 threshold or prove that no improvement exists. These findings do not validate the genetic score or establish an insomnia association. [Published Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC13264346/).

A systematic review of eight studies making same-population comparisons found broadly comparable but modest mortality discrimination across FI and Fried-based phenotype measures; the continuous versions tended to outperform categorical versions, but studies differed in case mix, models, definitions and follow-up. This is mortality evidence, not a direct insomnia comparison. [Review](https://pmc.ncbi.nlm.nih.gov/articles/PMC9438224/).

## What this says for insomnia

This cohort does not establish which frailty measure is more strongly associated with insomnia. It does show that the clinical categories can identify different sets of people even when measured in the same population, which supports examining both rather than treating the measures as interchangeable. In the project's genetic results, the auxiliary FI–FFS genetic correlation of approximately 0.76 and the higher insomnia–FI than insomnia–FFS genetic correlations address a different question. In a shared participant sample, compare continuous FI and FFS and their binary cutoffs, then model insomnia against both with a paired uncertainty estimate. Repeat FI after removing sleep and then sleep plus fatigue, explicitly recording the resulting denominator and checking whether physical items also overlap. Observational and GWAS sensitivity analyses should be labeled separately.

## An observed FI–Fried classification correlation in a large cohort

[Pasquet et al.](https://link.springer.com/article/10.1007/s40520-024-02706-w) measured both approaches in the same Canadian Longitudinal Study on Aging participants. Among 28,527 included adults aged 45–85, they reported a **phi correlation of 0.38** between binary frail/not-frail classifications and **92.2% overall agreement**. Their FI had 93 deficits and classified FI >0.20 as frail; their Fried-based phenotype classified 3–5 of five components as frail. Multiple imputation was used for missing data. Their article specifies phi for categorical measures; this is not the correlation between the continuous scores. People classified not frail include the prefrail group.

This cohort and its definitions differ from the UK Biobank/TwinGene GWAS sources used in our project. Its 0.38 phi describes the association between two observed binary classifications; the project's approximately 0.76 r_g describes genome-wide shared genetic effects. Their numerical difference has no simple interpretation as a gap between environmental and genetic causes. Overall agreement can be high when most people are classified not frail, as the Toledo contingency table above illustrates. The Canadian paper reports prognostic models including age, sex and multimorbidity with c-statistics of about 0.67–0.69 for hospitalization and 0.79–0.80 for mortality across measures; these are performance of full models in that setting, not of FI or FFS alone, nor an insomnia risk prediction.

A second study of 505 older adults undergoing major noncardiac surgery reported a Spearman correlation of 0.57 between its FI and frailty phenotype and weighted kappa 0.45 (95% CI 0.38–0.51). The paper does not clearly establish in the quoted result whether the rank correlation used continuous measures or ordinal categories, so that detail should not be assumed. Its surgical FI excludes cognition and the clinical population differs markedly from both the Canadian sample and our discovery GWAS. [Deiner et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11090994/)

## Correlation of FI and Fried scores before binary cutoffs

A Spanish cross-sectional study of older adults measured a 27-deficit FI fraction and the count of positive Fried criteria. Its [Figure 2](https://link.springer.com/article/10.1007/s11357-025-01931-z/figures/2) displays their **Spearman rank correlation as 0.542** (figure marker corresponds to p<0.01). This is a correlation of an FI score with the discrete 0–5 phenotype count, rather than a phi correlation of two frail/not-frail classifications. I read the value directly from the figure; the text confirms the variables and the Spearman method.

The abstract, methods and results report 150 participants, whereas one introductory sentence says 155; the exact number used for this particular correlation is not reported in the figure and should not be assumed. Participants were aged at least 65, recruited in Galicia and subject to multiple exclusions, including history of depression, cancer or autoimmune disease. The study is small and cross-sectional, its measures differ from the project's GWAS definitions, and the correlation has no direct implications for direction of effect or prediction. This measured rank correlation cannot be compared numerically with the project's genome-wide genetic correlation as if both measured the same thing. [Full article](https://link.springer.com/article/10.1007/s11357-025-01931-z).
