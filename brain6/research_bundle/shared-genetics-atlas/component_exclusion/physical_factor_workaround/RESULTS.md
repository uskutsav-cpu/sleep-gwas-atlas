# A sleep-item-free physical genetic factor: exploratory workaround

**Further sensitivity:** The [residual-covariance analysis](RESIDUAL_RESULTS.md) shows the latent-factor association is not identified if all indicator-specific links are free. The [explicit equal-weight genetic index](INDEX_RESULTS.md) supplies a separately defined physical-function sensitivity without that zero-residual restriction; its post-hoc and source limitations remain.

**Updated 6 November 2026:** A subsequent [joint fit](JOINT_RESULTS.md) re-estimates the loadings and both insomnia correlations together. It gives −0.087 (SE 0.045) for registry insomnia, with an approximate interval including zero, and −0.401 (0.032) for questionnaire insomnia. Prefer that sensitivity over the fixed-weight two-stage projection below; neither is a full validated genomic SEM or matched item-deletion GWAS.

## What this does and does not do

We constructed an **analytic latent genetic factor** from publicly downloadable GWAS of self-reported usual walking pace, left hand grip strength and self-reported days of moderate activity. None of these three input questions or measurements explicitly includes sleep or exhaustion. This is a way to ask whether insomnia shares genetic signals with a *common physical-function dimension*. It does not recreate a Fried frailty score, remove an item from the same people, measure individual frailty, or show that the inputs are biologically independent of sleep. Further, this is a calculation from genetic-correlation summaries, not a full multivariate factor GWAS.

Using [Neale Lab LDSC-formatted inputs](https://nealelab.github.io/UKBB_ldsc/downloads.html), we aligned the three indicators with the existing FinnGen registry-insomnia input on **1,053,614 nonpalindromic SNPs**. The three source file hashes matched the earlier audits. Every fit retained the same marker count and the 200 heritability delete vectors were identical across fits for the same trait. Direction assumes the published Neale processing conventions: higher values represent faster walking, greater left grip, and more moderate-activity days. [UKB's walking categories](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100313) run from slow to brisk; walking pace was not collected for people unable to walk. Activity's source-code direction audit is [here](../moderate_activity_ldsc/DIRECTION_AUDIT.md). Raw participant-level phenotype processing has not been independently reconstructed.

## Pairwise correlations and physical factor

| Comparison | Genetic correlation (SE) |
|---|---:|
| Walking pace–grip | 0.2931 (0.0222) |
| Walking pace–activity | 0.1306 (0.0298) |
| Grip–activity | 0.0848 (0.0265) |
| Registry insomnia–walking | −0.0631 (0.0319) |
| Registry insomnia–grip | −0.0896 (0.0352) |
| Registry insomnia–activity | 0.0051 (0.0421) |
| Questionnaire insomnia–walking | −0.2631 (0.0261) |
| Questionnaire insomnia–grip | −0.1959 (0.0249) |
| Questionnaire insomnia–activity | −0.0037 (0.0324) |

The three physical correlations imply standardized one-factor loadings of **0.672 walking, 0.436 grip, and 0.194 activity**, oriented toward greater function/activity. All 200 leave-block calculations were admissible under the chosen loading checks. Three indicators alone exactly determine a one-factor correlation model: they do **not** test its fit. The modest activity loading makes the construct largely driven by walking and grip; the estimated squared loading for activity is about 0.038, not an individual-level variance explained.

We divided each insomnia–indicator correlation by the corresponding loading and combined the three implied factor correlations using a generalized least-squares estimate based on paired 200-block jackknife covariance. Under the **strong assumption that each indicator's residual genetic component is uncorrelated with insomnia**, the exploratory registry-insomnia projection was **−0.112 (conditional SE 0.045)**. The three implied values were −0.094, −0.205 and +0.026; their exploratory two-degree-of-freedom heterogeneity statistic was 1.95 (p=0.377). Failure to reject heterogeneity does not establish the assumption or adequate power.

## Post-hoc alternative insomnia source

As a source sensitivity added after inspecting registry results, we downloaded [Neale Lab's ordinal UKB sleeplessness/insomnia GWAS](https://nealelab.github.io/UKBB_ldsc/h2_summary_1200.html). Its [question and response coding](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1200) describe increasingly frequent insomnia; exact processing was not rerun from individual-level records. It retained the same common markers. The same physical factor gave a projection of **−0.411 (conditional SE 0.033)**, with indicator-specific projections −0.392, −0.449 and −0.019. The exploratory heterogeneity statistic was 5.91 on two degrees of freedom (p=0.052); this is not a validated model-fit test or evidence that the common-factor restriction is true. The activity indicator received about 1.4% of the estimated weight, with walking and grip accounting for about 98.6%.

The post-hoc questionnaire-minus-registry contrast was **−0.299 (paired block SE 0.050)** when each source's full-data GLS weights were held fixed (unadjusted normal approximation p≈2.9×10⁻⁹). That uncertainty includes the shared physical loading estimates through the block calculations but does not include uncertainty from estimating the weights, model selection, GWAS source assumptions or the full research search. Questionnaire and physical GWAS share UK Biobank, while registry insomnia is Finnish; phenotype, ancestry, ascertainment, sample overlap and reference fit change together. The contrast cannot identify why the estimates differ and is not an independent replication. It is also not a causal effect or a test of item removal.

## Quality limits and reproducibility

The physical–physical cross-trait LDSC intercepts were 0.1635 for walking–grip, 0.1074 for walking–activity and 0.0604 for grip–activity; they may reflect shared samples or other model effects, not measured overlap counts. Registry insomnia's univariate intercept was 1.0778 (SE 0.0096). The UKB indicator GWAS are related and have uncertain source-level overlap. General-European LD may fit Finnish samples imperfectly. The physical indicators and questionnaire insomnia share UKB; no independent cohort validation was performed. These estimates use six and then three separately fitted LDSC pairs rather than a complete multivariate LDSC covariance model and Genomic SEM, and summary point correlations are rounded to four decimals before factor calculation. The conditional SEs hold estimated GLS weights fixed. The factor cannot be safely interpreted as a clean causal sleep-independent phenotype, and no local association or individual risk prediction follows.

The [protocol](PROTOCOL.md), [verification](VERIFICATION.md), [preparation](prepare.py), [questionnaire preparation](prepare_questionnaire.py), source and output hash receipts, [LDSC runs](run.sh), [questionnaire runs](run_questionnaire.sh), [factor evaluation](evaluate.py), [paired source contrast](compare_sources.py), logs and delete values document the analyses. Full source and prepared files are retained under `/Volumes/Extreme SSD/brain6-physical-factor-20261101/`. This comparison used no AlphaGenome output and did not train a model. A full documented same-person sleep/exhaustion-excluded frailty GWAS, source-valid multivariate fit and independent validation are still needed for the original question.
