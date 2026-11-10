# Public insomnia-score PheWAS: physical-function context

Checked 29 October 2026. The [ExPRSweb Michigan Genomics Initiative (MGI) PheWAS](https://csg.sph.umich.edu/larsf/exprsweb/scores/exprsweb-insomnia-1200-plink-mgi-freeze3-20211120/phewas/) publishes 1,685 diagnostic association rows for a previously developed insomnia polygenic score. This is [PGS003326](https://www.pgscatalog.org/score/PGS003326/), a 27,462-variant clumping-and-thresholding score using a UK Biobank insomnia GWAS, **not either project score**. The portal reports score tuning/evaluation in MGI; the PheWAS is not an independent held-out score-development cohort, and UK Biobank discovery overlaps sources used in our work. It does not provide individual genotypes or evaluation of our custom score.

Four diagnostic labels were selected **after inspecting** the table as frailty-component context. These are billing-code-derived disease/symptom phenotypes, not a measured Fried frailty score or deficit FI. In the full-cohort PheWAS, the reported score OR and p-value are:

| Phecode | Label | Cases / matched controls | OR (95% CI) | p |
|---|---|---:|---:|---:|
| 350.2 | Abnormality of gait | 3,920 / 30,520 | 1.01 (0.974–1.04) | 0.664 |
| 772.3 | Muscle weakness | 2,949 / 28,064 | 1.02 (0.985–1.06) | 0.225 |
| 798 | Malaise and fatigue | 11,539 / 35,196 | 1.02 (0.998–1.04) | 0.0773 |
| 798.1 | Chronic fatigue syndrome | 1,823 / 17,795 | 1.05 (1.00–1.11) | 0.0341 |

The reported confidence limits are rounded. The last result is nominal only and does not survive a conservative Bonferroni threshold of 0.05/1,685 ≈ 2.97×10⁻⁵; selection of these labels after viewing the data adds further interpretive caution. Absence of a statistically clear signal in the other rows is not proof of no effect. These diagnoses capture healthcare coding and ascertainment, not standardized gait speed, grip, exhaustion, or prospective frailty, and cannot establish causal effects. The precise score scaling and adjustment scheme should be confirmed in the source methods before comparing effect magnitudes with another study.

The source all-subject file is [here](https://csg.sph.umich.edu/larsf/prsweb-data/exprsweb/releases/exprsweb_20211120/phewas/ExPRSweb_insomnia_1200_PLINK_MGI-Freeze3_20211120/all.tsv), SHA-256 `fa887d69286e8887e592604ec874102655ee880736b1147543688c18266baf6e`; [extracted rows](stage6_exprsweb_insomnia_function_candidates.csv). The separate exclusion analysis is reviewed below. The resource is complementary published-score context, not a substitute for scoring our variants in an independently ascertained sample.


## Published exclusion sensitivity

The authors also provide an [exclusion PheWAS](https://csg.sph.umich.edu/larsf/prsweb-data/exprsweb/releases/exprsweb_20211120/phewas/ExPRSweb_insomnia_1200_PLINK_MGI-Freeze3_20211120/exclusion.tsv) that excludes participants with the score-defining insomnia diagnosis; the [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9606385/) explains the general exclusion design. Among the same post hoc selected diagnostic labels, the resulting odds ratios were 1.02 for gait abnormality (p=0.527; 1,641 cases), 1.01 for weakness (p=0.797; 1,306 cases), 1.01 for malaise/fatigue (p=0.589; 5,028 cases) and 1.01 for chronic fatigue syndrome (p=0.832; 607 cases). [Extracted rows](stage6_exprsweb_insomnia_function_exclusion.csv). The chronic-fatigue estimate is less precise after exclusion: 95% CI 0.929–1.10 versus 1.00–1.11 in the full table. The change is descriptive: these are nested, differently matched analyses with different case and control counts, not an independent replication or a formal test of the difference. Conditioning on a diagnosis may induce selection, and the result does not prove that insomnia mediated the association or that no genetic relationship exists.

The exclusion file has 1,685 rows and SHA-256 `0b771d3b0644eee828a2c4afb6b1d5b60c2f1333fe81e830e3bbad289549ac47`. Neither analysis estimates performance of the custom project scores.
