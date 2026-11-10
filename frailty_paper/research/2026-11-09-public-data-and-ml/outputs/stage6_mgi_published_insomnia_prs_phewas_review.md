# Published insomnia score in Michigan: secondary phenotype context

[ExPRSweb](https://csg.sph.umich.edu/larsf/exprsweb/scores/exprsweb-insomnia-1200-plink-mgi-freeze3-20211120/) publishes an existing insomnia polygenic score and a [phenome-wide association table](https://csg.sph.umich.edu/larsf/exprsweb/scores/exprsweb-insomnia-1200-plink-mgi-freeze3-20211120/phewas/) from the Michigan Genomics Initiative (MGI). This is a different score from both project prototypes: it contains 27,462 SNPs, uses UK Biobank insomnia GWAS weights and was tuned in MGI. The target cohort is a separate named US healthcare biobank; the project has not established absence of participant or relative overlap with all discovery and selection sources. This is a published secondary analysis, not a new test of either project's score.

The 1,685-row full PheWAS table contains no phenotype explicitly labelled frailty. Selected related results, with the portal's rounded odds ratios and intervals, are:

| EHR-derived phecode | Cases / matched controls | Score OR | 95% CI | p |
|---|---:|---:|---:|---:|
| 327.4 Insomnia | 5,524 / 30,960 | 1.08 | 1.05–1.11 | 1.42 × 10⁻⁷ |
| 772.3 Muscle weakness | 2,949 / 28,064 | 1.02 | 0.985–1.06 | 0.225 |
| 798 Malaise and fatigue | 11,539 / 35,196 | 1.02 | 0.998–1.04 | 0.0773 |

The [method documentation](https://csg.sph.umich.edu/larsf/exprsweb/method/) describes centered and scaled scores, phenotype-specific matching and adjusted models. The EHR labels are not measured grip, FI or FFS; uncertain weakness and fatigue estimates do not prove absence of an effect. The outcomes and score were inspected after seeing the results, and multiple correlated phenotypes and many score models exist. This is not a prespecified test of a shared genetic mechanism or proof of causality. The PheWAS sample and analyses should not be conflated with the score-detail page's separate 50/50 primary-trait evaluation; MGI was also used for tuning.

The [original MGI resource paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9932985/) states that the individual-level data are not public and that downloadable GWAS summaries require a data use agreement. Its interactive PheWeb has registration and anti-scraping conditions. The table above comes from ExPRSweb's explicitly downloadable published PheWAS output, not from scraping MGI's restricted browser. It does not supply arbitrary variant association results or FI/FFS outcomes for project-score validation.

For provenance, the [public all-subject TSV](https://csg.sph.umich.edu/larsf/prsweb-data/exprsweb/releases/exprsweb_20211120/phewas/ExPRSweb_insomnia_1200_PLINK_MGI-Freeze3_20211120/all.tsv) was downloaded to `work/stage6_mgi_insomnia_score_phewas.tsv`; its local SHA-256 is `fa887d69286e8887e592604ec874102655ee880736b1147543688c18266baf6e`. This records the retrieved copy, not an independent publisher checksum. No individual data were accessed.
