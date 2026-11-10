# Public HUNT and Partners follow-up of insomnia–frailty candidate variants

29 October 2026. This is an exploratory **variant-level insomnia follow-up**, not validation of the project's genetic scores or of a frailty outcome.

## Source and coverage

[Lane et al., Nature Genetics (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6415688/) published selected insomnia-locus results for the HUNT study (14,923 cases, 47,610 controls; self-reported insomnia symptoms) and Partners Biobank (2,217 cases, 14,240 controls; physician-diagnosed insomnia). Their publicly downloadable [supplementary workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC6415688.1/NIHMS1519714-supplement-2.xlsx), Tables S7 and S8, provides individual variant results, not participant genotypes or full-genome HUNT and Partners summary statistics. The published cohorts' relationships to all discovery and selection datasets used for these project scores have not been exhaustively established; treat independent replication as unconfirmed.

An exact-rsID match found six of the 33 unique candidates in both tables: five of the 15-marker prototype and two of the 22-marker equal-weight count; rs4886140 is in both. Absence from these selected-locus tables does not mean an association is null. The [complete 33-variant lookup](stage6_lane2019_hunt_partners_candidate_lookup.csv) includes nonmatches and the source and project allele pairs. Original nearest-gene labels describe proximity in the source and do not establish causal genes.

| Variant | Source effect / other allele | HUNT OR (95% CI), p | Partners OR (95% CI), p | Prototype membership |
|---|---|---|---|---|
| rs113851554 (near MEIS1) | T/G | 1.13 (1.06–1.20), 0.000182 | 1.28 (1.12–1.45), 0.0002891 | 15 |
| rs1942262 | A/G | 1.05 (1.01–1.08), 0.00414 | 1.06 (0.99–1.14), 0.081 | 15 |
| rs4886140 | G/A | 1.04 (1.01–1.07), 0.0119 | 1.02 (0.96–1.10), 0.508 | 15 and 22 |
| rs9845387 | C/A | 1.09 (1.02–1.16), 0.0108 | 1.07 (0.89–1.29), 0.468 | 22 |
| rs35881094 | G/T | 1.02 (0.99–1.05), 0.177 | 1.05 (0.99–1.12), 0.108 | 15 |
| rs10865954 | T/C | 1.01 (0.98–1.04), 0.462 | 0.98 (0.91–1.05), 0.532 | 15 |

Four HUNT associations and one Partners association have unadjusted p < 0.05. The published insomnia loci and these project candidates were selected before this comparison, so these numbers are descriptive, not a prespecified multiple-testing-corrected replication claim. Most estimates are imprecise; a nonsignificant result does not prove no effect. The two outcomes use different ascertainment. For rs4886140, the 15-marker prototype counts A and has a negative insomnia weight, while the 22-marker count uses G; the table's OR is for G. Some candidate rsIDs, including rs113851554 and rs4886140, have multiallelic reference representations: an rsID and matching letters alone do not establish a safe dosage conversion. Verify coordinates, build, strand, exact allelic contrast and participant overlap before any scoring or meta-analysis.

The MEIS1 result provides external published support for an insomnia association at this candidate, but does not independently establish a frailty association, a shared causal variant, or predictive performance of any score. Because only six selected variants are available, these tables cannot supply a valid evaluation of a 15- or 22-variant score, its calibration, discrimination, or the comparison between FI- and FFS-based weights. A genuinely held-out cohort with linked genotype and sleep and frailty phenotypes remains necessary; free-access registration or controlled access may be required.

## Reproduction

Downloaded workbook: `work/lane2019/supplement.xlsx`, SHA-256 `957ae2c2a3d15bcf9a29e313d9283a6c03d73439159f96116e27fa726d8e6919`. Run `python3 work/lane2019/extract_candidate_replication.py` with `openpyxl` installed. The script reads both prototype lists, checks identical source allele coding across tables, and exports every candidate including nonmatches. ORs and p-values are copied from the published supplementary tables, not re-estimated. Results are observational associations, not causal conclusions.
