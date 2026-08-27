# Selective source-provenance harvest — 2026-08-26

The remote `claude/neuro-immune-cardio` branch expands the project to 149
traits, so it cannot be merged into the locked atlas-v1.0 scope. Commit
`8023d143f2bde431ce0f54fa827302cd9d37adcd` nevertheless contains useful
file-level provenance from complete downloads that were later evicted to save
disk. This review considered only the 28 source-pending traits in the locked
45-trait manifest.

## Promoted exact matches

Eleven selected releases matched the locked publication, phenotype stratum,
sample-size intent, European ancestry, and GRCh37 requirement. Their official
GWAS Catalog landing/download URLs, exact byte counts, SHA-256 checksums, and
branch-recorded build-anchor validation are now in
`config/public_gwas_sources.tsv`:

- de Lange 2017 combined IBD, Crohn's disease, and ulcerative colitis;
- Ishigaki 2022 European-only rheumatoid arthritis;
- Han 2020 European-only asthma;
- Yengo 2018 GIANT+UKB BMI;
- Graham 2021 European HDL, LDL, and log-triglycerides strata;
- Mishra 2022 GIGASTROKE European any stroke; and
- Nielsen 2018 atrial fibrillation.

This branch harvest advances source verification from 17/45 to 28/45. A
separate direct verification of Campos 2020 snoring brings the current total to
29/45. Subsequent direct verification of the selected colorectal- and
lung-cancer releases brings the current total to 31/45. These counts do not
claim that the archives are currently materialized. Direct verification of the
FinnGen R9 sleep-apnoea endpoint subsequently brings the total to 32/45; its
GRCh38 build remains a separate harmonization blocker. These counts also do not
claim that every header is supported by the harmonizer, binary prevalence is
resolved, or LDSC is ready. Shared samples are
retained as explicit analysis concerns: IBD/CD/UC overlap, the three lipid
traits share cohorts, and BMI overlaps UKB-derived sleep phenotypes.

## Deferred branch candidates

These records were not promoted:

| Trait | Reason |
|---|---|
| insomnia | Branch candidate is Watanabe 2022 rather than the locked Jansen 2019 UKB-only release; the selected source/phenotype decision must be reviewed explicitly. |
| chronotype | Branch candidate is a continuous BOLT chronotype measure, while the locked trait is binary morningness with different sample counts. |
| Alzheimer disease | Branch substituted smaller Kunkle 2019 GRCh37 data for Bellenguez 2022 because the latter is GRCh38. Build standardization alone is not grounds to discard the selected larger GWAS. |
| Parkinson disease | The branch could not independently establish whether the accession is the appropriate non-23andMe/proxy-case release. |
| major depression | Branch substituted Wray 2018 for the selected, larger Howard 2019 release because the latter is controlled. Ease of access is not the selection criterion. |
| schizophrenia | Accession/publication/phenotype notes conflict and require primary-source reconciliation. |
| bipolar disorder | Branch substituted a roughly five-fold smaller 2016 GWAS for Mullins 2021. |
| ADHD | Branch substituted a later accession whose phenotype and relationship to the locked Demontis 2023 source require review. |
| multiple sclerosis | Branch substituted a roughly ten-fold smaller 2016 German GWAS for IMSGC 2019. |
| type 2 diabetes | Branch substituted Scott 2017; although genome-wide, its bare coordinate marker lacks the rsID mapping needed by the current LDSC path. |
| coronary artery disease | The branch describes the selected file as EUR plus ancestry-not-reported, so it does not yet satisfy the EUR-only contract. |
| systolic blood pressure | Branch substituted Keaton 2024 for Evangelou 2018; overlap, phenotype transformation, and selection implications require review first. |

No branch candidate record existed for snoring, sleep apnea, colorectal cancer,
lung cancer, or melanoma. Snoring, colorectal cancer, and lung cancer were
resolved independently from official GWAS Catalog sources, and sleep apnea was
resolved from FinnGen's R9 manifest, endpoint definition, and data dictionary.
The selected Landi melanoma accession `GCST010304` has `fullPvalueSet=false` in the official API
and no directory in the Catalog summary-statistics FTP tree, so it remains
pending rather than being represented by an invented download. Sleep apnea
requires a validated build-conversion path; melanoma still requires new
file-level primary-source work.

The rule for the next pass is unchanged: prefer the largest appropriate and
scientifically compatible release, even when it requires controlled access,
documented liftover, or explicit variant-ID mapping. A smaller public file is a
candidate or sensitivity dataset, not an automatic replacement.
