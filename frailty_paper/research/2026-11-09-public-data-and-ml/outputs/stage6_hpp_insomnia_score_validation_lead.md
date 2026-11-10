# Human Phenotype Project: potential insomnia-score validation cohort

**Status, 29 October 2026:** The Human Phenotype Project (HPP) has a documented individual-level genetic-score and sleep-questionnaire overlap. This is a **controlled-access lead, not verified free public data**. Its suitability for FI or five-component Fried frailty remains unestablished, and no project score has been evaluated.

## Evidence for a linked sleep analysis

* A published [obstructive-sleep-apnoea polygenic-score study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12171572/) used HPP data to test its own scores against several self-reported sleep measures, explicitly including a question on trouble falling asleep. It also studied sleep monitoring. This establishes that a genetic-score/sleep linkage was used in that particular study, but does not establish our exact score-marker coverage or the usable sample size for the falling-asleep question. A sleep symptom is not by itself a clinical insomnia diagnosis, and an OSA score is not either of our scores.
* The [EGA study record](https://ega-archive.org/studies/EGAS00001008040) describes an imputed autosomal genotype dataset, EGAD00010002714, for 8,958 participants, supplied as PLINK `.bed/.bim/.fam` files and imputed to a 1000 Genomes Phase 3 reference from low-pass sequencing. This is a dataset count, not a genotype-and-sleep complete-case count. The described PLINK 1 binary files are hard-call representation; ascertain whether dosage, variant INFO and exact multiallelic handling are available separately. The EGA entry is a controlled-access catalogue, not an open individual-level download.
* The OSA paper reports 5,950 HPP participants in its population description and gives a different subset for objectively measured apnoea; those counts should not be used as this project's eligible sample or as the count with all of our score markers. Distinguish objective sleep measurements, self-report and diagnosed insomnia.

## Access and important limits

[Pheno.AI's open innovation page](https://www.pheno.ai/open-innovation) describes competitive non-commercial academic access through a trusted research environment. Its timeline lists applications **1 September–31 December 2026**, review January–February 2027 and onboarding March–April 2027; elsewhere the same page refers to March–April 2026, so the timetable should be confirmed directly. The page does **not** state whether data access, computing, or other charges are zero. A separate [HPP access page](https://humanphenotypeproject.org/data-access) directs applicants to request information. Selection and institutional eligibility are not assured; no application or contact has been made.

The [sleep-study manuscript](https://weizmann.elsevierpure.com/ws/files/169042407/es_NatMed_PhenomeWideAssociationsofSleep_AM2024.pdf) says individual data cannot be made publicly available for privacy reasons. That statement and the TRE requirements must be respected; study-level research results and public documentation are distinct from permission to download or redistribute individual records.

## Practical checks before applying

1. Confirm zero total fees for this project and exact application timetable, eligible academic institution/supervisor and permitted analysis environment.
2. Request available questionnaire wording, coding and frequency distributions for trouble falling asleep and other insomnia symptoms, date/age and missingness; verify how many genotyped participants have each outcome and whether outcome collection is independent of score discovery and selection.
3. Check availability and quality of all 33 distinct project variants, exact GRCh37/38 coordinates and two-allele contrasts, sample-level calls/dosages, imputation accuracy, ancestry, relatedness, and overlap with discovery datasets. Do not infer variant coverage from the genome-wide dataset label.
4. Check whether a validated deficit FI or five-component FFS can be constructed, including temporal alignment and sleep-item overlap. Until verified, treat HPP primarily as an **insomnia-symptom** evaluation opportunity, not a frailty cohort. If frailty cannot be constructed, test the prespecified sleep outcome separately in a properly approved and independent subset.
5. Freeze the experimental score and comparator, select the primary outcome and model before accessing outcome associations, and report uncertainty, incremental performance and any unavailable markers. Assess ancestry transferability and sample selection; avoid clinical or causal claims.

The [CLS NCDS route](stage6_ncds_sleep_frailty_genetics_feasibility.md) has a published no-charge data-access statement and remains a stronger no-fee lead; HPP adds a distinct genotype/sleep opportunity whose cost remains unknown.
