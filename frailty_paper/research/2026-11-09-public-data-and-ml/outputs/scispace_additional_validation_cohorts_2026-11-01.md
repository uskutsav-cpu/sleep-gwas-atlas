# Additional cohorts for sleep–frailty genetic-score validation

Reviewed 1 November 2026 using SciSpace searches, original cohort publications and data-custodian documentation. These are **leads**, not approved access, verified marker coverage, or completed score validations. The project's scores were developed largely in European-ancestry sources, so a new target requires ancestry-aware evaluation; scores cannot be assumed transportable.

A [cohort-specific access and analysis plan](nicola_halst_access_and_analysis_plan.md) specifies the proposed data, feasibility checks, and what must be frozen before any target association analysis.

## 1. Northern Ireland Cohort for the Longitudinal Study of Ageing (NICOLA)

**Why it matters:** NICOLA collected questionnaires and clinical measurements from community-dwelling adults predominantly aged 50+. A [genotype data paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13174777/) describes genetic data after initial QC for 2,978 participants, with 551,839 directly typed and 18,148,478 imputed variants. These are cohort-level data, not the number with complete sleep, frailty and genotypes. The paper describes two genotyping batches and relatedness filtering; verify the actual release, build, dosage/format and quality for every candidate. [EGA record EGAD00010002762](https://ega-archive.org/studies/EGAS00001007915) lists the genetic release.

A [published 30-deficit FI comparison](https://academic.oup.com/eurpub/article/34/4/710/7632123) uses NICOLA data and explicitly includes daytime sleepiness. Thus, the full FI is not free of sleep overlap; inspect the original item list and establish whether an appropriate no-sleep/no-fatigue FI can be rebuilt in the same people. The [NICOLA cohort and health-assessment profile](https://pmc.ncbi.nlm.nih.gov/articles/PMC9999338/) documents physical and health measurements. A [mental-health measure catalogue](https://www.cataloguementalhealth.ac.uk/NICOLA) lists sleep problems; this is not proof of a specific insomnia definition in every genotyped participant. Confirm question names, visit, case counts and data-quality rules.

**Access:** NICOLA's [data-access policy](https://www.qub.ac.uk/sites/NICOLA/FileStore/DataAccessPolicy2025_23092025%20V2.pdf) states that EGA molecular data are not linked to phenotype information, and describes use through an approved safe setting; its [researcher page](https://nicola.qub.ac.uk/FileStore/) says remote secure access can be facilitated. An EGA download alone is not a linked score-validation dataset. Ask the custodian whether approved genotype–phenotype linkage can occur in their safe setting, actual joint counts, project eligibility, and total cost; its policy allows cost-recovery charges, so no-fee use is not established. No application or contact was made. Independence and close-relative overlap with all discovery/follow-up sources are unverified.

**Priority:** High scientific relevance, medium feasibility until linkage and fees are confirmed. 2,978 is a ceiling for this genotyped release, not a power calculation or a joint analysis N.

## 2. Healthy Aging Longitudinal Study in Taiwan (HALST)

**Why it matters:** A [published HALST study](https://ir.nhri.edu.tw/handle/3990099045/14916?locale=en-US) examined self-reported insomnia and Fried-criteria frailty longitudinally in 4,744 participants, mean age about 69.4. This demonstrates those phenotypes were studied together; it does **not** show how many of those participants have genotypes. The [HALST cohort profile](https://pmc.ncbi.nlm.nih.gov/articles/PMC5837206/) describes older adults and clinical assessments. A separate [HALST genetics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12788161/) documents genome-wide TWBv2 array data, while a [different study](https://diabetesjournals.org/diabetes/article/75/3/582/164189/Integrating-Polygenic-Scores-and-Family-History) reports 1,968 HALST samples for diabetes score replication. None of these counts establishes the insomnia–frailty–genotype intersection or coverage of the present scores.

**Access:** The Taiwan NHRI's [HALST resource description](https://ghd.nhri.org.tw/Home/Datasets_card?result=Key_161) describes a data application, research plan and IRB documentation; its [cohort profile](https://pmc.ncbi.nlm.nih.gov/articles/PMC5837206/) encourages collaboration. Availability to this team, linkage, cost, genotype release and use conditions are not verified. No request has been submitted. The physical phenotype includes exhaustion; a no-exhaustion sensitivity is a different construct and requires prespecification. As an East Asian target, this could assess score portability but cannot be assumed to reproduce European-ancestry effects.

**Priority:** High phenotype relevance, potentially useful non-European follow-up; cost and practical access are unknown. Compare outcome and discovery ascertainment, relatedness and independence before any replication claim.

## 3. Mexican Health and Aging Study (MHAS)

SciSpace surfaced MHAS's [genetic-data profile](https://pmc.ncbi.nlm.nih.gov/articles/PMC11653003/), reporting 7,204 people with genome-wide data across study collections, a number that is not an eligible joint analysis count. MHAS documents sleep and published frailty measures. Genomic access differs for blood and saliva and requires approved routes; exact sample linkage, costs, variant coverage and study independence remain unresolved. **The [MHAS data-use notice](https://www.mhasweb.org/DataProducts/Linkages.aspx) prohibits LLMs and other AI tools from managing, processing or analyzing MHAS-distributed microdata.** A human-led analysis may be possible under the relevant agreements; do not put MHAS microdata into this AI environment. [Prior detailed review](mhas_sleep_frailty_genetics_feasibility.md).

## 4. oriGen, Mexico: a non-free future/collaboration option

The [oriGen project](https://tec.mx/es/investigacion/proyecto-origen/investigadores) describes 100,000 adult participants and emerging genomic resources. A [2026 WGS data paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13642253/) says the 1,427-genome dataset is restricted, requires collaboration and fees during its embargo through August 2028, and has possible additional computing costs. These exact terms pertain to that WGS dataset; other oriGen products can have different policies. A published sleep–frailty–genotype intersection, candidate-marker coverage and a free route for this team are not established. Exome sample counts cannot be assumed to contain noncoding score variants. **Not recommended as a current free-data route.**

## Practical ranking and next action

1. Check NICOLA and HALST with a qualified researcher/custodian for genotype–sleep–FI/physical-frailty overlap, individual score-variant coverage, cohort/relative overlap, access terms and *total* cost. NICOLA's EGA data do not by themselves include the linked survey phenotypes.
2. Compare with the existing no-data-access-fee [Framingham](framingham_mentor_access_proposal.md), [NHATS](stage6_nhats_genetic_validation_feasibility.md) and [NCDS](stage6_ncds_sleep_frailty_genetics_feasibility.md) requests, all conditional on approval and linkage; no study has yet been shown usable for our exact analysis.
3. Treat MHAS as human-led under its actual use terms. Do not use oriGen to satisfy a no-fee requirement without a verified exemption.

No data access request, outreach, participant-data download or score analysis was performed. Do not present sample sizes from separate publications as the number of people with all required observations. An independent association, discrimination or calibration assessment remains outstanding.
