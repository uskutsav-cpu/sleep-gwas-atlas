# Independent outcome validation: cohort feasibility

## GEMINI public summary statistics

The [GEMINI release check](stage6_gemini_insomnia_public_release_review.md) identifies a free insomnia GWAS, but its metadata describe it as UK Biobank-only. It is therefore not independent validation of the project's UK Biobank-derived insomnia discovery score; no frailty-labeled file was identified in its release inventory.

## Canadian Open Parkinson Network lead

[C-OPN](stage6_copn_parkinson_validation_feasibility.md) reports sleep history, genotyping and a project to develop a frailty index. A 2025 announcement described free data access for academic members, subject to network conditions; the exact accessible genotypes, FI, linked sample, permitted use and current total cost are unverified. Its disease-enriched sample is a specialized follow-up rather than a general-population validation cohort.

## GEHAS circadian cohort lead

[GEHAS](stage6_gehas_sleep_frailty_genetic_feasibility.md) describes a 596-person Chinese circadian subcohort with sleep monitoring and subcohort genotyping, alongside separately reported frailty-related research. It offers a potentially relevant non-European research lead, but exact linkage, outcome and candidate coverage, access and costs remain unverified; its data are available on request, not established as an immediately open validation dataset.

## Czech HAPIEE lead

The 2026 [HAPIEE-CZ assessment](stage6_hapiee_cz_validation_feasibility.md) identifies documented sleep measures, a published frailty index, genotyping and a cohort statement that it does not charge for data preparation and release. Access is by reviewed request and agreement, and the genetic release, linkage, exact project costs, candidate coverage and usable joint sample remain unverified. It is not an immediately open dataset or a completed validation.

## HCP-Aging / AABC lead

The 2026 [AABC genetic release](stage6_hcp_aging_aabc_validation_feasibility.md) reports TOPMed-imputed genotype data for 1,024 participants, while sleep, gait and grip assessments are described separately. This creates a potentially useful sleep/physical-function follow-up, but individual genotypes require approval, linkage and exact score coverage are unverified, and no full FI or Fried phenotype has been established. Download costs are covered for many NIAGADS files, but project-specific total costs have not been established. It is not an immediately available open score-validation dataset.

The separate [enhanced NKI-Rockland lead](stage6_enki_sleep_motor_genetics_feasibility.md) has published combined sleep-quality and motor data, but a currently available linked genotype release and full frailty outcome have not been established; phenotypic access requires a data-use agreement.

## Current project constraint: public, free data only

The user has confirmed that the team does not have access to an independent participant-level cohort and wants publicly available free data. The controlled-access cohort discussion below is background from an earlier feasibility assessment, **not the current recommended path**. Continue with public GWAS summary statistics and openly available reference panels; these support association and robustness analyses, but they do not establish individual-level prediction, calibration, or clinical validity. A suitable fully public, individual-level dataset with genotypes, the relevant outcomes, and adequate independent sample size has not been established.

An [HRS precomputed-score check](stage6_hrs_precomputed_scores_feasibility.md) found no sleep, insomnia or frailty label in the current public included-score list. The product page contains both a public-release description and access-grant/order-form instructions, so current conditions need verification. Existing unrelated scores cannot reconstruct either project prototype.

The [Wisconsin Longitudinal Study](stage6_wls_validation_feasibility.md) has published sleep and physical-frailty composite analyses, a genome-wide genotyped cohort and a version-2 insomnia-score lead. The published four-component composite is not a deficit-accumulation FI or the full Fried phenotype. Linked genetic data require approval and an institutional agreement; costs, exact shared participants and candidate coverage are unverified. WLS prohibits uploading its data to public LLMs or AI platforms. It is not an unrestricted open validation dataset.

The newly released Sleep Health Data Commons and Human Sleep Project are credentialed sleep and health data leads, but their public descriptions do not establish linked genotypes or FI/FFS for score validation. The Human Sleep Project also describes EEG embeddings that might support a separate future ML question. See the [access and feasibility review](stage6_bdsp_sleep_commons_feasibility.md).

## MrOS sleep and genetics lead

The [MrOS Sleep Study](https://sleepdata2.org/datasets/mros/pages) reports 3,135 older male participants with polysomnography and actigraphy, and the [parent MrOS study](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000373.v1.p1) reports baseline DNA and QC-cleaned PLINK files. This makes the cohort scientifically relevant for investigating sleep and frailty, but the dbGaP genetic files are under authorized access; its public FTP provides manifests and release notes, not open linked individual-level genotypes. [MrOS Online](https://mrosonline.ucsf.edu/) requires an account and data-use agreement for clinical downloads. The [dbGaP study-specific data-use addendum](https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?stacc=phs000373.v1.p1&view_pdf=) lists sleep disorders and decline in physical function including frailty among permitted research domains; this is a scope indication, not access approval. A no-data-fee dbGaP route is documented below, but participant-level access, linked phenotype coverage, candidate-variant coverage, discovery-sample independence and total project costs remain unverified. MrOS is not an immediately usable open validation dataset.

**Consent-group clarification, 29 October 2026:** The linked addendum has different restrictions by consent group. Groups 1 and 3 list sleep and physical-function or frailty-related research, while group 2 is limited to osteoporosis/fracture, osteoarthritis and obesity. The project's eligibility for each group and any intersection with sleep/frailty data must be decided through the data-access process, not inferred for all MrOS participants from the study-wide description. A [published analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC12089207/) reports 5,143 imputed-genotype records and 5,130 phenotype records, but these counts do not establish the number eligible, overlapping the sleep visit, or usable for this project.

**Phenotype clarification:** A [published MrOS Sleep analysis](https://pubmed.ncbi.nlm.nih.gov/19793160/) assessed both sleep and a modified Fried physical frailty phenotype in the *same* 3,133 men aged at least 67; 437 were classified as frail. It used self-reported and objectively measured sleep. This corrects the earlier uncertainty about whether sleep and a physical-frailty phenotype were ever jointly measured; genotype linkage and consent-specific coverage for this sample remain unverified. The definition includes poor energy and differs from our discovery FFS, and the male-only study does not establish generalizability to women. A [separate data-availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC9881675/) describes NSRR MrOS polysomnography and clinical data as free subject to a data-use agreement and IRB. The [NSRR study page](https://sleepdata2.org/datasets/mros/pages) states its MrOS release does not include actigraphy raw data or results; the published study's complete phenotype set is not thereby guaranteed in that release. An authorized linkage remains unverified. MrOS also contributed to the project's previously inspected CHARGE gait meta-analysis, requiring a decision-history and participant-overlap audit.

**Longitudinal option:** A [follow-up study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3449012/) assessed baseline sleep in 2,505 men who were not frail, and frailty again about 3.4 years later. Its five components were weight loss since baseline, low grip strength, low reported energy, slow gait, and low physical activity; it retained the sleep-exam cutpoints at follow-up. Its primary published models include death in the outcome (for example frail **or dead** versus robust/intermediate), so their estimates must not be quoted as incident-frailty-only associations. If the relevant observations are available with approval, prespecify a follow-up endpoint and how death and missing follow-up will be handled, rather than assuming the cross-sectional phenotype or composite is equivalent to incident FFS. The 2,505 figure is not the genotyped or consent-eligible intersection and does not establish score validation. The release-specific availability of follow-up components, and the omission of actigraphy results from the NSRR release, must be resolved.

**Genetics and insomnia measures were linked in published MrOS analyses:** [Foldager et al.](https://pubmed.ncbi.nlm.nih.gov/34170227/) included 1,026 MrOS men from the **second sleep visit** in a multicohort analysis combining genome-wide array data, polysomnography and questionnaires. The [methods](https://jcsm.aasm.org/doi/pdf/10.5664/jcsm.9468) describe an insomnia-severity measure and a separate reported restless-legs question for that visit; the study used an existing insomnia genetic score, not our prototypes. Its 1,026 is not the number who also have the project's frailty definitions, the needed SNPs, or consent for our proposed use. This strengthens feasibility of genotype–sleep linkage but does not establish genotype–sleep–frailty linkage. Sleep visit 2 was 2009–2012, whereas the frailty/sleep study above used the earlier sleep exam; require careful visit timing and participant-level linkage rather than treating these as the same measurement.

That study found the candidate rs113851554 at MEIS1 associated with periodic limb movements. A [separate multicohort GWAS](https://pubmed.ncbi.nlm.nih.gov/35670608/) also reported this locus for periodic limb movement in sleep and included MrOS in discovery. The project 15-marker prototype contains rs113851554, which warrants a prespecified sensitivity to restless-legs/limb-movement phenotype where permitted, and caution when interpreting it as insomnia-specific. These studies do not establish this variant causes frailty or validate our score; do not use the MrOS-containing association as an independent follow-up for the same participants.

**Important release and timing constraint:** The [MrOS analyst guide](https://agingresearchbiobank.nia.nih.gov/studies/mros/documents/download/Visit_Information.pdf/) says the calculated CHS-like frailty variables were created where possible through sleep visit 2 but **are not released** in the dataset it documents; the guide directs investigators to the coordinating center for those calculated variables and their definitions. This statement does not prove that all raw components or every other release are unavailable; verify the current release, permissions and any applicable fees before choosing a reconstruction or request. The guide identifies the Insomnia Severity Index (`SLISISCR`, `SLISICAT`) at sleep visit 2. A [published visit timeline](https://academic.oup.com/sleep/article/45/11/zsac206/6678752) dates sleep visit 2 to November 2009–March 2012, after the main visit 3 window in 2007–2009. Thus the visit-2 insomnia measure cannot be assumed to precede the earlier visit-3 frailty outcome; examine individual dates and any later frailty components before making a prospective claim. No request to the coordinating center has been sent.

**Later-visit reconstruction feasibility:** The same guide lists weight, grip strength, six-meter walk, the PASE activity score (`PASCORE`), and Geriatric Depression Scale data at **visit 4**, as well as at sleep visit 2. A [published visit description](https://www.sciencedirect.com/science/article/abs/pii/S2352721825002232) dates visit 4 to May 2014–May 2016, following sleep visit 2. This makes a later modified physical-frailty outcome plausible. However, the guide's inventory alone does not establish that the exact energy item, intentional weight-loss determination, visit-2 reference cutpoints or all needed person-level values are distributed or linkable to genotype and insomnia data. The guide's non-release warning applies to its calculated frailty variables, and it does not list a calculated visit-4 frailty value. Participant-level dates, missingness, death before follow-up and inability to perform tests require explicit handling. This remains a proposed reconstruction, not an existing validated visit-4 FFS or a completed prediction analysis.

**Visit-4 phenotype confirmed in published research:** A [2024 MrOS study](https://pubmed.ncbi.nlm.nih.gov/38819605/) assessed a five-component physical frailty phenotype at visit 4 in 490 men selected for its muscle and kidney biomarker analysis. A [2026 prospective MrOS study](https://pubmed.ncbi.nlm.nih.gov/42844230/) also reports visit-4 frailty following visit 3 and treats deaths as a separate outcome category. Thus visit-4 phenotypic frailty has in fact been constructed in published work; the previous sentence refers specifically to the *project's proposed reconstruction* and does not suggest such a phenotype has never existed. These published samples are not confirmed to overlap the 1,026-person sleep/genetics study, and neither paper establishes distribution rights, identical cutpoints, marker coverage or a usable joint sample. The analyst guide's non-release statement and the need to inspect current components still apply.

**No-fee dbGaP route and eligibility:** An [NCBI educational guide](https://www.nlm.nih.gov/ncbi/workshops/2023-05_genetic-disorder-mechanisms/supplemental-material/Flyer_Pharmocogenetics_CaseStudy.pdf) describes dbGaP data as provided at no cost, with an application required for participant-level data. MrOS [phs000373](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000373.v1.p1) lists QC-cleaned PLINK data through dbGaP. This is evidence of a no-data-fee **controlled-access** genetic route, not proof that the complete sleep, later frailty and genotype linkage is included, authorized or free of computing, institutional or coordinating-center charges. The [NIH application guidance](https://grants.nih.gov/policy-and-compliance/policy-topics/sharing-policies/accessing-data/dbgap) requires an eligible institutional investigator and institutional signing official; trainees may participate under an eligible investigator but cannot submit the request themselves. Consent-group limitations continue to apply. This clarifies the earlier statement that no-fee genetics access had not been established, without changing the unverified joint-data and project-specific cost status.

## VIKING Genes access check

The [VIKING access page](https://viking.ed.ac.uk/our-data-and-samples/access) offers an openly accessible data dictionary and research-specific summary data, and describes VIKING II/III questionnaire data through DPUK. It says exome and whole-genome sequence data on subsets are in EGA with access managed by a data access committee; clinic exomes are shared under contract. This does not establish open linked individual-level genetics and outcomes for the project. Sleep and frailty definitions, linkage, score-variant coverage, overlap, independence and cost remain unverified. Open documentation and summary results should not be confused with an openly downloadable validation cohort.

## Genes & Health access check

The [Genes & Health resource](https://www.genesandhealth.org/researchers/data/) includes genetic and health data for British Bangladeshi and Pakistani volunteers, making it potentially useful for future ancestry-diverse research. However, its [current application page](https://www.genesandhealth.org/researchers/apply-for-access/) lists a mandatory £5,000 plus VAT institutional setup charge and £1,000 plus VAT per user per year for academic participant-level data-only access, with approved use in its secure environment. Those published terms do not satisfy this project's free-data constraint; actual sleep/frailty phenotype coverage and project-marker availability have not been assessed. The resource also offers public summary results, and a [blood-test GWAS paper](https://www.nature.com/articles/s41467-024-53091-x) supplies downloadable summary statistics; these are molecular/biomarker associations, not linked sleep–frailty score validation. No application, payment or participant-data download was made.

The resource's [data page](https://www.genesandhealth.org/researchers/data/) describes precomputed GWAS results as free under CC BY-SA and points to the public `genesandhealth_publicdatasets` cloud bucket. A direct unauthenticated bucket-list request on 29 October 2026 returned HTTP 403 with an organizational-policy restriction; the provider-recommended `gcloud` tool is not installed in this workspace. We therefore did not establish the bucket's phenotype inventory or whether it includes sleep or frailty results. This response does not prove the advertised results are unavailable to other users. Any accessible summary results would permit association follow-up only, not direct validation of the scores in people.

## Simons Sleep Project access check

A [published sleep study](https://www.nature.com/articles/s41593-025-02146-3) describes wearable sleep data and linked exome sequencing for 200 participants aged 10–17, including sibling pairs. Despite the description as an open-science resource, the paper specifies institutional review board approval and a data-transfer agreement for sleep data, with separate exome access through SFARI. It is not an open linked adult-frailty validation cohort; its age range, relatedness, selected autism/sibling sampling and exome coverage also limit applicability. No fee or project-specific score-variant coverage was verified.

## Sukoyaka Health Survey sleep/genetics lead

The Japanese [Sukoyaka Health Survey dataset JGAD000809](https://humandbs-production.ddbj.nig.ac.jp/en/dataset/JGAD000809) describes adult cohort questionnaires, sleep diaries and sleep heart-rate/brain-wave data together with SNP-array genotypes for 962 participants; a related [phenotype release](https://humandbs-production.ddbj.nig.ac.jp/en/dataset/JGAD000808) is listed separately. Both are **controlled-access Type I**. The [data-use instructions](https://humandbs-production.ddbj.nig.ac.jp/en/data-use) require an account and application; fee, institutional eligibility, overseas use terms, marker coverage, usable genotype–sleep overlap and discovery overlap have not been verified. The listing describes healthy adult volunteers but does not establish an FI, Fried phenotype, age distribution, or frailty case numbers. It is therefore a possible future sleep-genetics or methods lead, not a demonstrated no-fee open cohort or a validation dataset for the sleep–frailty scores. No data or account were obtained.

The provider's [public data summary](https://humandbs-production.ddbj.nig.ac.jp/en/hum0395-data-summary) clarifies that the phenotype files include an **age field**, body measurements and health-history questionnaires, plus a diary containing awakenings and sleep measurements listing duration and number of awakenings. The summary lists 625 measured people in 2019 summer and 319 in 2020 summer, but these counts do not establish usable overlap with genotypes or the number of older adults. It does not explicitly list grip strength, gait speed or a constructed FI/FFS; the short questionnaire description is not a complete item inventory, so their absence from this page is not proof of absence from the files. These details strengthen the sleep-phenotype lead while leaving frailty, data access and score validation unresolved.

## ELSA access wording recheck

The [SSGAC access page](https://www.thessgac.org/data-access-procedures) describes ELSA genotype data as publicly available. The more specific [ELSA genetics page](https://www.elsa-project.ac.uk/genetics) says directly genotyped data are accessed via EGA, while QC and imputed genotypes are available directly from the study team; its precomputed PGS is distributed via the UK Data Service. These are different products and the generic statement does not establish an unrestricted download of individual genotypes. Precomputed scores from other GWAS cannot reconstruct the two project prototypes. The ELSA page reports 7,183 samples and 7,179,780 imputed variants after QC, but project-variant coverage, relevant linked outcomes and discovery/relative overlap still require verification. The current [ELSA genetic-data access policy](https://www.elsa-project.ac.uk/genetic-data-access) explicitly charges all successful applications for genetic data linked to phenotypes on a project-specific cost-recovery basis, so ELSA does not satisfy the requested free linked-data route. No application or data download was made.

## Personal Genome Project Canada status

[PGP Canada](https://personalgenomes.ca/) says its former data-hosting application has been retired and that data will be uploaded to a public repository. The page does not provide a currently usable repository link. Its consent for public sharing is promising but does not establish present access, linked sleep/frailty outcomes or adequate sample size. No repository availability date was verified.

## Historical ELSA DNA Repository candidate panel

A [2015 frailty paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4877432/) calls the ELSA DNA Repository public and reports a 620-SNP candidate panel in 3,160 people, with an adapted Fried phenotype. However, an [ELSA cohort profile](https://pmc.ncbi.nlm.nih.gov/articles/PMC3900867/) states that genetic access required an application; an [older access document](https://www.ucl.ac.uk/population-health-sciences/sites/population_health_sciences/files/5050_elsa_genetics_data_access_v3.pdf) listed fees for linked candidate-genotype and phenotype requests. Those are historical conditions, not verified current prices. The published supplementary table includes 590 analyzed SNPs, none with an exact rsID match to either project score; the remaining panel markers, current availability and linked sleep outcomes were not established. See the [coverage check](stage6_elsa_ednar_590_marker_review.md). This historical panel is therefore a lead, not proof of an unrestricted free validation dataset.

## PAAWS sleep and movement resource

[PAAWS Release 2](https://www.paawsstudy.org/data.html) reports sleep polysomnography, multi-site accelerometry and expert scoring for 156 participants (282 nights), plus a smaller released free-living movement subset. This could support future wearable sleep/activity methods, but the released-data description does not establish linked genome-wide genotypes, insomnia diagnosis, FI/FFS, or sufficient independent participants to validate these genetic scores. The study's repository download link returned HTTP 403 in this environment; the authors note that repository throttling can cause this, so the result does not establish that the data are generally unavailable. License and exact participant linkage still need checking before use. The authors' [benchmark repository](https://github.com/mHealth-Research-Group/paaws-benchmarking) separately exposes activity-model code and example data; those are not a linked genetic-score validation set. No PAAWS participant data were downloaded.

## openSNP availability

The [openSNP maintainers announced](https://opensnp.wordpress.com/2025/03/31/sunsetting-opensnp-a-personal-retrospective/) that the service and stored data would be removed in April 2025. Historical descriptions of openSNP as a public genotype–phenotype resource are therefore not evidence of a currently accessible live cohort. No third-party archive was obtained or analyzed.

## COFITAGE open-data lead

A 2026 preprint describes sleep, neuroimaging and genotyping in 101 adults aged 50–69. Its full text describes 84 with APOE genotyping, but its repository link is a placeholder, and we have not verified public genome-wide variant files, FI/FFS, or score-marker coverage. The two APOE polymorphisms explicitly named in the curated variables occur in neither project score. The small healthy sample and reported sleep-disorder exclusions in a related COFITAGE study limit its potential as a validation cohort. [Feasibility assessment](stage6_cofitage_public_validation_feasibility.md).

## Berlin Aging Study lead

BASE public documentation describes an insomnia complaint and observed signs of need for help; a study page reports DNA analysis for about 380 participants. The available sources do not establish a publicly downloadable linked genome-wide dataset, the project score variants, or FI/FFS. BASE-II is a separate cohort. [Detailed feasibility](stage6_base_sleep_frailty_genetics_feasibility.md).

## Scripps sleep-clinic targeted-panel lead

A published clinic study measured sleep phenotypes and used a 768-marker panel, but its full panel coverage and access to linked participant-level genotypes and frailty outcomes were not established. Its article's 16 named rsIDs do not match either project score; this is not a screen of its complete panel. [Feasibility and source check](stage6_scripps_sleep_genetics_feasibility.md).

## ADNI no-fee application-based lead

ADNI has published frailty-index and sleep studies and provides genotyping to authorized investigators. It requires an application and data-use agreement; linked sample size, project-marker coverage, costs for this use and discovery overlap remain unverified. Three published FI definitions all contain explicitly named sleep or fatigue-related items; a separate published 11-item comorbidity-focused index has no explicitly named sleep or fatigue item but differs from the project’s FI and FFS. [Detailed feasibility](stage6_adni_sleep_frailty_genetics_feasibility.md).

## NACC and NIAGADS lead

NACC describes its requested phenotype data as free, and has published frailty-index and informant-reported sleep analyses. NIAGADS lists linked ADRC array and imputed-genotype data under controlled access; download costs for most files are covered according to its general guidance, but the specific project, permissions and costs remain unverified. A published NACC FI contains nighttime behavior and low-energy items, so overlap with sleep must be addressed. The genotyped intersection and score-marker coverage are unknown. [NACC feasibility](stage6_nacc_sleep_frailty_genetics_feasibility.md).

## FRAILPOL openly licensed frailty data

[FRAILPOL](stage6_frailpol_gait_baseline_review.md) provides 668 older adults' gait features and Fried-phenotype classes under CC BY 4.0, but the downloaded gait file contains no genetic or sleep variables. An internal ML baseline and [exploratory transfer to GSTRIDE](stage6_gstride_transfer_review.md) are useful for physical-frailty methods, not genetic-score validation; the outcome and gait inputs share slowness-related content.

## PEGS genotype and sleep lead

The [Personalized Environment and Genes Study](https://doi.org/10.1038/s41597-026-07011-x) describes sleep surveys and whole-genome data, but states that participant data are not in a public repository. Its [2026 collaboration guidelines](https://www.niehs.nih.gov/research/atniehs/labs/crb/studies/pegs/collaboration/guidelines) require a PEGS investigator, proposal review and agreements, and say reimbursement may be required case by case. Linked frailty outcomes, sample independence and score-marker coverage remain unverified. It is not an unrestricted free download.

## NCDS no-fee application-based lead

The 1958 British birth cohort has a published age-50 FI and sleep questions, and CLS reports 6,396 imputed genotypes available without charge following application. The FI includes sleep and fatigue items; participant linkage, exact marker coverage and discovery independence remain unverified. It is not an unrestricted open download. [Detailed assessment](stage6_ncds_sleep_frailty_genetics_feasibility.md).

## ASPREE/ALSOP controlled-access lead

ASPREE/ALSOP has published sleep and adapted Fried measures in the same older cohort, with separate large-scale genotyping efforts. It is a potentially valuable target, but linking the actual genotype, sleep and frailty files, no-fee eligibility and access conditions remain unverified. It is not a currently usable openly available participant-level dataset. See the [feasibility audit](stage6_aspree_alsop_validation_feasibility.md).

## Open Personal Genome Project lead

The Harvard Personal Genome Project describes public linked genomic and trait data. In an older, accessible snapshot, only nine people have an insomnia-labeled condition, demographic metadata and a listed genetic download of the types checked; usability, current sharing status, score-allele coverage and a frailty outcome are unverified. This does not support a reliable score test or rule out other PGP resources. See [public-source feasibility](stage6_pgp_validation_feasibility.md).

## Open Humans metadata lead

A broad metadata scan had unreliable pagination. A later project-specific scan found 32 owners with both genetic-upload and wearable-source files in the projects checked, with one wearable metadata description mentioning sleep. A limited keyword scan of 11 PGP-source health/trait JSON files found none of the selected sleep or frailty terms; this does not exclude data represented elsewhere or with other wording. Frailty outcomes remain unverified, and these counts are not a complete census or a score-validation cohort. [Metadata feasibility and quality limits](stage6_openhumans_feasibility.md).

## CHARLS access check

The China Health and Retirement Longitudinal Study has been used in public frailty research and is a possible non-European setting, but access to survey data must not be confused with access to genetic data. Its [2026 genetic-data policy](https://charls.charlsdata.com/pages/Data/2026-charls-geneticdata/en.html) states that individual-level human genomic data are controlled-access for qualified research. Its availability therefore does not establish an immediately accessible public genotype–frailty dataset for this project. No application was made, and exact score-variant coverage, phenotype comparability, participant independence and cost have not been established.

## Public VETSA outcome-data lead

The public VETSA concordance lists sleep-onset and grip-strength measures in male veteran twins, but an accessible copy of the exact genetic variants or an evaluable frailty outcome has not been established. It is a lead for further file-level verification, not a validated score cohort. [VETSA feasibility](stage6_vetsa_validation_feasibility.md).

## NHANES access check

The public 2013–2014 survey can nonetheless provide a same-participant **phenotypic** sleep-question/grip-strength comparison; see the [exploratory analysis](stage6_nhanes_sleep_grip_review.md). It has neither a formal insomnia diagnosis nor full FI/FFS and does not test the genetic scores.

[NHANES](https://www.cdc.gov/nchs/nhanes/biospecimen/dna-specimens.html) has publicly released survey data useful for frailty research, but CDC says most genetic datasets are restricted. The limited anonymized genetic datasets require a data-use agreement and cannot be linked to public-use NHANES files. Thus the public survey and these anonymized files do not provide an immediately usable linked genotype–phenotype cohort for evaluating either project score. This check does not establish whether any restricted NHANES data contain the particular variants or outcomes; no application was made.

## Medical Genome Reference Bank access check

The Medical Genome Reference Bank includes older adults with genome sequencing, and a [published analysis](https://www.nature.com/articles/s41467-019-14079-0.pdf) measured grip strength and gait speed. Its [data access policy](https://sgc.garvan.org.au/assets/MGRB_DAP_v_1_0_3.pdf) distinguishes an open portal of aggregated variant frequencies from controlled individual-level genomic data; detailed phenotypes require additional cohort approval. The [cohort design paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336775/) confirms this separation. Thus the public portal cannot establish individual score prediction or provide a free, immediately open linked genotype–FI/FFS validation cohort. Exact candidate coverage, insomnia phenotype, relevant outcome and discovery overlap have not been evaluated; no application was made.

## Earlier controlled-access assessment (not the current plan)

All of Us is a plausible first setting for an independently assessed, EHR-based frailty outcome if the team has Controlled Tier access and can establish that the evaluated participants and close relatives were not included in discovery or candidate selection. BioVU is a possible further evaluation setting for Vanderbilt-affiliated collaborators. Neither cohort is automatically independent of the other or of all discovery datasets; nor do the published resources establish that the original five-component Fried frailty score (FFS) and the Atkins frailty index (FI) can be reconstructed in the same people. Confirm actual overlap under applicable data-use rules before describing a replication as independent.

A recent [preprint](https://www.medrxiv.org/content/10.64898/2026.09.24.26363888v1.full) applied the Hospital Frailty Risk Score (HFRS), an ICD-based measure, in UK Biobank, All of Us, and BioVU. It demonstrates a feasible phenotype in those research settings, not that the 15- or 22-variant germline scores predict frailty. Its genetic exposure is somatic clonal hematopoiesis, a different biological and statistical question. HFRS is not a substitute for FI or FFS and should not be assumed sleep-free.

[All of Us](https://www.researchallofus.org/data-tools/workbench/) makes genomic data available in its Controlled Tier to registered researchers at institutions with the appropriate data-use agreement, and charges apply for computing and storage. Its published access guidance restricts importing identifying information and row-level linkage of external data; a proposed cross-cohort overlap check must therefore be reviewed with the program and each data steward, rather than assumed technically or legally possible. [Vanderbilt](https://victr.vumc.org/data-resource-biovu/) describes BioVU as linked to longitudinal EHR and genotyping or whole-genome data, with access limited to Vanderbilt-affiliated investigators following application. These descriptions do not establish this project's access, eligible sample size, available variants, or the presence of particular measures.

## Test to preregister

1. Freeze the candidate list, alleles, build, and the three trait-specific weight columns before opening outcome results. The [15-marker prototype](unvalidated_15_variant_score_prototype.csv) and [22-marker prototype](unvalidated_22_variant_score_prototype.csv) are experimental, not validated or interchangeable. Do not combine weights from different phenotypes as though on the same scale. Resolve multiallelic encoding and FI cohort sign issues first.
2. Define a baseline date, required prior EHR coverage, ancestry/relatedness handling, genotype quality and missingness. Check exact allele and genomic-build matches, including multiallelic sites; quantify usable markers before calculating a score. Independently assess discovery, selection, and comparator-score sample overlap as permitted; if it cannot be excluded, describe the limitation and avoid an independent-replication claim. If All of Us and BioVU cannot be de-duplicated under their data-use rules, do not treat them as independent replications.
3. Specify an incident frailty outcome and follow-up time using a documented, implementable HFRS algorithm. Distinguish pre-existing frailty from incident events, and ensure all baseline covariates precede the outcome. Investigate EHR availability, healthcare utilization, diagnosis coding, age, and competing risk of death. Report feasibility and count exclusions. Verify coding and whether sleep-related codes enter the score.
4. Prespecify primary model and comparisons: age, sex, ancestry principal components, and relevant assessment variables in a baseline model; add each frozen genetic score separately. Report effect per within-analysis standard deviation with confidence intervals and discrimination/calibration plus incremental performance over baseline. Compare with a published HFRS polygenic score only after verifying weights, target data, allowed use, coverage, ancestry applicability and independence. Multiple candidate scores and outcomes require an explicit multiple-testing plan.
5. Treat measured insomnia association and gene-level or cell-level analyses as separately defined secondary aims. Evaluate actual FI or five-component FFS only if each can be operationalized and validated using available participant measures. Stratified or cross-ancestry results require sufficient sample sizes and should not be interpreted as causal.

The small EUR reference-panel calculation in [the 22-versus-15 diagnostic](reference_22_vs_15_score_diagnostic.md) tests score construction and sensitivity to marker choice, not association with patient outcomes. Without eligible participant-level genotypes and phenotype outcomes, an out-of-sample predictive result cannot be calculated.

## Local NLSY79 lead checked

A separate local sleep-aging project contains an `nlsy79_all_1979-2022.zip` archive. Its directory listing contains survey data and documentation in CSV, DAT, codebook and statistical-program formats. The [NLS cohort download documentation](https://www.nlsinfo.org/content/access-data-investigator/accessing-data-cohorts) describes these as survey-data files, not a genotype release. A filename search for common PLINK, BGEN and compressed VCF extensions in three relevant visible local project trees found no candidate genotype files. This is a limited search: it does not establish that genotypes are unavailable elsewhere, that no genetic variables exist in a survey file, or that any particular phenotype can be derived.

Consequently, the presence of the local survey archive does not establish an executable validation of the 22-marker genetic count. Before using this cohort, establish permitted participant-level genotype access, exact variant coverage and linkage, an appropriate measured frailty and insomnia phenotype, sample size, relatedness and discovery-cohort independence. No individual survey records were extracted or linked in this check.

## PGP-UK open metadata check

The separate [PGP-UK public API metadata screen](stage6_pgpuk_open_cohort_feasibility.md) returned 118 records and VCF-labelled variation metadata for 100, but its limited questionnaire search did not establish insomnia and FI/FFS outcomes. These are metadata counts, not verified linked genotypes or a representative cohort; no score validation was possible. The screen is distinct from the Harvard PGP lead above.

## ELSA linked-genetic access clarification

**Distinction:** Precomputed PGSs (including insomnia) are a separate [safeguarded UK Data Service product](stage6_elsa_precomputed_scores_access_review.md), with registration and an end-user licence. The fee and approval rules for custom linked SNP data do not by themselves settle whether linkage of this PGS product to standard survey data is free or allowed. Its existing scores also cannot compute this project’s 15- or 22-marker scores.

The [ELSA genetics page](https://www.elsa-project.ac.uk/genetics) describes individual direct genotypes through EGA, quality-controlled/imputed data through the ELSA team, and a separate process for linkage to phenotypes. The [EGA dataset listing](https://www.ega-archive.org/dacs/EGAC00001000270) is under a data-access committee, and ELSA's [FAQ](https://www.elsa-project.ac.uk/faqs) states that a fee applies to genetic data linked to phenotypic data. A third-party statement that genetic resources are publicly available does not establish unrestricted, no-cost linked individual data. Fees and procedures may change; no application or request was made. This is outside the user's public-free data constraint as currently documented.

## Human Phenotype Project lead

The [public catalog](https://data-browser.humanphenotypeproject.org/) lists human genetics, sleep monitoring and hand grip, but a [sleep study's data statement](https://www.nature.com/articles/s41591-024-03481-x) says the individual data cannot be publicly released and directs eligible researchers to a managed-access process. Its [access page](https://humanphenotypeproject.org/data-access) describes a trusted research environment. Current no-cost eligibility, matching genotype/phenotype coverage, frailty definitions and independence are not established. The catalog alone is not public individual-level data and no application or external message was sent.

## All of Us public association browser

The [All by All public-browser screen](stage6_aou_public_variant_screen.md) found insomnia and generalized-weakness phenotype analyses and tested 22 exact candidate mappings in a selected-association endpoint. No lookup returned an association row; positive controls succeeded. The endpoint is restricted to its selected association table, so blank results cannot be treated as null or nonsignificant associations and do not provide the summary data required for score analysis. The browser is distinct from controlled-access participant-level data.

## China Kadoorie Biobank public browser

A [phenotype-code inventory screen](stage6_ckb_catalog_screen.md) did not identify a directly labelled insomnia, FI or FFS GWAS in the CKB PheWeb catalog, but codes are abbreviated and the search cannot establish absence of relevant measures or other CKB studies. The site says its hosted summary files are encrypted and keys require an application; some separate CKB summaries are publicly hosted elsewhere. Browser and download genome builds differ for some studies. No participant data or relevant full summary GWAS was obtained through this check.

## Canadian public summary-data lead

A [CLSA and CARTaGENE feasibility review](stage6_canadian_public_gwas_feasibility.md) identifies a 2026 Canadian ordinal frailty GWAS and public PheWeb summary portals. The portals reported maintenance until November 1; the exact GWAS, relevant phenotypes, usable variant data and independence remain unverified. The publicly described summary portal does not establish open linked individual-level data.

## Additional NHANES survey cycle

The separately sampled public 2011–2012 cycle contains the same trouble-sleeping question and grip measure as 2013–2014. The [fixed-method follow-up](stage6_nhanes_2011_sleep_grip_followup.md) produced an imprecise adjusted difference of +0.19 kg (95% CI −1.43 to 1.80). It is an additional phenotype-level check, not FI/FFS, clinical insomnia, or genetic score validation; sequence-ID separation alone does not prove no participant overlap.

## TOPMed handgrip study

A [TOPMed grip GWAS feasibility check](stage6_topmed_grip_validation_feasibility.md) found public selected supplementary tables but the paper describes full summary results as controlled access. No current 22-candidate rsID was found in the selected tables; missing entries are not null results. Grip is a frailty component, and this source does not yet provide an independent genetic-score validation.

## African-ancestry grip summary GWAS

The newly available Hansen 2026 public handgrip GWAS permitted a [22-candidate lookup](stage6_hansen_grip_candidate_followup.md): 16 exact position-and-allele matches, with one nominal association that did not pass a 22-test Bonferroni threshold. It offers cross-population component-phenotype evidence, not a participant-level score validation or a full FI/FFS analysis; source effect units remain unresolved.

## Frailty ontology inventory

A [current GWAS Catalog frailty-term inventory](stage6_frailty_efo_catalog_inventory.md) found 11 curated study entries, including nine with Catalog full-summary flags; those nine correspond to sources already examined. This bounded query does not rule out studies under different phenotype terms or files hosted outside the Catalog.

## CHARGE gait-speed study

A [public selected-variant screen](stage6_charge_gait_candidate_feasibility.md) checked an older-adult, multi-cohort gait-speed GWAS. None of the 22 candidate rsIDs appeared in its 536 selected supplementary rows; absence from that list is not a null association. Catalog metadata does not identify full summaries as available there, and participant-level score validation remains unavailable.

**Update to CHARGE gait-speed lead:** the [full GRASP results and follow-up](stage6_charge_gait_full_candidate_review.md) are now available and were scanned; 11 of 22 candidate rsIDs and allele pairs matched, with no nominal p<0.05. This resolves the full-file availability question but does not provide linked genotypes and FI/FFS outcomes for score validation.

## NCI PLCO summary-data atlas

A bounded [PLCO atlas phenotype-label screen](stage6_nci_plco_atlas_feasibility.md) returned no results for selected sleep, insomnia, frailty or fatigue search terms. This is not a comprehensive survey of the underlying cohort, and it did not provide score validation data.

### Additional search: autoMR example files (2026-10-28)

The public autoMR package contains 10-instrument FI and Fried-frailty association examples, not participant-level genotypes with frailty outcomes. None of its FI instrument rsIDs exactly matches either current score candidate list; source-cohort independence is not established. It is therefore not a score-validation cohort. Details and reproducibility receipt: [autoMR feasibility audit](stage6_automr_frailty_example_feasibility.md).

### Additional search: public CHARLS sleep–psychological-frailty release (2026-10-28)

A [public processed Sleep-PFI CHARLS dataset](https://data.mendeley.com/datasets/76zz3d74r6/2) is described as containing sleep duration and psychological frailty, and the checksum-verified main R data file contains no genotype or genetic-score columns; linkage permissions remain unverified. It is a phenotypic background and reproducibility lead, not a validated genetic-score dataset. See [source audit](stage6_charls_sleep_pfi_open_data_lead.md). ELSA's current [genetic data access page](https://www.elsa-project.ac.uk/genetic-data-access) independently confirms cost-recovery charges for approved genetic–phenotype linkage.

### Additional search: Stanford ADRC SPHERE synthetic release (2026-10-28)

The public [SPHERE portal](https://adrc-sphere.stanford.edu/) describes multimodal synthetic records with selected AD-associated variants. Its public variant dictionary has no exact rsID overlap with either project prototype, and reports carrier coding rather than additive 0/1/2 dosage. Synthetic records do not constitute independent patients. Registration and a research-use agreement apply to the data. See [feasibility audit](stage6_stanford_sphere_validation_feasibility.md).

### Additional search: PISA (2026-10-28)

PISA has published linked genotype–FI analysis and an onsite sleep-monitoring subset, but individual data require a project-specific request and data-sharing agreement; cost, permitted linkage, sample overlap and marker coverage are unverified. Its published FI explicitly contains multiple sleep-related deficits. It is a relevant controlled-access lead, not a confirmed free public validation set. See [PISA feasibility and FI item audit](stage6_pisa_sleep_frailty_validation_feasibility.md).

### Additional search: South Texas Aging Registry and Repository (2026-10-28)

[STARR](https://sapepper.barshop.uthscsa.edu/cores/starr/) lists performance-based frailty measures, DNA specimens and genetic data generated in subsequent research. This does not establish an existing cohort with linked genome-wide genotypes, insomnia, FI/FFS and adequate marker coverage. Requests require approval and an agreement; its policy assigns distribution costs to the recipient. It is therefore not an established no-fee, public score-validation set. No request was made.

### Additional search: HANDLS public PBMC frailty expression (2026-10-28)

[GSE129534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129534) provides 16 public frailty-labelled PBMC expression files. They contain expression counts, not the linked genotype dosages and insomnia outcomes required for score validation. A cautious four-gene detection screen is in the [HANDLS review](stage6_handls_candidate_expression_review.md).

### Additional search: ROS/MAP (2026-10-28)

ROS/MAP is a substantive potential validation lead: separate publications describe frailty-index construction and genome-wide genotyping, and RADC lists sleep questionnaires and actigraphy. RADC’s aggregate availability tool reports 401 MAP participants with both a sleep-questionnaire flag and at least one listed genotyping/WGS flag; usable individual insomnia and frailty outcomes and exact marker coverage remain unverified. Phenotype access requires an approved [RADC request](https://www.radc.rush.edu/requests/data.htm) and institutional agreement; Synapse genomic assets have their own controls. No specific linked-data cost was established. See [ROS/MAP feasibility](stage6_rosmap_sleep_frailty_genetics_feasibility.md).

### Additional search: CIMA-Q (2026-10-28)

The [CIMA-Q catalog](https://portal.conp.ca/dataset?id=projects/CIMA-Q) lists insomnia severity and a Fried frailty measure, but describes restricted consortium access and does not establish available genome-wide genotype data for our markers. Cost and the jointly usable sample are unknown. It is a phenotype lead, not an immediately usable public score-validation dataset. [Details](stage6_cimaq_sleep_frailty_feasibility.md).

### Additional search: NSHAP (2026-10-28)

NSHAP has a published actigraphy–physical-frailty analysis and a separate publication describing genome-wide genotyping and imputation; usable genotype–sleep–frailty overlap, candidate-marker coverage, access and costs are unknown. Its [AI policy](https://www.norc.org/content/dam/norc-org/pdf2026/NSHAP-AI-policy%20statement2.pdf) restricts upload of NSHAP data to AI platforms. No individual data were accessed. See [NSHAP feasibility audit](stage6_nshap_sleep_frailty_genetics_feasibility.md).

### Additional search: LASA (2026-10-28)

LASA has published work using both FI and physical frailty, and its cohort profile describes 2,279 unrelated European-ancestry genotyped participants with imputation. A sleep-calendar subset is described, but a usable insomnia–frailty–genotype intersection, independence, candidate variant coverage and cost have not been established. Access is by proposal and agreement. See [LASA feasibility review](stage6_lasa_sleep_frailty_genetics_feasibility.md).

### Additional search: PPMI (2026-10-29)

PPMI has published insomnia-symptom and mobility measures and an authorized genomic-data route. Genetic data require a request and agreement; current cost, linkage, variant coverage and an FI/FFS outcome remain unverified. Its Parkinson-focused sampling makes it a secondary disease-specific lead, not an established public score-validation cohort. See the [PPMI feasibility review](stage6_ppmi_sleep_genetics_validation_feasibility.md).

The [expanded GEPITEL array/reference search](stage6_gepitel_public_validation_feasibility.md) found nearby reference tags with r² ≥ 0.8 for 6/15 and 6/22 score candidates under unverified platform-build assumptions. It does not establish target imputation quality or a linkage to frailty outcomes.

### LBC1936 wave-aligned candidate (2026-10-29)

LBC1936 has published Wave 3 Fried-type frailty and public documentation of Wave 3 sleep, grip and walking variables, plus genome-wide genotyping. Linked data require a request; fees, usable intersection, variant coverage and discovery/selection overlap remain unverified. See [LBC1936 feasibility](lbc1936_sleep_frailty_genetics_validation_lead.md).

A further [E-Risk review](stage6_erisk_sleep_aging_feasibility.md) identifies free managed access for eligible researchers and an observed association of sleep quality with a methylation aging measure at age 18. Genotype linkage and our exact markers are unverified, and aging biomarkers in young adults are not FI or Fried frailty; it is a lower-priority exploratory lead rather than a score-validation cohort.

The [VIKING II dictionary audit](stage6_viking_ii_phenotype_dictionary_review.md) now identifies sleep-onset latency, tiredness and activity fields in the public documentation. It finds no explicitly labelled grip field in the scanned workbook and does not establish an FI or FFS; mood-related sleep questions have unverified routing. Individual genetics, linkage, fees and joint coverage remain unverified.

[COSFAMM circadian study](stage6_cosfamm_clock_sleep_frailty_review.md) measured insomnia and Fried frailty in older Mexican adults. Its underlying data are available on request rather than public; genome-wide genotype availability, linkage, project eligibility and cost were not established. Its published clock-gene expression results are hypothesis-generating and cannot validate the genetic scores.

[Ansan/Ansung KoGES linkage lead](stage6_koges_ansan_ansung_sppb_linkage_lead.md): insomnia/genetics and SPPB/genetics have been separately reported in the same named subcohorts, with controlled access and participant-level linkage still unverified. The reported physical-function case count is small, and SPPB is not FI or FFS.

An actual [PREDICT-NEURODEGEN archive audit](stage6_predict_neurodegen_open_data_feasibility.md) found freely downloadable smartwatch sleep/activity and EEG sheets, but no named genotype, frailty or diagnosis columns in the released workbook headers. It cannot currently evaluate the genetic scores; the broader project description alone does not establish availability or linkage of those other data.

[Mexico City Prospective Study](stage6_mcps_sleep_frailty_genetic_feasibility.md): genotype and sleep-duration data are documented, and a frailty-index project is listed, but their individual linkage and the FI definition remain unverified. The current access policy describes a no-data-charge route for eligible applicants from low- and middle-income countries; high-income applicants face application/usage charges. No fee exemption or eligibility for this project is established, and sleep duration is not insomnia.

## OPPERA sleep-genetics lead

A [sleep-quality GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7953222/) reports genotyped participants from the OPPERA pain cohort and follow-up sleep cohorts. Its discovery outcome is Pittsburgh Sleep Quality Index, not a diagnosis of insomnia; the article explicitly says insomnia was not assessed in the OPPERA discovery cohort. The separately described OPPERA study enrolled adults aged 18–44 ([study catalog](https://datacatalog.hshsl.umaryland.edu/dataset/62)). Its [GWAS dbGaP record](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000796.v1.p1) describes controlled-access genetic data, not a freely downloadable linked genotype–sleep–frailty cohort. It is therefore a sleep-biology or genetic-association lead, not an established free aging-cohort score-validation dataset. The paper also describes a 757-person genotyped insomnia-questionnaire sample in a different replication cohort; do not confuse that sample with OPPERA discovery. Exact candidate-variant coverage and access conditions for those replication studies remain unchecked.

## MIVAC sleep and physical-function lead

The [MIVAC review](stage6_mivac_sleep_grip_genetics_feasibility.md) identifies a longitudinal study with sleep, grip and walking measurements, and genetic biospecimens sent to the NIMH Genetics Repository. No linked genotype release, access fee or candidate coverage has been verified, and detailed phenotypes require approval and a data-use agreement. It is not yet an established free score-validation dataset.

## Brain Resilience Study prospective release

The [Brain Resilience Study](stage6_brain_resilience_study_validation_lead.md) plans an open resource with genotyping, sleep questionnaires and actigraphy in older adults. A linked individual-level data release, complete FI/FFS and actual marker coverage have not been confirmed; its BC Generations parent has a separate controlled application process. It is not yet an executable free validation cohort.

## Korean Welfare Genome Project access check

The [Welfare Genome Project](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.633731/full) reports sequencing and health examinations for 1,000 healthy Korean volunteers. Despite the word “free” in its title referring to participant health checks and genetic reports, the article's data-availability statement says its datasets are not readily available and directs researchers to request them from Ulsan University Hospital Biobank. No anonymous individual-level download, no-cost researcher access, exact candidate coverage, linked sleep and FI/FFS definitions, or adequate sample for this task has been established. No request was submitted. This is not currently a usable free public validation cohort.

## SABE Brazil

The [SABE assessment](stage6_sabe_brazil_genetics_sleep_frailty_feasibility.md) identifies an older admixed Brazilian whole-genome cohort and separate descriptions of sleep and frailty measurements. Individual-level sequencing and phenotypes are controlled access, and no-cost eligibility, linked sample and candidate variant coverage are unverified; aggregated ABraOM frequencies cannot test scores.

## Semmelweis Study lead

The [Semmelweis cohort protocol](stage6_semmelweis_validation_feasibility.md) plans sleep-quality assessment, grip strength and genetic work, but no presently accessible linked genotype–sleep–frailty release or no-fee access was verified. Grip strength and PSQI do not by themselves establish a full frailty or insomnia phenotype. This is a future lead rather than a validation dataset.

## Adult Changes in Thought lead

The [ACT study](stage6_act_sleep_genetics_validation_lead.md) documents older participants, sleep questionnaires and device measures, and controlled-access genetic datasets. A historical dbGaP description names 2,613 ACT participants with viable genotyping specimens, not a confirmed linked sleep–frailty–genotype sample. A published ACT analysis constructed a modified five-component Fried phenotype; its public dictionary also lists an Insomnia Severity Index potentially available since 2022 and a longitudinal CES-D restless-sleep item; the linked genotyped sample, current access and exact marker coverage and project-specific costs remain unverified.

NIAGADS lists 1,397 ACT participants in its ADSP cohort inventory and states that it covers download costs for most joint-called files. ACT data access requires proposal review and external collaboration; linkage, variant coverage and total project cost remain unverified. [ACT route](stage6_act_sleep_genetics_validation_lead.md).

## Baltimore Longitudinal Study of Aging lead

A published [BLSA sleep-genetics analysis](stage6_blsa_sleep_frailty_genetics_feasibility.md) measured an existing sleep score and insomnia symptoms in a selected older sample. A separate BLSA study reports FI, but the joint genotype–insomnia–frailty sample, variant quality, linkage and cost remain unverified; participant data require approval.

## Understanding Society (UKHLS) lead

The [UKHLS review](stage6_ukhls_sleep_frailty_genetics_feasibility.md) documents a 36-item FI, sleep questions, a published joint sleep-duration/genetics analysis and a 9,920-sample controlled EGA genotype release. The FI includes energy but no explicitly named sleep item in its published list. The relevant linked age-50-plus complete sample, variants, independence, eligibility and total cost are unverified; no score validation has been run.

## Cleveland Family Study lead

The [CFS review](stage6_cfs_insomnia_genetics_feasibility.md) identifies TOPMed whole-genome data and a publicly documented insomnia field, but the selected clinical cohort has only 35 positive responses before any genotype or frailty restriction. Genomic access is controlled, a full frailty phenotype and joint sample are unverified, and family relationships and sleep-apnea ascertainment limit generalizability.

## Bogalusa middle-age physical-function lead

BHS has a published sleep–genotype subset and separately published grip and walking measurements from the 2013–2016 visit. Linkage, complete frailty definitions and a free access route are unverified, and the sleep-genetic subset is middle-aged. [Feasibility assessment](stage6_bogalusa_sleep_function_genetics_feasibility.md).

## CARDIA year-35 lead

CARDIA has documented year-35 insomnia-symptom and physical-function assessments, plus controlled-access genetic data. BioLINCC says its datasets carry no access cost and offers exam data through year 35, but sleep/function ancillary availability, genetics linkage, project variant coverage, total costs and FI/FFS construction remain unverified. [Detailed assessment](stage6_cardia_sleep_function_genetics_feasibility.md).
